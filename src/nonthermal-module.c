/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include <slang.h>

#include "_nonthermal.h"
#include "interp.h"
#include "version.h"

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390  /* 2/sqrt(pi) */
#endif

typedef struct
{
   Particle_Type particle;
   double kT;              /* thermal e- temperature [keV] */
   double n_th;            /* thermal e- density */
   double n_GeV;           /* non-thermal e- density @ 1 GeV */
}
Density_Info;
#define NULL_DENSITY_INFO {NULL_PARTICLE_TYPE,0.0,0,0}
/* Note that fitted norm = n_GeV * V / (4*pi*D^2) */

static int pop_pdf_info (char **name, SLang_Array_Type **pars) /*{{{*/
{
   *pars = NULL;
   *name = NULL;

   if ((-1 == SLang_pop_array_of_type (pars, SLANG_DOUBLE_TYPE))
       || (*pars == NULL))
     return -1;

   if (-1 == SLpop_string (name))
     {
        SLang_free_array (*pars);
        *pars = NULL;
        return -1;
     }

   return 0;
}

/*}}}*/

int init_pdf (Particle_Type *p, unsigned int type) /*{{{*/
{
   static char msg[] = "*** nonthermal: init_pdf: failed initializing particle distribution function\n";
   SLang_Array_Type *sl_pars;
   char *pdf_name;
   int status;

   if (-1 == pop_pdf_info (&pdf_name, &sl_pars))
     {
        fprintf (stderr, msg);
        return -1;
     }

   status = init_particle_spectrum (p, type, pdf_name, sl_pars);
   if (status) fprintf (stderr, msg);

   SLang_free_array (sl_pars);
   SLfree(pdf_name);

   return status;
}

/*}}}*/

static int pop_density_info (Density_Info *di) /*{{{*/
{
   Particle_Type *p = &di->particle;
   int particle_type;

   if (-1 == SLang_pop_double (&di->n_th)
       || -1 == SLang_pop_double (&di->kT)
       || -1 == SLang_pop_double (&di->n_GeV)
       || -1 == SLang_pop_integer (&particle_type))
     {
        return -1;
     }

   if (-1 == init_pdf (p, particle_type))
     {
        free_particle_spectrum (p);
        return -1;
     }

   return 0;
}
/*}}}*/

static double thermal_distrib (double *mv, double *pkT_keV, double *pmass) /*{{{*/
{
   /* Maxwellian momentum distribution */
   double p = *mv;
   double m = *pmass;
   double kT = *pkT_keV * KEV;
   double a, x, f;

   if (kT <= 0.0 || p <= 0.0 || m <= 0.0)
     return 0.0;

   /* p is non-relativistic momentum, p = m*v */
   a = 2 * m * kT;
   x = (p*p) / a;

   if (x < 0.0 || 500.0 < x)
     return 0.0;

   f = 2 * M_2_SQRTPI * x * exp (-x) / sqrt (a);

   return f;
}

/*}}}*/

static double particle_distrib (double *pc, unsigned int *type) /*{{{*/
{
   Particle_Type pt;
   double f;

   if (-1 == init_pdf (&pt, *type))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        free_particle_spectrum (&pt);
        return -1.0;
     }

   /* dn/d(Pc) */
   (void) (*pt.spectrum)(&pt, *pc, &f);

   free_particle_spectrum (&pt);

   return f;
}
/*}}}*/

static double root_func (double mv, void *cd) /*{{{*/
{
   Density_Info *di = (Density_Info *)cd;
   Particle_Type *pt = &di->particle;
   double f_th, f_nth, gamma, gamma2, beta, pc;

   if (di->n_th <= 0.0 || di->n_GeV <= 0.0)
     {
        fprintf (stderr, "invalid input:  n_th = %e  n_GeV = %e\n",
                 di->n_th, di->n_GeV);
        /* SLang_set_error (SL_INTRINSIC_ERROR); */
        return 0.0;
     }

   /* non-relativistic momentum, p = m*v */
   beta = mv / (pt->mass * GSL_CONST_CGSM_SPEED_OF_LIGHT);
   if (beta >= 1.0)
     {
        fprintf (stderr, "crazy momentum value:  beta = %e\n", beta);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return 0.0;
     }

   /* dn/dp */
   f_th = thermal_distrib (&mv, &di->kT, &pt->mass);
   f_th *= di->n_th;

   gamma2 = 1.0 /((1.0 + beta) * (1.0 - beta));
   gamma = sqrt(gamma2);
   pc = gamma * mv * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   /* dn/d(Pc) */
   (void) (*pt->spectrum) (pt, pc, &f_nth);
   f_nth *= di->n_GeV * GSL_CONST_CGSM_SPEED_OF_LIGHT / GEV;
   /* dn/dp = dP/dp * dn/dP;   p = m*v, P = gamma*mv */
   f_nth *= gamma2 * gamma;

   return f_th - f_nth;
}

/*}}}*/

static int find_momentum_min (Density_Info *di, double *momentum) /*{{{*/
{
   Particle_Type *pt = &di->particle;
   double mc, p_th, pmax, p, f_lo, f_hi;
   int status = 0;

   /* FIXME: This won't handle pathological cases, e.g. where
    * the two curves cross at a tangent point or where there are
    * two roots above the thermal peak.  I'm not worried about
    * that right now though. */

   p_th = sqrt (2 * pt->mass * di->kT * KEV);
   mc = pt->mass * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   f_lo = root_func (p_th, di);
   pmax = p_th;

   do
     {
        pmax *= 2;
        if (pmax >= mc)
          break;
        f_hi = root_func (pmax, di);
     }
   while (f_lo * f_hi >= 0);

   if (pmax < mc)
     {
        /* root should be bracketed */
        status = bisection (&root_func, p_th, pmax, di, &p);
     }

   /* No solution, use thermal peak */
   if (status || p > mc)
     p = p_th;

   *momentum = p;

   return 0;
}

/*}}}*/

static double _find_momentum_min (void) /*{{{*/
{
   Density_Info di = NULL_DENSITY_INFO;
   double p;

   if (-1 == pop_density_info (&di))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return 0.0;
     }

   if (-1 == find_momentum_min (&di, &p))
     {
        fprintf (stderr, "find_momentum_min: no solution -- using thermal peak momentum\n");
     }

   free_particle_spectrum (&di.particle);

   return p;
}

/*}}}*/

static int Failed_Finding_Momentum_Min;

static int eval_nontherm_integral (gsl_function *f, double *value) /*{{{*/
{
   Density_Info *di = (Density_Info *)f->params;
   Particle_Type *pt = &di->particle;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   double epsabs, epsrel, abserr;
   double p_min, beta, gamma, pc_min, pc_max;
   size_t limit;
   int status;

   *value = 0.0;

   epsabs = 0.0;
   epsrel = 1.e-7;
   limit = MAX_QAG_SUBINTERVALS;

   Failed_Finding_Momentum_Min = 0;
   if (-1 == find_momentum_min (di, &p_min))
     {
        Failed_Finding_Momentum_Min = 1;
     }

   beta = (p_min / pt->mass) / GSL_CONST_CGSM_SPEED_OF_LIGHT;
   gamma = 1.0 / sqrt ((1.0 + beta)*(1.0 - beta));
   pc_min = gamma * p_min * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   /* For simpler continuity into the relativistic regime,
    * we use the relativistic momentum,
    *    Pc = gamma * pc
    * where p = mv.
    */
   pc_max = (*pt->momentum_max)(pt);

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();
   status = gsl_integration_qag (f, log(pc_min), log(pc_max), epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS31,
                                 work, value, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     {
        if (status != GSL_EROUND)
          fprintf (stderr, "status = %d:  %s\n",
                   status,
                   gsl_strerror (status));
     }

   *value /= GEV;

   return 0;
}

/*}}}*/

static int nontherm_integral (Density_Info *di, double (*func)(double, void *), double *x) /*{{{*/
{
   gsl_function f;

   *x = 0.0;

   f.function = func;
   f.params = di;

   if (-1 == eval_nontherm_integral (&f, x))
     return -1;

   return 0;
}

/*}}}*/

static double nontherm_energy_density_integrand (double pc, void *cd) /*{{{*/
{
   Density_Info *di = (Density_Info *)cd;
   Particle_Type *pt = &di->particle;
   double mc2 = pt->mass * C_SQUARED;
   double ne, x, gamma, e_kinetic;

   /* changed integration variable */
   pc = exp(pc);

   /* dn/d(Pc) */
   (void) (*pt->spectrum)(pt, pc, &ne);

   /* relativistic momentum = Pc = gamma * m * v */
   x = pc / mc2;
   gamma = sqrt (1.0 + x * x);
   e_kinetic = (gamma - 1.0) * mc2;

   /* changed integration variable */
   ne *= pc;

   return ne * e_kinetic;
}

/*}}}*/

static double nontherm_density_integrand (double pc, void *cd) /*{{{*/
{
   Density_Info *di = (Density_Info *)cd;
   Particle_Type *pt = &di->particle;
   double ne;

   /* changed integration variable */
   pc = exp(pc);

   /* relativistic momentum = Pc = gamma * m * v */

   /* dn/d(Pc) */
   (void) (*pt->spectrum)(pt, pc, &ne);

   /* changed integration variable */
   ne *= pc;

   return ne;
}

/*}}}*/

static double nontherm_energy_density (void) /*{{{*/
{
   Density_Info di = NULL_DENSITY_INFO;
   double ne;

   if (-1 == pop_density_info (&di))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   if (-1 == nontherm_integral (&di, &nontherm_energy_density_integrand, &ne))
     ne = 0.0;

   free_particle_spectrum (&di.particle);

   return ne;
}
/*}}}*/

static double nontherm_density (void) /*{{{*/
{
   Density_Info di = NULL_DENSITY_INFO;
   double n;

   if (-1 == pop_density_info (&di))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   if (-1 == nontherm_integral (&di, &nontherm_density_integrand, &n))
     n = 0.0;

   free_particle_spectrum (&di.particle);

   return n;
}
/*}}}*/

static double pc_inj (Particle_Type *pt) /*{{{*/
{
   double T_inj, gamma, mc2, pc;

   /* Injection kinetic energy.
    * Presumably non-relativistic for both e- and protons.
    */
   T_inj = 50.0 * KEV;

   mc2 = pt->mass * C_SQUARED;
   gamma = T_inj/ mc2 + 1.0;
   pc = gamma * mc2 * sqrt ((1.0 + 1.0/gamma)*(1.0 - 1.0/gamma));

   return pc;
}

/*}}}*/

static double CR_Electron_Density;

static double charge_error (double p_norm, void *cl) /*{{{*/
{
   Density_Info *dp = (Density_Info *)cl;
   double fp;

   dp->n_GeV = p_norm;
   if (-1 == nontherm_integral (dp, &nontherm_density_integrand, &fp))
     {
        fprintf (stderr, "failed integrating proton distribution\n");
        /* SLang_set_error (SL_INTRINSIC_ERROR); */
        return -1.0;
     }

   return 1.0 - fp * p_norm / CR_Electron_Density;
}

/*}}}*/

static double equal_integrated_nonthermal_densities (Density_Info *de, Density_Info *dp) /*{{{*/
{
   double fe, p_norm;

   if (-1 == nontherm_integral (de, &nontherm_density_integrand, &fe))
     {
        fprintf (stderr, "failed integrating electron distribution\n");
        /* SLang_set_error (SL_INTRINSIC_ERROR); */
        return -1.0;
     }

   CR_Electron_Density = de->n_GeV * fe;

   /* Plausible first guess */
   dp->n_GeV = 100.0 * de->n_GeV;

   /* The CR proton density @1 GeV should be larger than the
    * CR electron density @1 GeV (de.n_GeV) but definitely
    * smaller than the total thermal electron density (de.n_th).
    */
   if (bisection (&charge_error, de->n_GeV, de->n_th, dp, &p_norm) < 0)
     {
        fprintf (stderr, "Error:  couldn't find c.r. proton norm: bisection failed\n");
        /* SLang_set_error (SL_INTRINSIC_ERROR); */
        return -1.0;
     }

   return p_norm;
}

/*}}}*/

static double equal_injection_densities (Density_Info *de, Density_Info *dp) /*{{{*/
{
   double mr = sqrt (GSL_CONST_CGSM_MASS_ELECTRON / GSL_CONST_CGSM_MASS_PROTON);
   Particle_Type *e = &de->particle;
   Particle_Type *p = &dp->particle;
   double fe, fp, p_norm;

   /* Assume equal numbers of protons and electrons are
    * injected with kinetic energy T_inj
    *
    *                      n(P_electron) dP_electron
    * ratio (T_inj) \equiv -------------------------
    *                        n(P_proton) dP_proton
    *
    *  dP_electron/dP_proton = sqrt (m_e/m_p)
    *
    * See Bell (1978) paper II
    */

   (void) (*e->spectrum)(e, pc_inj(e), &fe);
   (void) (*p->spectrum)(p, pc_inj(p), &fp);

   p_norm = mr * de->n_GeV * fe / fp;

   return p_norm;
}

/*}}}*/

static double conserve_charge (int *method) /*{{{*/
{
   Density_Info de=NULL_DENSITY_INFO, dp=NULL_DENSITY_INFO;
   double p_norm;

   if (-1 == pop_density_info (&dp)
       || -1 == pop_density_info (&de))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1.0;
     }

   switch (*method)
     {
      case 0:
        /* drop */
      default:
        p_norm = equal_injection_densities (&de, &dp);
        break;

      case 1:
        p_norm = equal_integrated_nonthermal_densities (&de, &dp);
        break;
     }

   free_particle_spectrum (&dp.particle);
   free_particle_spectrum (&de.particle);

   return p_norm;
}

/*}}}*/

/* DUMMY_CLASS_TYPE is a temporary hack that will be modified to the true
  * id once the interpreter provides it when the class is registered.  See below
  * for details.  The reason for this is simple: for a module, the type-id
  * must be assigned dynamically.
  */
#define DUMMY_CLASS_TYPE   255
#define MT DUMMY_CLASS_TYPE
#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define UI SLANG_UINT_TYPE
#define V SLANG_VOID_TYPE
#define S SLANG_STRING_TYPE

static SLang_Intrin_Fun_Type Intrinsics [] =
{
   MAKE_INTRINSIC_3("thermal_distrib", thermal_distrib, D, D, D, D),
   MAKE_INTRINSIC_2("particle_distrib", particle_distrib, D, D, UI),
   MAKE_INTRINSIC_1("conserve_charge", conserve_charge, D, I),
   MAKE_INTRINSIC("_find_momentum_min", _find_momentum_min, D, 0),
   MAKE_INTRINSIC("_nontherm_density", nontherm_density, D, 0),
   MAKE_INTRINSIC("_nontherm_energy_density", nontherm_energy_density, D, 0),
   MAKE_INTRINSIC_2("bspline_open_intrin", bspline_open_intrin, V, I, I),
   MAKE_INTRINSIC_3("bspline_eval_intrin", bspline_eval_intrin, D, MT, D, D),
   MAKE_INTRINSIC_1("bspline_init_info_intrin", bspline_init_info_intrin, V, MT),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef I
#undef UI
#undef V
#undef S

static SLang_IConstant_Type Intrin_Const [] =
{
   MAKE_ICONSTANT("PROTON", PROTON),
   MAKE_ICONSTANT("ELECTRON", ELECTRON),
   MAKE_ICONSTANT("_nonthermal_module_version", MODULE_VERSION_NUMBER),
   SLANG_END_ICONST_TABLE
};

static char *Module_Install_Prefix = INSTALL_PREFIX ;
extern double Min_Curvature_Pc;

static SLang_Intrin_Var_Type Intrin_Variables [] =
{
   MAKE_VARIABLE("_nonthermal_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("_nonthermal_install_prefix", &Module_Install_Prefix, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("_nonthermal_curvature_momentum_gev", &Min_Curvature_Pc, SLANG_DOUBLE_TYPE, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static void patchup_intrinsic_table (unsigned int valid_id) /*{{{*/
{
   SLang_Intrin_Fun_Type *f;
   unsigned int dummy_id = DUMMY_CLASS_TYPE;

   f = Intrinsics;
   while (f->name != NULL)
     {
        unsigned int i, nargs;
        SLtype *args;

        nargs = f->num_args;
        args = f->arg_types;
        for (i = 0; i < nargs; i++)
          {
             if (args[i] == dummy_id)
               args[i] = valid_id;
          }

        /* For completeness */
        if (f->return_type == dummy_id)
          f->return_type = valid_id;

        f++;
     }
}

/*}}}*/

SLANG_MODULE(nonthermal);
int init_nonthermal_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;
   SLang_Class_Type *cl;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return -1;

   if (Bspline_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Bspline_Type")))
          return -1;
        (void) SLclass_set_destroy_function (cl, destroy_bspline_type);

        /* By registering as SLANG_VOID_TYPE,
         * slang will dynamically allocate a type
         */
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Bspline_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return -1;

        Bspline_Type_Id = SLclass_get_class_id (cl);
        patchup_intrinsic_table (Bspline_Type_Id);
     }

   if ((-1 == SLns_add_intrin_fun_table (ns, Intrinsics, NULL))
        || (-1 == SLns_add_intrin_var_table (ns, Intrin_Variables, NULL))
        || (-1 == SLns_add_iconstant_table (ns, Intrin_Const, NULL)))
     return -1;

   return 0;
}

