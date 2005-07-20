/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include <slang.h>

#include "_nonthermal.h"
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
/* Note that fitted norm = n_GeV * V / (4*pi*D^2) */

static int pop_density_info (Density_Info *di) /*{{{*/
{
   Particle_Type *p = &di->particle;
   int particle_type;

   if (-1 == POP_DOUBLE (&di->n_th)
       || -1 == POP_DOUBLE (&di->kT)
       || -1 == POP_DOUBLE (&di->n_GeV)
       || -1 == POP_DOUBLE (&p->cutoff_energy)
       || -1 == POP_DOUBLE (&p->curvature)
       || -1 == POP_DOUBLE (&p->index)
       || -1 == SLang_pop_integer (&particle_type))
     {
        return -1;
     }

   (void) init_particle_spectrum (p);

   switch (particle_type)
     {
      case PROTON:
        p->mass = GSL_CONST_CGSM_MASS_PROTON;
        break;

      default:
        p->mass = GSL_CONST_CGSM_MASS_ELECTRON;
        break;
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

   /* match units of nonthermal electron distribution function */
   return f * (GEV / GSL_CONST_CGSM_SPEED_OF_LIGHT);
}

/*}}}*/

static double particle_distrib (double *pc, double *index, /*{{{*/
                                double *curvature, double *cutoff_energy,
                                double *mass)
{
   Particle_Type pt;
   double f;

   (void) init_particle_spectrum (&pt);

   pt.index = *index;
   pt.curvature = *curvature;
   pt.cutoff_energy = *cutoff_energy;
   pt.mass = *mass;

   /* dn/d(Pc) */
   (void) (*pt.spectrum)(&pt, *pc, &f);

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

   f_th = thermal_distrib (&mv, &di->kT, &pt->mass);
   f_th *= di->n_th;

   gamma2 = 1.0 /((1.0 + beta) * (1.0 - beta));
   gamma = sqrt(gamma2);
   pc = gamma * mv * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   /* dn/d(Pc) */
   (void) (*pt->spectrum) (pt, pc, &f_nth);
   f_nth *= di->n_GeV;
   /* dn/dp = d(Pc)/dp * dn/d(Pc);   p = m*v, P = gamma*mv */
   f_nth *= gamma2 * gamma;
   
   /* fprintf (stderr, "f_th = %g  f_nth = %g\n", f_th, f_nth); */

   return f_th - f_nth;
}

/*}}}*/

static int find_momentum_min (Density_Info *di, double *momentum) /*{{{*/
{
   Particle_Type *pt = &di->particle;
   double mc, p_th, pmax, p;
   unsigned int n = 8;
   int status = 1;

   /* FIXME: This won't handle pathological cases, e.g. where
    * the two curves cross at a tangent point or where there are
    * two roots above the thermal peak.  I'm not worried about
    * that right now though. */

   p_th = sqrt (2 * pt->mass * di->kT * KEV);
   mc = pt->mass * GSL_CONST_CGSM_SPEED_OF_LIGHT;
   
   pmax = p_th * 4;

   while (n-- > 0 && pmax < mc)
     {
        if (0 == (status = bisection (&root_func, p_th, pmax, di, &p)))
          break;

        pmax *= 2;
     }

   /* No solution, use thermal peak */
   if (status || p > mc)
     p = p_th;

   *momentum = p;

   return status;
}

/*}}}*/

static double _find_momentum_min (void) /*{{{*/
{
   Density_Info di;
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
   status = gsl_integration_qag (f, pc_min, pc_max, epsabs, epsrel, limit,
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

   /* dn/d(Pc) */
   (void) (*pt->spectrum)(pt, pc, &ne);

   /* relativistic momentum = Pc = gamma * m * v */
   x = pc / mc2;
   gamma = sqrt (1.0 + x * x);
   e_kinetic = (gamma - 1.0) * mc2;

   return ne * e_kinetic;
}

/*}}}*/

static double nontherm_density_integrand (double pc, void *cd) /*{{{*/
{
   Density_Info *di = (Density_Info *)cd;
   Particle_Type *pt = &di->particle;
   double ne;

   /* relativistic momentum = Pc = gamma * m * v */

   /* dn/d(Pc) */
   (void) (*pt->spectrum)(pt, pc, &ne);

   return ne;
}

/*}}}*/

static double nontherm_energy_density (void) /*{{{*/
{
   Density_Info di;
   double ne;

   if (-1 == pop_density_info (&di))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   if (-1 == nontherm_integral (&di, &nontherm_energy_density_integrand, &ne))
     ne = 0.0;

   return ne;
}
/*}}}*/

static double nontherm_density (void) /*{{{*/
{
   Density_Info di;
   double n;

   if (-1 == pop_density_info (&di))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   if (-1 == nontherm_integral (&di, &nontherm_density_integrand, &n))
     n = 0.0;

   return n;
}
/*}}}*/

static double Electron_Density;
static double charge_balance (double p_norm, void *cd) /*{{{*/
{
   Density_Info *dp = (Density_Info *) cd;
   double fp;

   dp->n_GeV = p_norm;

   if (-1 == nontherm_integral (dp, &nontherm_density_integrand, &fp))
     return -1.0;

   return 1.0 - fp * p_norm / Electron_Density;
}

/*}}}*/

static double conserve_charge (void) /*{{{*/
{
   Density_Info de, dp;
   double p_norm, e_norm, fe, lo, hi;
   int k, status;
   char *msg = "";

   if (-1 == pop_density_info (&dp)
       || -1 == pop_density_info (&de))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1.0;
     }

   e_norm = de.n_GeV;
   if (-1 == nontherm_integral (&de, &nontherm_density_integrand, &fe))
     {
        if (Failed_Finding_Momentum_Min) msg = ":  failed finding min momentum";
        fprintf (stderr, "failed computing e- density%s\n", msg);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1.0;
     }

   Electron_Density = e_norm * fe;

   lo = e_norm;
   hi = e_norm * 2.e3;
   k = 0;
   do
     {
        status = bisection (&charge_balance, lo, hi, &dp, &p_norm);
        hi *= 1.e2;
        lo *= 1.e-2;
        k++;
     }
   while (status != 0 && k < 8);

   if (status || Failed_Finding_Momentum_Min)
     {
        if (Failed_Finding_Momentum_Min) msg = ":  failed finding min momentum";
        fprintf (stderr, "failed computing proton density%s\n", msg);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1.0;
     }   

   return p_norm;
}

/*}}}*/

static void synchrotron1 (void) /*{{{*/
{
   SLang_Array_Type *sl_x = NULL;
   SLang_Array_Type *sl_f = NULL;
   double *x, *y;
   int i, n;

   if ((-1 == SLang_pop_array_of_type (&sl_x, SLANG_DOUBLE_TYPE))
       || (sl_x == NULL))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   x = (double *)sl_x->data;
   n = sl_x->num_elements;

   if (NULL == (sl_f = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
     {
        SLang_free_array (sl_x);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   y = (double *)sl_f->data;

   for (i = 0; i < n; i++)
     {
        gsl_sf_result f;
        (void) gsl_sf_synchrotron_1_e (x[i], &f);
        y[i] = f.val;
     }

   SLang_free_array (sl_x);

   if (n == 1)
     {
        SLang_push_double (y[0]);
        SLang_free_array (sl_f);
     }
   else SLang_push_array (sl_f, 1);
}

/*}}}*/

static double gamma_function (double *x) /*{{{*/
{
   gsl_sf_result y;
   gsl_error_handler_t *gsl_error_handler;

   gsl_error_handler = gsl_set_error_handler_off ();
   if (-1 == gsl_sf_gamma_e (*x, &y))
     SLang_set_error (SL_INTRINSIC_ERROR);
   gsl_set_error_handler (gsl_error_handler);

   return y.val;
}

/*}}}*/

#define D SLANG_DOUBLE_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Fun_Type Intrinsics [] =
{
   MAKE_INTRINSIC_1("_gamma", gamma_function, D, D),
   MAKE_INTRINSIC("sync_F", synchrotron1, V, 0),
   MAKE_INTRINSIC_3("thermal_distrib", thermal_distrib, D, D, D, D),
   MAKE_INTRINSIC_5("particle_distrib", particle_distrib, D, D, D, D, D, D),
   MAKE_INTRINSIC("conserve_charge", conserve_charge, D, 0),
   MAKE_INTRINSIC("_find_momentum_min", _find_momentum_min, D, 0),
   MAKE_INTRINSIC("_nontherm_density", nontherm_density, D, 0),
   MAKE_INTRINSIC("_nontherm_energy_density", nontherm_energy_density, D, 0),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef V

static SLang_IConstant_Type Intrin_Const [] =
{
   MAKE_ICONSTANT("GSL", GSL),
   MAKE_ICONSTANT("FAST", FAST),
   MAKE_ICONSTANT("PROTON", PROTON),
   MAKE_ICONSTANT("ELECTRON", ELECTRON),
   MAKE_ICONSTANT("_nonthermal_module_version", MODULE_VERSION_NUMBER),
   SLANG_END_ICONST_TABLE
};

static char *Module_Install_Prefix = INSTALL_PREFIX ;

static SLang_Intrin_Var_Type Intrin_Variables [] =
{
   MAKE_VARIABLE("_nonthermal_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("_nonthermal_install_prefix", &Module_Install_Prefix, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

SLANG_MODULE(nonthermal);
int init_nonthermal_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return -1;

   if ((-1 == SLns_add_intrin_fun_table (ns, Intrinsics, NULL))
        || (-1 == SLns_add_intrin_var_table (ns, Intrin_Variables, NULL))
        || (-1 == SLns_add_iconstant_table (ns, Intrin_Const, NULL)))
     return -1;

   return 0;
}

