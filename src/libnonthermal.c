/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <slang.h>

#include "isis.h"

#include "photons.h"
#include "synchrotron.h"
#include "syn_table.h"
#include "inverse_compton.h"
#include "ic_table.h"
#include "ntbrems.h"
#include "ntb_table.h"
#include "pizero.h"
#include "pizero_table.h"

#define EV_ANGSTROM \
   ((GSL_CONST_CGSM_PLANCKS_CONSTANT_H * GSL_CONST_CGSM_SPEED_OF_LIGHT) \
       / (1.e-8 * GSL_CONST_CGSM_ELECTRON_VOLT))

typedef int Bin_Integral_Method_Type (double, double, void *,
                                      int (*)(void *, double, double *),
                                      double *);

static void *ic_client_data = NULL;
static void *sync_client_data = NULL;
static void *ntb_client_data = NULL;

static int Syn_Interpolate;
static int Syn_Bin_Integral_Method = FAST;

static int IC_Interpolate;
static int IC_Complain_On_Extrapolation;
static int IC_Bin_Integral_Method = FAST;

static int Ntb_Interpolate = 1;

static int Pizero_Interpolate = 0;

#define X_HE (0.1)
#define X_H  (1.0 - X_HE)
static double Ntb_ee_weight = (X_H + 2*X_HE);
static double Ntb_ep_weight = (X_H + 4*X_HE);   /* ep goes like Z^2 */

#if (SLANG_VERSION<20000)
int SLang_get_error (void) /*{{{*/
{
   return SLang_Error;
}

/*}}}*/
#endif

#if (SLANG_VERSION<20000)
void SLang_set_error (int err) /*{{{*/
{
   SLang_Error = err;
}

/*}}}*/
#endif

static int (*Eval_Continuum)(void *, double, double *);
static double integrand (double e, void *cl) /*{{{*/
{
   double y;
   (void)(*Eval_Continuum)(cl, e, &y);
   return y;
}

/*}}}*/

static int bin_integral_gsl (double el, double eh, void *cl, /*{{{*/
                             int (*eval_contin)(void *, double, double *),
                             double *area)
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_function f;
   double abserr;
   size_t neval;
   int ret;

   f.function = &integrand;
   f.params = cl;
   Eval_Continuum = eval_contin;

   gsl_error_handler = gsl_set_error_handler_off ();
   ret = gsl_integration_qng (&f, el, eh, 0.0, 1.e-4, area, &abserr, &neval);
   gsl_set_error_handler (gsl_error_handler);

   if (ret < 0)
     fprintf (stderr, "*** area = %12.4e  abserr = %12.4e  neval = %3d\r",
              *area, abserr, neval);

   return ret;
}

/*}}}*/

static int bin_integral_fast (double el, double eh, void *cl, /*{{{*/
                              int (*eval_contin)(void *, double, double *),
                              double *area)
{
   double de, sl, sm, sh;

   de = eh - el;

   (void) (*eval_contin) (cl, el, &sl);
   (void) (*eval_contin) (cl, eh, &sh);

#if 0
   (void) sm;
   /* trapezoid rule */
   *area = de * (sl + sh) / 2.0;
#else
   /* Simpson's rule */
   (void) (*eval_contin) (cl, 0.5*(el+eh), &sm);
   *area = de * (sl + 4*sm + sh) / 6.0;
#endif

   return 0;
}

/*}}}*/

static Bin_Integral_Method_Type *Bin_Integral_Method = &bin_integral_fast;

static void set_bin_integral_method (int method) /*{{{*/
{
   switch (method)
     {
      case GSL:
        Bin_Integral_Method = &bin_integral_gsl;
        break;

      default:
        Bin_Integral_Method = &bin_integral_fast;
        break;
     }
}

/*}}}*/

static int _nt_binned_contin (void *cl, /*{{{*/
                              int (*eval_contin)(void *, double, double *),
                              double *val, Isis_Hist_t *g, double *par, unsigned int npar)
{
   double norm = par[0];
   int i;

   (void) npar;

   if (Bin_Integral_Method == NULL)
     return -1;

    for (i=0; i < g->n_notice; i++)
     {
        double el, eh, area;
        int n;

        n = g->notice_list[i];

        /* Angstrom -> eV */
        el = EV_ANGSTROM / g->bin_hi[n];
        eh = EV_ANGSTROM / g->bin_lo[n];

        (void) (*Bin_Integral_Method)(el, eh, cl, eval_contin, &area);

        /* norm = (V/4*pi*D^2) * n0 (cm^-2 GeV^-1)
         * s(E)dE = s(y)dy = photons cm^-2 s^-1
         */
        val[i] = norm * (area * 1.e-9);
     }

   return 0;
}

/*}}}*/

static int _nt_contin (void *cl, /*{{{*/
                       int (*eval_contin)(void *, double, double *),
                       double *val, Isis_User_Grid_t *g, double *par, unsigned int npar)
{
   double norm = par[0];
   int i;

   (void) npar;

   /* g->x[i] in eV
    * val in photons /s /cm^2 /GeV
    */

    for (i=0; i < g->npts; i++)
     {
        double s;
        (void) (*eval_contin) (cl, g->x[i], &s);
        val[i] = norm * s;
     }

   return 0;
}

/*}}}*/

static double _sync_angular_integral (double *x, int *interpolate) /*{{{*/
{
   double y;

   if (*interpolate)
     (void) syn_interp_angular_integral (sync_client_data, *x, &y);
   else
     (void) syn_angular_integral (*x, &y);

   return y;
}

/*}}}*/

static void init_sync (double *par, Synchrotron_Type *s, Particle_Type *elec) /*{{{*/
{
   (void) init_particle_spectrum (elec);

   elec->index = par[2];
   elec->curvature = par[3];
   elec->cutoff_energy = par[4];
   elec->mass = GSL_CONST_CGSM_MASS_ELECTRON;

   s->B_tot = par[1] * 1.e-6;  /* Gauss */
   s->electrons = elec;
   s->client_data = sync_client_data;
   s->interpolate = Syn_Interpolate;

   set_bin_integral_method (Syn_Bin_Integral_Method);
}

/*}}}*/

static int binned_sync (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Synchrotron_Type s = NULL_SYNCHROTRON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;

   init_sync (par, &s, &elec);

   return _nt_binned_contin ((void *)&s, &syn_calc_synchrotron,
                              val, g, par, npar);
}

/*}}}*/

static int unbinned_sync (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Synchrotron_Type s = NULL_SYNCHROTRON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;

   init_sync (par, &s, &elec);

   return _nt_contin ((void *)&s, &syn_calc_synchrotron,
                       val, g, par, npar);
}

/*}}}*/

static void _invc_make_table (double *t) /*{{{*/
{
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   int status;

   if (*t <= 0) *t = CBR_TEMPERATURE;

   fprintf (stdout, "Generating IC table for T = %g K\n", *t);
   fflush (stdout);

   set_incident_photon_kelvin_temperature (*t);

   ic.electrons = NULL;
   ic.incident_photons = &incident_photon_spectrum;
   ic.incident_photon_max_energy = incident_photon_max_energy ();

   status = ic_push_table (&ic);

   if (status)
     fprintf (stderr, "*** failed computing IC table for T = %g K\n", *t);
}

/*}}}*/

static double _invc_photon_integral (double *gamma, double *energy_final_photon, /*{{{*/
                                    double *t, int *interpolate)
{
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   double x;

   if (*t <= 0) *t = CBR_TEMPERATURE;
   set_incident_photon_kelvin_temperature (*t);

   ic.energy_final_photon = *energy_final_photon;  /* E/(mc^2) */
   ic.electron_gamma = *gamma;
   ic.electrons = NULL;
   ic.incident_photons = &incident_photon_spectrum;
   ic.incident_photon_max_energy = incident_photon_max_energy();
   ic.client_data = ic_client_data;

   if (*interpolate)
     (void) ic_interp_photon_integral (&ic, &x);
   else
     (void) ic_integral_over_incident_photons (&ic, &x);

   return x;
}

/*}}}*/

static void init_invc (double *par, Inverse_Compton_Type *ic, Particle_Type *elec) /*{{{*/
{
   (void) init_particle_spectrum (elec);

   elec->index = par[1];
   elec->curvature = par[2];
   elec->cutoff_energy = par[3];
   elec->mass = GSL_CONST_CGSM_MASS_ELECTRON;

   set_incident_photon_kelvin_temperature (par[4]);

   ic->electrons = elec;
   ic->incident_photons = &incident_photon_spectrum;
   ic->incident_photon_max_energy = incident_photon_max_energy ();
   ic->client_data = ic_client_data;
   ic->interpolate = IC_Interpolate;
   ic->complain_on_extrapolate = IC_Complain_On_Extrapolation;

   set_bin_integral_method (IC_Bin_Integral_Method);
}

/*}}}*/

static int binned_invc (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;

   init_invc (par, &ic, &elec);

   return _nt_binned_contin ((void *)&ic, &ic_calc_inverse_compton,
                              val, g, par, npar);
}

/*}}}*/

static int unbinned_invc (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;

   init_invc (par, &ic, &elec);

   return _nt_contin ((void *)&ic, &ic_calc_inverse_compton,
                       val, g, par, npar);
}

/*}}}*/

static void sync_free_client_data (void) /*{{{*/
{
   syn_free_table (sync_client_data);
   sync_client_data = NULL;
}

/*}}}*/

#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Var_Type Sync_Intrin_Vars [] =
{
   MAKE_VARIABLE("Syn_Interpolate", &Syn_Interpolate, I, 0),
   MAKE_VARIABLE("Syn_Bin_Integral_Method", &Syn_Bin_Integral_Method, I, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Fun_Type Sync_Intrinsics [] =
{
   MAKE_INTRINSIC_2("_sync_angular_integral", _sync_angular_integral, D, D,I),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef I
#undef S
#undef V

static int sync_init_client_data (char *file) /*{{{*/
{
   int status = 0;

   if (file == NULL)
     return -1;

   if ((-1 == SLns_add_intrin_var_table (NULL, Sync_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_fun_table (NULL, Sync_Intrinsics, NULL)))
     return -1;

   if (sync_client_data != NULL)
     sync_free_client_data ();

   if (NULL == (sync_client_data = syn_load_table (file)))
     status = -1;

   if (status == 0)
     Syn_Interpolate = 1;

   return status;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(sync,p,options) /*{{{*/
{
   static const char *parameter_names[] = {"norm", "B_tot", "index", "curvature", "cutoff", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", "microgauss", "",          "",    "TeV", NULL};
   static double default_max[]   = {10.0, 10.0, -1.0,  1.0, 100.0};
   static double default_value[] = { 1.0,  1.0, -2.0,  0.0, 10.0};
   static double default_min[]   = { 0.0,  0.0, -3.0, -1.0,  1.0};
   static unsigned int default_freeze[] = {0, 0, 0, 1, 0};
   static unsigned int norm_indexes[] = {0};

   /*           V
    * norm = --------  (dn/dE)_0
    *        4 pi D^2
    * where (dn/dE)_0 = relativistic e- density @ 1 GeV
    *                   [cm^(-3) GeV^(-1)]
    */

   p->function_exit = sync_free_client_data;
   p->binned = binned_sync;
   p->unbinned = unbinned_sync;
   p->parameter_names = (char **) parameter_names;
   p->parameter_units = (char **) parameter_units;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_min = default_min;
   p->default_freeze = default_freeze;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   sync_init_client_data (options);

   return 0;
}

/*}}}*/

static void invc_free_client_data (void) /*{{{*/
{
   ic_free_client_data (ic_client_data);
   ic_client_data = NULL;
}

/*}}}*/

static void _invc_set_dilution_factors (void) /*{{{*/
{
   SLang_Array_Type *a;

   if ((-1 == SLang_pop_array_of_type (&a, SLANG_DOUBLE_TYPE))
       || (ic_client_data == NULL))
     SLang_set_error (SL_INTRINSIC_ERROR);

   ic_set_dilution_factors (ic_client_data, (double *)a->data, a->num_elements);

   SLang_free_array (a);
}

/*}}}*/

#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Var_Type Invc_Intrin_Vars [] =
{
   MAKE_VARIABLE("IC_Interpolate", &IC_Interpolate, I, 0),
   MAKE_VARIABLE("IC_Complain_On_Extrapolation", &IC_Complain_On_Extrapolation, I, 0),
   MAKE_VARIABLE("IC_Bin_Integral_Method", &IC_Bin_Integral_Method, I, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Fun_Type Invc_Intrinsics [] =
{
   MAKE_INTRINSIC("_invc_set_dilution_factors", _invc_set_dilution_factors, V, 0),
   MAKE_INTRINSIC_1("_invc_make_table", _invc_make_table, V, D),
   MAKE_INTRINSIC_4("_invc_photon_integral", _invc_photon_integral, D, D,D,D,I),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef I
#undef S
#undef V

static int invc_init_client_data (char *file) /*{{{*/
{
   int status = 0;

   if (file == NULL)
     return -1;

   if ((-1 == SLns_add_intrin_var_table (NULL, Invc_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_fun_table (NULL, Invc_Intrinsics, NULL)))
     return -1;

   if (IC_Interpolate != 0)
     fprintf (stderr, " IC_Interpolate=%d ignoring temperature parameter\n",
              IC_Interpolate);

   if (ic_client_data != NULL)
     invc_free_client_data ();

   if (NULL == (ic_client_data = ic_init_client_data (file)))
     status = -1;

   if (status == 0)
     IC_Interpolate = 1;

   return status;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(invc,p,options) /*{{{*/
{
   static const char *parameter_names[] = {"norm", "index", "curvature", "cutoff", "T_photon[K]", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", "",      "",    "TeV",        "K", NULL};
   static double default_max[]   = {10.0, -1.0,  1.0, 100.0, 1.e10};
   static double default_value[] = { 1.0, -2.0,  0.0, 10.0,  CBR_TEMPERATURE};
   static double default_min[]   = { 0.0, -3.0, -1.0,  1.0,   0.0};
   static unsigned int default_freeze[] = {0, 0, 1, 0, 1};
   static unsigned int norm_indexes[] = {0};

   p->function_exit = invc_free_client_data;
   p->binned = binned_invc;
   p->unbinned = unbinned_invc;
   p->parameter_names = (char **)parameter_names;
   p->parameter_units = (char **) parameter_units;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_min = default_min;
   p->default_freeze = default_freeze;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   invc_init_client_data (options);

   return 0;
}

/*}}}*/

static void _ntb_free_client_data (void) /*{{{*/
{
   ntb_free_client_data (ntb_client_data);
   ntb_client_data = NULL;
}

/*}}}*/

static void init_brem (double *par, Brems_Type *b, Particle_Type *elec) /*{{{*/
{
   (void) init_particle_spectrum (elec);

   elec->index = par[1];
   elec->curvature = par[2];
   elec->cutoff_energy = par[3];
   elec->mass = GSL_CONST_CGSM_MASS_ELECTRON;

   b->electrons = elec;
   b->ee_weight = Ntb_ee_weight;
   b->ep_weight = Ntb_ep_weight;
   b->client_data = ntb_client_data;
   b->interpolate = Ntb_Interpolate;
}

/*}}}*/

static int binned_brem (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Brems_Type b = NULL_BREMS_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;

   init_brem (par, &b, &elec);

   return _nt_binned_contin ((void *)&b, &ntb_brems,
                             val, g, par, npar);
}

/*}}}*/

static int unbinned_brem (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Brems_Type b = NULL_BREMS_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;

   init_brem (par, &b, &elec);

   return _nt_contin ((void *)&b, &ntb_brems,
                      val, g, par, npar);
}

/*}}}*/

static void _ntb_make_table (int *process) /*{{{*/
{
   int status;

   fprintf (stdout, "Generating nonthermal bremsstrahlung cross-section table\n");
   fflush (stdout);

   status = ntb_push_table (*process);

   if (status)
     fprintf (stderr, "*** failed computing table\n");
}

/*}}}*/

static int pop_2_matched_arrays (int type, SLang_Array_Type **ap, SLang_Array_Type **bp) /*{{{*/
{
   SLang_Array_Type *a, *b;

   *ap = *bp = NULL;

   if (-1 == SLang_pop_array_of_type (&b, type))
     return -1;
   if (-1 == SLang_pop_array_of_type (&a, type))
     {
        SLang_free_array (b);
        return -1;
     }

   if (a->num_elements == b->num_elements)
     {
        *ap = a;
        *bp = b;
        return 0;
     }

   fprintf (stderr, "inconsistent array sizes");
   SLang_set_error(SL_INTRINSIC_ERROR);
   SLang_free_array (a);
   SLang_free_array (b);

   return -1;
}

/*}}}*/

static void _ntb_sigma (double (*f)(double, double)) /*{{{*/
{
   SLang_Array_Type *sl_ekin, *sl_eph, *sl_sig;
   double *ekinetic, *ephoton, *sigma;
   int i, n;

   if ((f == NULL)
       || (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &sl_ekin, &sl_eph)))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   n = sl_ekin->num_elements;

   if (NULL == (sl_sig = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   ekinetic = (double *)sl_ekin->data;
   ephoton = (double *)sl_eph->data;
   sigma = (double *)sl_sig->data;

   for (i = 0; i < n; i++)
     {
        sigma[i] = (*f)(ekinetic[i], ephoton[i]);
     }

   SLang_free_array (sl_ekin);
   SLang_free_array (sl_eph);
   SLang_push_array (sl_sig, 1);
}

/*}}}*/

static void _ep_heitler (void) /*{{{*/
{
   _ntb_sigma (&_ntb_ei_sigma);
}

/*}}}*/

static void _ee_haug (void) /*{{{*/
{
   _ntb_sigma (&_ntb_ee_sigma_haug);
}

/*}}}*/

static double ee_interp_sigma (double ekinetic, double ephoton) /*{{{*/
{
   Brems_Type b = NULL_BREMS_TYPE;
   double sigma;

   /* other fields aren't used.. */
   b.client_data = ntb_client_data;
   (void) ntb_interp_sigma (&b, ekinetic, ephoton, &sigma);

   return sigma;
}

/*}}}*/

static void _ee_interp (void) /*{{{*/
{
   _ntb_sigma (&ee_interp_sigma);
}

/*}}}*/

static void _ntb_set_process_weights (double *ee, double *ep) /*{{{*/
{
   Ntb_ee_weight = *ee;
   Ntb_ep_weight = *ep;
}

/*}}}*/

#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define V SLANG_VOID_TYPE

static SLang_IConstant_Type Ntb_Intrin_Const [] =
{
   MAKE_ICONSTANT("NTB_ee", NTB_ee),
   MAKE_ICONSTANT("NTB_ep", NTB_ep),
   SLANG_END_ICONST_TABLE
};

static SLang_Intrin_Var_Type Ntb_Intrin_Vars [] =
{
   MAKE_VARIABLE("Ntb_Interpolate", &Ntb_Interpolate, I, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Fun_Type Ntb_Intrinsics [] =
{
   MAKE_INTRINSIC_1("_ntb_make_table", _ntb_make_table, V, I),
   MAKE_INTRINSIC_2("_ntb_set_process_weights", _ntb_set_process_weights, V, D,D),
   MAKE_INTRINSIC("_ep_heitler", _ep_heitler, V, 0),
   MAKE_INTRINSIC("_ee_haug", _ee_haug, V, 0),
   MAKE_INTRINSIC("_ee_interp", _ee_interp, V, 0),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef I
#undef S
#undef V

static int _ntb_init_client_data (char *options) /*{{{*/
{
   (void) options;

   if ((-1 == SLns_add_iconstant_table (NULL, Ntb_Intrin_Const, NULL))
       || (-1 == SLns_add_intrin_var_table (NULL, Ntb_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_fun_table (NULL, Ntb_Intrinsics, NULL)))
     return -1;

   if (ntb_client_data != NULL)
     _ntb_free_client_data ();

   /* NULL ok */
   ntb_client_data = ntb_init_client_data (options);

   return 0;

}

/*}}}*/

ISIS_USER_SOURCE_MODULE(ntbrem,p,options) /*{{{*/
{
   static const char *parameter_names[] = {"norm", "index", "curvature", "cutoff", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", "",       "",    "TeV", NULL};
   static double default_max[]   = {10.0, -1.0,  1.0, 100.0};
   static double default_value[] = { 1.0, -2.0,  0.0, 10.0};
   static double default_min[]   = { 0.0, -3.0, -1.0,  1.0};
   static unsigned int default_freeze[] = {0, 0, 1, 0};
   static unsigned int norm_indexes[] = {0};

   p->function_exit = NULL;
   p->binned = binned_brem;
   p->unbinned = unbinned_brem;
   p->parameter_names = (char **)parameter_names;
   p->parameter_units = (char **) parameter_units;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_min = default_min;
   p->default_freeze = default_freeze;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   _ntb_init_client_data (options);

   return 0;
}

/*}}}*/

static void pizero_diff_xsec_intrin (void);
static void pizero_distribution_intrin (void);
static double pizero_lidcs_intrin (double *T_p, double *T_pi, double *mu);

#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Var_Type Pizero_Intrin_Vars [] =
{
   MAKE_VARIABLE("Pizero_Method", &Pizero_Method, I, 0),
   MAKE_VARIABLE("Pizero_Interpolate", &Pizero_Interpolate, I, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Fun_Type Pizero_Intrinsics [] =
{
   MAKE_INTRINSIC("pizero_distribution", pizero_distribution_intrin, V, 0),
   MAKE_INTRINSIC("pizero_diff_xsec", pizero_diff_xsec_intrin, V, 0),
   MAKE_INTRINSIC_3("pizero_lidcs", pizero_lidcs_intrin, D, D, D, D),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef I
#undef V

static int pizero_init_client_data (void) /*{{{*/
{
   if ((-1 == SLns_add_intrin_var_table (NULL, Pizero_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_fun_table (NULL, Pizero_Intrinsics, NULL)))
     return -1;

   return 0;
}
/*}}}*/

static void init_pizero (double *par, Pizero_Type *p, Particle_Type *proton) /*{{{*/
{
   (void) init_particle_spectrum (proton);

   proton->index = par[1];
   proton->curvature = par[2];
   proton->cutoff_energy = par[3];
   proton->mass = GSL_CONST_CGSM_MASS_PROTON;
   
   p->protons = proton;
   p->interpolate = Pizero_Method ? Pizero_Interpolate : 0;
   p->client_data = pizero_alloc_table (PIZERO_TABLE_SIZE);
}

/*}}}*/

static void pizero_distribution_intrin (void) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   SLang_Array_Type *sl_par = NULL;
   SLang_Array_Type *sl_energies = NULL;
   SLang_Array_Type *sl_qpi = NULL;
   double *par, *energies, *qpi;
   int i, n;

   if ((-1 == SLang_pop_array_of_type (&sl_par, SLANG_DOUBLE_TYPE))
       || (sl_par == NULL)
       || (sl_par->num_elements != 4))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);        
        return;
     }   
   
   if ((-1 == SLang_pop_array_of_type (&sl_energies, SLANG_DOUBLE_TYPE))
       || (sl_energies == NULL))
     {
        SLang_free_array (sl_par);
        SLang_set_error (SL_INTRINSIC_ERROR);        
        return;
     }   
   
   n = sl_energies->num_elements;
   
   if (NULL == (sl_qpi = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
     {
        SLang_free_array (sl_energies);
        SLang_free_array (sl_par);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }   

   par = (double *)sl_par->data;
   init_pizero (par, &p, &proton);
   
   energies = (double *)sl_energies->data;
   qpi = (double *)sl_qpi->data;
   
   for (i = 0; i < n; i++)
     {
        p.energy = energies[i] * GSL_CONST_CGSM_ELECTRON_VOLT;
        if (-1 == pizero_distribution (&p, &qpi[i]))
          qpi[i] = 0.0;
     }

   pizero_free_table (p.client_data);   
   
   SLang_push_array (sl_qpi, 1);   
   SLang_free_array (sl_par);
   SLang_free_array (sl_energies);
}

/*}}}*/

static double pizero_lidcs_intrin (double *T_p, double *T_pi, double *mu)
{
   return pizero_lidcs (*T_p, *T_pi, *mu);
}

static void pizero_diff_xsec_intrin (void) /*{{{*/
{
   SLang_Array_Type *sl_tp, *sl_tpi, *sl_sig;
   double *T_p, *T_pi, *sigma;
   int i, n;

   if (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &sl_tp, &sl_tpi))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   n = sl_tp->num_elements;

   if (NULL == (sl_sig = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   T_p = (double *)sl_tp->data;
   T_pi = (double *)sl_tpi->data;
   sigma = (double *)sl_sig->data;

   for (i = 0; i < n; i++)
     {
        sigma[i] = pizero_differential_xsec (T_p[i], T_pi[i]);
     }

   SLang_free_array (sl_tp);
   SLang_free_array (sl_tpi);
   SLang_push_array (sl_sig, 1);
}

/*}}}*/

static int binned_pizero (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   int status;

   init_pizero (par, &p, &proton);

   status = _nt_binned_contin ((void *)&p, &pizero_decay,
                               val, g, par, npar);
   pizero_free_table (p.client_data);

   return status;
}

/*}}}*/

static int unbinned_pizero (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   int status;

   init_pizero (par, &p, &proton);

   status = _nt_contin ((void *)&p, &pizero_decay,
                        val, g, par, npar);
   pizero_free_table (p.client_data);

   return status;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(pizero,p,options) /*{{{*/
{
   static const char *parameter_names[] = {"norm", "index", "curvature", "cutoff", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", "",       "",    "TeV", NULL};
   static double default_max[]   = {10.0, -1.0,  1.0, 100.0};
   static double default_value[] = { 1.0, -2.0,  0.0, 10.0};
   static double default_min[]   = { 0.0, -3.0, -1.0,  1.0};
   static unsigned int default_freeze[] = {0, 0, 1, 0};
   static unsigned int norm_indexes[] = {0};

   (void) options;

   p->function_exit = NULL;
   p->binned = binned_pizero;
   p->unbinned = unbinned_pizero;
   p->parameter_names = (char **)parameter_names;
   p->parameter_units = (char **) parameter_units;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_min = default_min;
   p->default_freeze = default_freeze;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   pizero_init_client_data ();
   return 0;
}

/*}}}*/

