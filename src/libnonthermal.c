/* -*- mode: C; mode: fold -*- */
/*
  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008 John C. Houck

  This file is part of the nonthermal module

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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

/* Relative error tolerances for accepting numerical integral results */
double Sync_Epsrel = 1.e-6;
double Invc_Epsrel = 1.e-6;
double Ntb_Epsrel  = 1.e-6;
double Pizero_Epsrel = 1.e-6;

int Ntb_Integration_Method = 0;

static void *ic_client_data = NULL;
static void *sync_client_data = NULL;
static void *ntb_client_data = NULL;

static int Ntb_Intrinsics_Initialized;

static int Syn_Interpolate;
static int IC_Interpolate;
static int Ntb_Interpolate = 1;

static int IC_Complain_On_Extrapolation;

#define X_HE (0.1/1.1)
#define X_H  (1.0/1.1)
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

static int _nt_binned_contin (void *cl, /*{{{*/
                              int (*eval_contin)(void *, double, double *),
                              double *val, Isis_Hist_t *g, double *par, unsigned int npar)
{
   double norm = par[0];
   double saved_el, saved_sl;
   int i;

   (void) npar;

   if (norm == 0.0)
     {
        for (i = 0; i < g->n_notice; i++)
          {
             val[i] = 0.0;
          }

        return 0;
     }

   saved_el = DBL_MAX;
   saved_sl = 0.0;

    for (i=0; i < g->n_notice; i++)
     {
        double el, eh, sl, sm, sh, area;
        int n;

        if (isis_user_break() || (0 != SLang_get_error ()))
          return -1;

        n = g->notice_list[i];

        /* Angstrom -> eV */
        el = EV_ANGSTROM / g->bin_hi[n];
        eh = EV_ANGSTROM / g->bin_lo[n];

        (void) (*eval_contin)(cl, el, &sl);
        (void) (*eval_contin)(cl, 0.5*(el+eh), &sm);

        /* caching bin-edge value saves 30% */
        if (eh != saved_el)
          (void) (*eval_contin)(cl, eh, &sh);
        else sh = saved_sl;

        saved_el = el;
        saved_sl = sl;

        /* Simpson's rule */
        area = (eh - el) * (sl + 4*sm + sh) / 6.0;

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

   if (norm == 0.0)
     {
        for (i = 0; i < g->npts; i++)
          {
             val[i] = 0.0;
          }
        return 0;
     }

   /* g->x[i] in eV
    * val in photons /s /cm^2 /GeV
    */

    for (i=0; i < g->npts; i++)
     {
        double s;
        if (isis_user_break() || (0 != SLang_get_error ()))
          return -1;
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

static int init_sync (double *par, unsigned int npar, Synchrotron_Type *s, Particle_Type *elec) /*{{{*/
{
   (void) npar;

   if (-1 == init_pdf (elec, ELECTRON))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   s->B_tot = par[1] * 1.e-6;  /* Gauss */
   s->electrons = elec;
   s->client_data = sync_client_data;
   s->interpolate = Syn_Interpolate;

   return 0;
}

/*}}}*/

static int binned_sync (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Synchrotron_Type s = NULL_SYNCHROTRON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_sync (par, npar, &s, &elec))
     return -1;

   status = _nt_binned_contin ((void *)&s, &syn_calc_synchrotron,
                               val, g, par, npar);
   free_pdf (&elec);

   return status;
}

/*}}}*/

static int unbinned_sync (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Synchrotron_Type s = NULL_SYNCHROTRON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_sync (par, npar, &s, &elec))
     return -1;

   status = _nt_contin ((void *)&s, &syn_calc_synchrotron,
                        val, g, par, npar);
   free_pdf (&elec);

   return status;
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
   ic.incident_photon_min_energy = incident_photon_min_energy();
   ic.client_data = ic_client_data;

   if (*interpolate)
     (void) ic_interp_photon_integral (&ic, &x);
   else
     (void) ic_integral_over_incident_photons (&ic, &x);

   return x;
}

/*}}}*/

static int init_invc (double *par, unsigned int npar, Inverse_Compton_Type *ic, Particle_Type *elec) /*{{{*/
{
   (void) npar;

   if (-1 == init_pdf (elec, ELECTRON))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   set_incident_photon_kelvin_temperature (par[1]);

   ic->electrons = elec;
   ic->incident_photons = &incident_photon_spectrum;
   ic->incident_photon_max_energy = incident_photon_max_energy ();
   ic->client_data = ic_client_data;
   ic->interpolate = IC_Interpolate;
   ic->complain_on_extrapolate = IC_Complain_On_Extrapolation;

   return 0;
}

/*}}}*/

static int binned_invc (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_invc (par, npar, &ic, &elec))
     return -1;

   status = _nt_binned_contin ((void *)&ic, &ic_calc_inverse_compton,
                               val, g, par, npar);
   free_pdf (&elec);

   return status;
}

/*}}}*/

static int unbinned_invc (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_invc (par, npar, &ic, &elec))
     return -1;

   status = _nt_contin ((void *)&ic, &ic_calc_inverse_compton,
                        val, g, par, npar);
   free_pdf (&elec);

   return status;
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
   MAKE_VARIABLE("_sync_interpolate", &Syn_Interpolate, I, 0),
   MAKE_VARIABLE("_sync_epsrel", &Sync_Epsrel, D, 0),
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
   static const char *parameter_names[] = {"norm", "B_tot", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", "microgauss", NULL};
   static double default_max[]   = {1.0e10, 1.0e3};
   static double default_value[] = { 1.0,  1.0};
   static double default_min[]   = { 0.0,  0.0};
   static double default_relstep[] = {1.e-4, 1.e-4};
   static unsigned int default_freeze[] = {0, 0};
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
   p->default_relstep = default_relstep;
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

static double invc_knlimit_constant_intrin (double *p)
{
   return ic_knlimit_constant (*p);
}

static int invc_load_table (char *file) /*{{{*/
{
   if (ic_client_data != NULL)
     invc_free_client_data ();

   fprintf (stderr, "invc_load_table:  loading %s\n", file);
   if (NULL == (ic_client_data = ic_init_client_data (file)))
     return -1;

   return 0;
}

/*}}}*/

#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Var_Type Invc_Intrin_Vars [] =
{
   MAKE_VARIABLE("_invc_interpolate", &IC_Interpolate, I, 0),
   MAKE_VARIABLE("_invc_complain_on_extrapolation", &IC_Complain_On_Extrapolation, I, 0),
   MAKE_VARIABLE("_invc_epsrel", &Invc_Epsrel, D, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Fun_Type Invc_Intrinsics [] =
{
   MAKE_INTRINSIC("_invc_set_dilution_factors", _invc_set_dilution_factors, V, 0),
   MAKE_INTRINSIC_4("_invc_photon_integral", _invc_photon_integral, D, D,D,D,I),
   MAKE_INTRINSIC_1("_invc_knlimit_constant_intrin", invc_knlimit_constant_intrin, D, D),
   MAKE_INTRINSIC_1("invc_load_table_intrin", invc_load_table, I, S),
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
   static const char *parameter_names[] = {"norm", "T_photon", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", "K", NULL};
   static double default_max[]   = {1.0e10, 1.e10};
   static double default_value[] = { 1.0, CBR_TEMPERATURE};
   static double default_min[]   = { 0.0, 0.0};
   static double default_relstep[] = {1.e-4, 1.e-4};
   static unsigned int default_freeze[] = {0, 1};
   static unsigned int norm_indexes[] = {0};

   p->function_exit = invc_free_client_data;
   p->binned = binned_invc;
   p->unbinned = unbinned_invc;
   p->parameter_names = (char **)parameter_names;
   p->parameter_units = (char **) parameter_units;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_min = default_min;
   p->default_relstep = default_relstep;
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

static int init_brem (double *par, unsigned int npar, Brems_Type *b, Particle_Type *elec) /*{{{*/
{
   (void) par; (void) npar;

   if (-1 == init_pdf (elec, ELECTRON))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   b->electrons = elec;
   b->ee_weight = Ntb_ee_weight;
   b->ep_weight = Ntb_ep_weight;
   b->client_data = ntb_client_data;
   b->interpolate = Ntb_Interpolate;

   return 0;
}

/*}}}*/

static int binned_brem (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Brems_Type b = NULL_BREMS_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_brem (par, npar, &b, &elec))
     return -1;

   status = _nt_binned_contin ((void *)&b, &ntb_brems, val, g, par, npar);
   free_pdf (&elec);

   return status;
}

/*}}}*/

static int unbinned_brem (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Brems_Type b = NULL_BREMS_TYPE;
   Particle_Type elec = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_brem (par, npar, &b, &elec))
     return -1;

   status = _nt_contin ((void *)&b, &ntb_brems, val, g, par, npar);
   free_pdf (&elec);

   return status;
}

/*}}}*/

static int pop_2_matched_arrays (SLtype type, SLang_Array_Type **ap, SLang_Array_Type **bp) /*{{{*/
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

static double _ep_heitler1 (double *ekinetic, double *ephoton) /*{{{*/
{
   return _ntb_ei_sigma (*ekinetic, *ephoton);
}

/*}}}*/

static double _ee_haug1 (double *ekinetic, double *ephoton) /*{{{*/
{
   return _ntb_ee_sigma_haug (*ekinetic, *ephoton);
}

/*}}}*/

static double _ee_haug1_lab (double *ekinetic, double *ephoton) /*{{{*/
{
   return _ntb_ee_sigma_haug_lab (*ekinetic, *ephoton);
}

/*}}}*/

static double _ee_interp1 (double *ekinetic, double *ephoton)
{
   Brems_Type b = NULL_BREMS_TYPE;
   double sigma;

   /* other fields aren't used.. */
   b.client_data = ntb_client_data;
   (void) ntb_interp_sigma (&b, *ekinetic, *ephoton, &sigma);

   return sigma;
}

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
   MAKE_VARIABLE("_ntbrem_interpolate", &Ntb_Interpolate, I, 0),
   MAKE_VARIABLE("_ntbrem_integration_method", &Ntb_Integration_Method, I, 0),
   MAKE_VARIABLE("_ntbrem_epsrel", &Ntb_Epsrel, D, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Fun_Type Ntb_Intrinsics [] =
{
   MAKE_INTRINSIC_2("_ntb_set_process_weights", _ntb_set_process_weights, V, D,D),
   MAKE_INTRINSIC("_ep_heitler", _ep_heitler, V, 0),
   MAKE_INTRINSIC("_ee_haug", _ee_haug, V, 0),
   MAKE_INTRINSIC("_ee_interp", _ee_interp, V, 0),
   MAKE_INTRINSIC_2("_ep_heitler1", _ep_heitler1, D, D,D),
   MAKE_INTRINSIC_2("_ee_haug1", _ee_haug1, D, D,D),
   MAKE_INTRINSIC_2("_ee_haug1_lab", _ee_haug1_lab, D, D,D),
   MAKE_INTRINSIC_2("_ee_interp1", _ee_interp1, D, D,D),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef D
#undef I
#undef S
#undef V

static int _ntb_init_client_data (char *options) /*{{{*/
{
   (void) options;

   if (Ntb_Intrinsics_Initialized)
     return 0;

   if ((-1 == SLns_add_iconstant_table (NULL, Ntb_Intrin_Const, NULL))
       || (-1 == SLns_add_intrin_var_table (NULL, Ntb_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_fun_table (NULL, Ntb_Intrinsics, NULL)))
     return -1;

   if (ntb_client_data != NULL)
     _ntb_free_client_data ();

   /* NULL ok */
   ntb_client_data = ntb_init_client_data (options);

   Ntb_Intrinsics_Initialized = 1;

   return 0;

}

/*}}}*/

ISIS_USER_SOURCE_MODULE(ntbrem,p,options) /*{{{*/
{
   static const char *parameter_names[] = {"norm", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", NULL};
   static double default_max[]   = {1.e10};
   static double default_value[] = { 1.0};
   static double default_min[]   = { 0.0};
   static double default_relstep[] = {1.e-4};
   static unsigned int default_freeze[] = {0};
   static unsigned int norm_indexes[] = {0};

   p->function_exit = NULL;
   p->binned = binned_brem;
   p->unbinned = unbinned_brem;
   p->parameter_names = (char **)parameter_names;
   p->parameter_units = (char **) parameter_units;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_min = default_min;
   p->default_relstep = default_relstep;
   p->default_freeze = default_freeze;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   _ntb_init_client_data (options);

   return 0;
}

/*}}}*/

#ifdef USE_CPARAMLIB
#include "pizero_cparamlib.c"
#else

static int init_pizero (double *par, unsigned int npar, Pizero_Type *p, Particle_Type *proton);

static void pizero_distribution_intrin (void) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   SLang_Array_Type *sl_par = NULL;
   SLang_Array_Type *sl_energies = NULL;
   SLang_Array_Type *sl_qpi = NULL;
   double *par, *energies, *qpi;
   unsigned int npar;
   int i, n;

   if ((-1 == SLang_pop_array_of_type (&sl_par, SLANG_DOUBLE_TYPE))
       || (sl_par == NULL)
       || (sl_par->num_elements != 4))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   npar = sl_par->num_elements;

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
   if (-1 == init_pizero (par, npar, &p, &proton))
     {
        SLang_free_array (sl_energies);
        SLang_free_array (sl_par);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   energies = (double *)sl_energies->data;
   qpi = (double *)sl_qpi->data;

   for (i = 0; i < n; i++)
     {
        p.energy = energies[i] * GSL_CONST_CGSM_ELECTRON_VOLT;
        if (-1 == pizero_distribution (&p, &qpi[i]))
          qpi[i] = 0.0;
     }

   pizero_free_table (p.client_data);
   free_pdf (&proton);

   SLang_push_array (sl_qpi, 1);
   SLang_free_array (sl_par);
   SLang_free_array (sl_energies);
}

/*}}}*/

static double pizero_lidcs_intrin (double *T_p, double *T_pi, double *mu) /*{{{*/
{
   return pizero_lidcs (*T_p, *T_pi, *mu);
}

/*}}}*/

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

static int Pizero_Interpolate = 0;

#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Var_Type Pizero_Intrin_Vars [] =
{
   MAKE_VARIABLE("_pizero_method", &Pizero_Method, I, 0),
   MAKE_VARIABLE("_pizero_interpolate", &Pizero_Interpolate, I, 0),
   MAKE_VARIABLE("_pizero_epsrel", &Pizero_Epsrel, D, 0),
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
       || (-1 == SLns_add_intrin_fun_table (NULL, Pizero_Intrinsics, NULL))
       )
     return -1;

   return 0;
}
/*}}}*/

static int init_pizero (double *par, unsigned int npar, Pizero_Type *p, Particle_Type *proton) /*{{{*/
{
   (void) par; (void) npar;

   if (-1 == init_pdf (proton, PROTON))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   p->protons = proton;
   p->interpolate = Pizero_Method ? Pizero_Interpolate : 0;
   p->client_data = pizero_alloc_table (PIZERO_TABLE_SIZE);

   return 0;
}

/*}}}*/

static int binned_pizero (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_pizero (par, npar, &p, &proton))
     return -1;

   status = _nt_binned_contin ((void *)&p, &pizero_decay, val, g, par, npar);
   pizero_free_table (p.client_data);
   free_pdf (&proton);

   return status;
}

/*}}}*/

static int unbinned_pizero (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   int status;

   if (-1 == init_pizero (par, npar, &p, &proton))
     return -1;

   status = _nt_contin ((void *)&p, &pizero_decay, val, g, par, npar);
   pizero_free_table (p.client_data);
   free_pdf (&proton);

   return status;
}

/*}}}*/

#endif

ISIS_USER_SOURCE_MODULE(pizero,p,options) /*{{{*/
{
   static const char *parameter_names[] = {"norm", NULL};
   static const char *parameter_units[] = {"cm^-2 GeV^-1", NULL};
   static double default_max[]   = {1.e10};
   static double default_value[] = { 1.0};
   static double default_min[]   = { 0.0};
   static double default_relstep[] = {1.e-4};
   static unsigned int default_freeze[] = {0};
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
   p->default_relstep = default_relstep;
   p->default_freeze = default_freeze;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   pizero_init_client_data ();
   return 0;
}

/*}}}*/

