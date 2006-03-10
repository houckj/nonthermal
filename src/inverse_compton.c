/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 *
 * Implementation based on
 *   Baring, M.G., et al, 1999, ApJ, 513, 311
 *   Blumenthal, G.R. & Gould, R.J., 1970, Rev. Mod. Phys., 43, 237
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "_nonthermal.h"
#include "inverse_compton.h"
#include "ic_table.h"

/*{{{ physical constants */

#define KLEIN_NISHINA_COEF \
   (2 * M_PI * ELECTRON_RADIUS * ELECTRON_RADIUS)

/*}}}*/

static int sigma_klein_nishina (double omega_i, double gamma_e, /*{{{*/
                                double omega, double *sigma)
{
   double q, g, eg, qmin, x;

   *sigma = 0.0;

   /* energies in units of me*c^2 */

   if ((omega_i <= 0.0) || (gamma_e <= 0.0) || (omega <= 0.0))
     return -1;

   if (omega_i > omega)
     return 0;

   eg = omega_i * gamma_e;
   g = 4.0 * eg;
   x = omega / gamma_e;
   q = (x/g) / (1.0 - x);

   qmin = 0.25 / gamma_e / gamma_e;
   if (q < qmin || 1.0 < q)
     return 0;

   *sigma = (2.0 * q * log(q)
             + (1.0 - q) * (1.0 + (2.0/g + 0.5*x) * x / (1.0 - x)));

   *sigma *= (KLEIN_NISHINA_COEF / (eg * gamma_e));

   return 0;
}

/*}}}*/

static double incident_photons_integrand (double lg_ei, void *p) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)p;
   double sigma_kn, num_photons, omega_i, val;

   /* changed integration variable */
   omega_i = exp (lg_ei);

   (void) (*ic->incident_photons) (omega_i, &num_photons);
   if (num_photons == 0.0)
     return 0.0;

   (void) sigma_klein_nishina (omega_i,
                               ic->electron_gamma,
                               ic->energy_final_photon,
                               &sigma_kn);

   val = num_photons * sigma_kn;

   /* changed integration variable */
   return val * omega_i;
}

/*}}}*/

int ic_integral_over_incident_photons (Inverse_Compton_Type *ic, /*{{{*/
                                       double *val)
{
   double epsabs, epsrel, abserr, e_min, e_max;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   size_t limit;

   f.function = &incident_photons_integrand;
   f.params = ic;
   epsabs = 0.0;
   epsrel = 1.e-12;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   e_max = ic->incident_photon_max_energy;
   e_max *= (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);
   e_min = ic->incident_photon_min_energy;
   e_min *= (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);

   gsl_error_handler = gsl_set_error_handler_off ();

   (void) gsl_integration_qag (&f, log(e_min), log(e_max), epsabs, epsrel, limit,
                               GSL_INTEG_GAUSS31, work,
                               val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   return 0;
}

/*}}}*/

static double electrons_integrand (double pc, void *pt) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)pt;
   Particle_Type *elec = ic->electrons;
   double pcomc2, x, ne;

   pcomc2 = pc / ELECTRON_REST_ENERGY;
   ic->electron_gamma = sqrt (1.0 + pcomc2*pcomc2);

   if (ic->interpolate == 0)
     (void) ic_integral_over_incident_photons (ic, &x);
   else
     (void) ic_interp_photon_integral (ic, &x);

   (void) (*elec->spectrum) (elec, pc, &ne);

   return ne * (x / ELECTRON_REST_ENERGY);
}

/*}}}*/

static int integral_over_electrons (Inverse_Compton_Type *ic, /*{{{*/
                                    double *val)
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double pc_min, pc_max;
   size_t limit;
   int status;

   f.function = &electrons_integrand;
   f.params = ic;
   epsabs = 0.0;
   epsrel = 1.e-11;
   limit = MAX_QAG_SUBINTERVALS;

   pc_min = (*ic->electrons->momentum_min) (ic->electrons);
   pc_max = (*ic->electrons->momentum_max) (ic->electrons);

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qagiu (&f, pc_min, epsabs, epsrel, limit, work,
                                   val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     {
        if (! ((ic->interpolate != 0) && (status == GSL_EROUND)))
          fprintf (stderr, "*** %s\n", gsl_strerror (status));
     }

   return 0;
}

/*}}}*/

int ic_calc_inverse_compton (void *vic, double omega, double *emissivity) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)vic;
   double integral = 0.0;

   if (ic == NULL)
     return -1;

   /* eV -> dimensionless */
   omega *= (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);
   ic->energy_final_photon = omega;

   if (-1 == integral_over_electrons (ic, &integral))
     return -1;

   *emissivity = GSL_CONST_CGSM_SPEED_OF_LIGHT * integral;

   return 0;
}

/*}}}*/

