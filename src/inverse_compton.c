/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 *
 * Implementation based on
 *   Baring, M.G., et al, 1999, ApJ, 513, 311
 *   Blumenthal, G.R. & Gould, R.J., 1970, Rev. Mod. Phys., 43, 237
 */

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

static int sigma_klein_nishina (double energy_initial_photon, /*{{{*/
                                double gamma_electron,
                                double energy_final_photon, double *sigma)
{
   double q, g, eg;

   *sigma = 0.0;

   /* energies in units of me*c^2 */

   if ((energy_initial_photon <= 0.0)
       || (gamma_electron <= 0.0)
       || (energy_final_photon <= 0.0))
     return -1;

   if (energy_initial_photon > energy_final_photon)
     return 0;

   eg = energy_initial_photon * gamma_electron;
   g = 4.0 * eg;
   q = (energy_final_photon
        / (g * (gamma_electron - energy_final_photon)));

   if (q < 0.0 || 1.0 < q)
     return 0;

   *sigma = (2.0 * q * log(q)
             + (1.0 - q) * (1.0 + q * (2.0 + 0.5*g*g*q/(1.0 + g*q))));

   *sigma *= (KLEIN_NISHINA_COEF / (eg * gamma_electron));

   return 0;
}

/*}}}*/

static double incident_photons_integrand (double energy_initial_photon, void *p) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)p;
   double sigma_kn, num_photons;

   (void) sigma_klein_nishina (energy_initial_photon,
                               ic->gamma_electron,
                               ic->energy_final_photon,
                               &sigma_kn);

   (void) (*ic->incident_photons) (energy_initial_photon, &num_photons);

   return num_photons * sigma_kn;
}

/*}}}*/

int ic_integral_over_incident_photons (Inverse_Compton_Type *ic, /*{{{*/
                                       double *val)
{
   double epsabs, epsrel, abserr, e_max;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   size_t limit;

   f.function = &incident_photons_integrand;
   f.params = ic;
   epsabs = 0.0;
   epsrel = 1.e-10;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   e_max = ic->incident_photon_max_energy;
   e_max *= (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);

   gsl_error_handler = gsl_set_error_handler_off ();

   (void) gsl_integration_qag (&f, 0.0, e_max, epsabs, epsrel, limit,
                               GSL_INTEG_GAUSS31, work,
                               val, &abserr);

   gsl_set_error_handler (gsl_error_handler);

   gsl_integration_workspace_free (work);

   return 0;
}

/*}}}*/

static double electrons_integrand (double gamma, void *p) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)p;
   Particle_Type *elec = ic->electrons;
   double ne_gamma, x;

   ic->gamma_electron = gamma;

   (void) (*elec->spectrum) (elec, gamma, &ne_gamma);

   if (ic->interpolate == 0)
     (void) ic_integral_over_incident_photons (ic, &x);
   else
     (void) ic_interp_photon_integral (ic, &x);

   return ne_gamma * x;
}

/*}}}*/

static int integral_over_electrons (Inverse_Compton_Type *ic, /*{{{*/
                                    double *val)
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double gamma_min, gamma_max;
   size_t limit;
   int status;

   f.function = &electrons_integrand;
   f.params = ic;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;
   
   gamma_min = (*ic->electrons->gamma_min) (ic->electrons);
   gamma_max = (*ic->electrons->gamma_max) (ic->electrons);   

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, gamma_min, gamma_max, epsabs, epsrel,
                                 limit, GSL_INTEG_GAUSS31, work,
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

int ic_calc_inverse_compton (void *vic, /*{{{*/
                             double energy_final_photon,
                             double *emissivity)
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)vic;
   double integral = 0.0;

   if (ic == NULL)
     return -1;

   /* eV -> dimensionless */
   energy_final_photon *= (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);
   ic->energy_final_photon = energy_final_photon;

   if (-1 == integral_over_electrons (ic, &integral))
     return -1;

   *emissivity = GSL_CONST_CGSM_SPEED_OF_LIGHT * integral;

   return 0;
}

/*}}}*/

