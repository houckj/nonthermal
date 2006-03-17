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

static double incident_photons_integrand (double q, void *p) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)p;
   double omega = ic->energy_final_photon;
   double gamma = ic->electron_gamma;
   double x, omega_i, gq, f_kn, num_photons, val;

   /* changed integration variable */
   q = exp(q);

   x = omega/gamma;
   gq = x / (1.0 - x);
   omega_i = gq / (4 * gamma) / q;
   f_kn = 2 * q * log(q) + (1.0 - q) * (1.0 + 2*q + 0.5*gq*gq/(1.0 + gq));

   (void) (*ic->incident_photons) (omega_i, &num_photons);
   if (num_photons == 0.0)
     return 0.0;

   val = num_photons * f_kn / q;

   /* changed integration variable */
   val *= q;

   return val;
}

/*}}}*/

int ic_integral_over_incident_photons (Inverse_Compton_Type *ic, /*{{{*/
                                       double *val)
{
   double epsabs, epsrel, abserr;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   size_t limit;
   double x, h, min_omega_i, max_omega_i, q_min, q_max;
   double gamma = ic->electron_gamma;
   double omega = ic->energy_final_photon;

   *val = 0.0;

   x = omega / gamma;
   if (x > 1.0)
     return 0;

   h = 0.25 * x / (1.0 - x) / gamma;
   max_omega_i = ic->incident_photon_max_energy;
   max_omega_i *= (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);
   q_min = h / max_omega_i;

   f.function = &incident_photons_integrand;
   f.params = ic;
   epsabs = 0.0;
   epsrel = 1.e-12;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   (void) gsl_integration_qag (&f, log(q_min), 0.0, epsabs, epsrel, limit,
                               GSL_INTEG_GAUSS31, work,
                               val, &abserr);

   *val *= KLEIN_NISHINA_COEF / gamma / gamma;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   return 0;
}

/*}}}*/

static double Exponent = 3.0;

static double electrons_integrand (double pc, void *pt) /*{{{*/
{
   Inverse_Compton_Type *ic = (Inverse_Compton_Type *)pt;
   Particle_Type *elec = ic->electrons;
   double pcomc2, x, ne;
   
   /* changed integration variable */
   double arg = pc;
   pc = pow (pc, Exponent);

   pcomc2 = pc / ELECTRON_REST_ENERGY;
   ic->electron_gamma = sqrt (1.0 + pcomc2*pcomc2);

   if (ic->interpolate == 0)
     (void) ic_integral_over_incident_photons (ic, &x);
   else
     (void) ic_interp_photon_integral (ic, &x);

   (void) (*elec->spectrum) (elec, pc, &ne);

   /* changed integration variable */
   x *= Exponent * pc / arg;
   
   return ne * x;
}

/*}}}*/

static int ic_integral_over_electrons (Inverse_Compton_Type *ic, /*{{{*/
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
   epsrel = 1.e-12;
   limit = MAX_QAG_SUBINTERVALS;

   pc_max = (*ic->electrons->momentum_max) (ic->electrons);
   pc_min = (*ic->electrons->momentum_min) (ic->electrons);

   {
      double pcm, gamma_min, max_omega_i;
      double omega = ic->energy_final_photon;
#if 0
      /* lowest possible threshold occurs for largest omega_i */
      max_omega_i = (incident_photon_max_energy ()
                     * GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY);
#else
      /* FIXME:  Omega0 from lookup table */
      max_omega_i = 1.48981689186384e-08;
#endif

      gamma_min = 0.5 * (omega + sqrt (omega * (omega + 1.0/max_omega_i)));
      pcm = ELECTRON_REST_ENERGY * sqrt((gamma_min + 1.0) * (gamma_min - 1.0));
      if (pcm > pc_min)
        pc_min = pcm;
   }

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qagiu (&f, pow(pc_min, 1.0/Exponent), epsabs, epsrel, limit, work,
                                   val, &abserr);

   *val /= ELECTRON_REST_ENERGY;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     {
        if (status != GSL_EROUND)
          fprintf (stderr, "*** invc:  %s\n", gsl_strerror (status));
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

   if (-1 == ic_integral_over_electrons (ic, &integral))
     return -1;

   *emissivity = GSL_CONST_CGSM_SPEED_OF_LIGHT * integral;

   return 0;
}

/*}}}*/

