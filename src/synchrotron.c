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
#include <gsl/gsl_sf_synchrotron.h>

#include <gsl/gsl_spline.h>

#include "_nonthermal.h"
#include "syn_table.h"

/*{{{ physical constants */

#define ELECTRON_CHARGE_OVER_ME_C \
      ((ELECTRON_CHARGE \
       / (GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT)))

#define SYNCHROTRON_CRIT_FREQ_COEF \
      ((3 / (4 * M_PI)) * ELECTRON_CHARGE_OVER_ME_C)

#define SYNCHROTRON_CRIT_ENERGY_COEF \
      ((GSL_CONST_CGSM_PLANCKS_CONSTANT_H * SYNCHROTRON_CRIT_FREQ_COEF) \
         / GSL_CONST_CGSM_ELECTRON_VOLT)  /* eV */

#define SYNCHROTRON_COEF \
      ((M_SQRT3 / (2 * M_PI)) * GSL_CONST_NUM_FINE_STRUCTURE \
          * ELECTRON_CHARGE_OVER_ME_C)

/*}}}*/

static double eval_x (double gamma, Synchrotron_Type *s) /*{{{*/
{
   /* x = (photon energy) / (critical energy) */
   return (s->photon_energy
           / (SYNCHROTRON_CRIT_ENERGY_COEF * s->B_tot * gamma * gamma));
}

/*}}}*/

static int eval_angular_integral (double gamma, void *p, double *y) /*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)p;
   double x, yx;

   x = eval_x (gamma, s);

   if (s->interpolate)
     (void) syn_interp_angular_integral (s->client_data, x, &yx);
   else
     (void) syn_angular_integral (x, &yx);

   *y = yx;

   return 0;
}

/*}}}*/

static double synchrotron_integrand (double gamma, void *p) /*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)p;
   Particle_Type *elec = s->electrons;
   double y, ne_gamma;

   (void) eval_angular_integral (gamma, p, &y);
   (void) (*elec->spectrum) (elec, gamma, &ne_gamma);

   return ne_gamma * y;
}

/*}}}*/

int syn_calc_synchrotron (void *vs, /*{{{*/
                          double photon_energy, double *emissivity)
{
   Synchrotron_Type *s = (Synchrotron_Type *)vs;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr, integral, eomc2;
   double gamma_min, gamma_max;
   size_t limit;
   int status;

   /* eV */
   s->photon_energy = photon_energy;

   f.function = &synchrotron_integrand;
   f.params = s;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;
   
   gamma_min = (*s->electrons->gamma_min) (s->electrons);
   gamma_max = (*s->electrons->gamma_max) (s->electrons);

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, gamma_min, gamma_max, epsabs, epsrel,
                                 limit, GSL_INTEG_GAUSS31, work,
                                 &integral, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     {
        if (status != GSL_EROUND)
          fprintf (stderr, "status = %d:  %s\n",
                   status,
                   gsl_strerror (status));
     }

   /* E/(m c^2) */
   eomc2 = ((photon_energy * GSL_CONST_CGSM_ELECTRON_VOLT)
        / ELECTRON_REST_ENERGY);

   *emissivity = (SYNCHROTRON_COEF * s->B_tot * integral) / eomc2;

#if 0
   fprintf (stdout, "%15.7e  %15.7e\n", photon_energy, *emissivity);
#endif

   return 0;
}

/*}}}*/

