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

/* Implementation based on
 *   Baring, M.G., et al, 1999, ApJ, 513, 311
 *   Blumenthal, G.R. & Gould, R.J., 1970, Rev. Mod. Phys., 43, 237
 */

#include "config.h"
#include <float.h>
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

#undef DO_ANGULAR_INTEGRAL
#ifdef DO_ANGULAR_INTEGRAL
static double angular_integrand (double w, void *p) /*{{{*/
{
   gsl_sf_result f;
   double s, t, x = *(double *)p;

   t = sqrt ((1.0 + w)*(1.0 - w));
   (void) gsl_sf_synchrotron_1_e (x/t, &f);

   s = f.val * t;

   return s;
}

/*}}}*/

int syn_angular_integral (double x, double *y) /*{{{*/
{
   double epsabs, epsrel, abserr;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   size_t limit;
   int status;

   *y = 0.0;

   f.function = &angular_integrand;
   f.params = &x;
   epsabs = 0.0;
   epsrel = 1.e-12;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, 0.0, 1.0, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS31,
                                 work, y, &abserr);
   if (status)
     {
        if (status != GSL_EROUND)
          fprintf (stderr, "status = %d:  %s\n",
                   status,
                   gsl_strerror (status));
     }

   gsl_set_error_handler (gsl_error_handler);

   gsl_integration_workspace_free (work);

   return 0;
}

/*}}}*/
#else  /* DO_ANGULAR_INTEGRAL */

#include <gsl/gsl_sf_hyperg.h>

static double uu (double kappa, double mu, double x) /*{{{*/
{
   return gsl_sf_hyperg_U (0.5+mu-kappa, 1+2*mu, x);
}

/*}}}*/

static double R_factored (double x) /*{{{*/
{
   double a, u1, u2, s;

   u1 = uu(0.0, 4.0/3, x) * uu( 0.0, 1.0/3, x);
   u2 = uu(0.5, 5.0/6, x) * uu(-0.5, 5.0/6, x);

   if (u1 == 0.0 && u2 == 0.0)
     return 0.0;

   if (fabs(u2) < fabs(u1))
     {
        s = u1 * (1.0 - u2/u1);
     }
   else
     {
        s = u2 * (u1/u2 - 1.0);
     }

   a = 0.5 * M_PI * exp(-x) * pow(x, 11.0/3);

   return a * s;
}

/*}}}*/

int syn_angular_integral (double x, double *y) /*{{{*/
{
   *y = R_factored (x);
   return 0;
}
/*}}}*/

#endif /* DO_ANGULAR_INTEGRAL */

static int eval_angular_integral (double x, void *pt, double *y) /*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)pt;
   double yx;

   if (s->interpolate)
     (void) syn_interp_angular_integral (s->client_data, x, &yx);
   else
     (void) syn_angular_integral (x, &yx);

   *y = yx;

   return 0;
}

/*}}}*/

static double Coef;

static double synchrotron_integrand (double t, void *pt) /*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)pt;
   Particle_Type *elec = s->electrons;
   double x, pc, y, ne, val;

   /* changed integration variable */
   x = exp(t);

   /* Note that x/Coef = 1/gamma^2 */
   pc = ELECTRON_REST_ENERGY * sqrt (Coef/x - 1.0);

   (void) eval_angular_integral (x, pt, &y);
   (void) (*elec->spectrum) (elec, pc, &ne);

   val = ne * y;

   /* changed integration variable */
   val /= sqrt(x);

   /* Including the sqrt(1-x/Coef) = sqrt(1-1/gamma^2) factor
    * is "exact" but is inconsistent with the assumptions used
    * in deriving the analytic solution and therefore leads to
    * disagreement with it at low-frequencies and
    * and especially for large B_tot.  That's because the
    * analytic solution assumes ultra-relativistic electrons
    * but also extends the gamma integral to zero.
    * Letting sqrt(1-1/gamma^2)=1 is necessary to match
    * the analytic solution to "machine precision" at all
    * frequencies.
    */
#if 0
   val /= sqrt (1.0 - x/Coef);
#endif

   return val;
}

/*}}}*/

int syn_calc_synchrotron (void *vs, double photon_energy, double *emissivity)/*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)vs;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr, integral, eph;
   double xmin, xmax;
   size_t limit;
   int status = 0;

   /* eV */
   s->photon_energy = photon_energy;
   Coef = s->photon_energy / SYNCHROTRON_CRIT_ENERGY_COEF / s->B_tot;

   f.function = &synchrotron_integrand;
   f.params = s;
   epsabs = 0;
   epsrel = 1.e-10;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;
   gsl_error_handler = gsl_set_error_handler_off ();

   /* FIXME -- x range from lookup table */
   xmax = 30.0;
   xmin = 1.e-40;

   /* Changed momentum integration variable to
    * x = (photon energy) / (critical energy)
    */
   status = gsl_integration_qag
     (&f, log(xmin), log(xmax),
      epsabs, epsrel, limit, GSL_INTEG_GAUSS31, work, &integral, &abserr);

   if (status)
     {
        /* if (status != GSL_EROUND) */
        if ((fabs(integral) > 0) && (abserr > Sync_Epsrel * fabs(integral)))
          {
             Nonthermal_Error_Type e;
             e.error_msg = gsl_strerror (status);
             e.value = integral;
             e.estimated_abserr = abserr;
             e.allowed_abserr = Sync_Epsrel * fabs(integral);
             e.allowed_epsrel = Sync_Epsrel;
             nonthermal_error_hook (&e, s->electrons, __FILE__, __LINE__);
          }
     }

   /* constant coefficient from change of integration variable */
   integral *= 0.5 * sqrt (Coef) * ELECTRON_REST_ENERGY;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   eph = photon_energy * GSL_CONST_CGSM_ELECTRON_VOLT;
   *emissivity = (SYNCHROTRON_COEF * s->B_tot * integral) / eph;

   return 0;
}

/*}}}*/

