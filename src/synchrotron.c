/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 *
 * Implementation based on
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

static double angular_integrand (double t, void *p) /*{{{*/
{
   double x = *(double *)p;
   gsl_sf_result f;

   if (t == 0.0 || t == 1.0)
     return 0.0;

   (void) gsl_sf_synchrotron_1_e (x/t, &f);

   return f.val * t*t / sqrt ((1.0 + t)*(1.0 - t));
}

/*}}}*/

int syn_angular_integral (double x, double *y) /*{{{*/
{
   double epsabs, epsrel, abserr;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   size_t limit;

   *y = 0.0;

   f.function = &angular_integrand;
   f.params = &x;
   epsabs = 0.0;
   epsrel = 1.e-10;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

#if 0   
   (void) gsl_integration_qag (&f, 0.0, 1.0, epsabs, epsrel, limit,
                               GSL_INTEG_GAUSS31,
                               work, y, &abserr);
#else
   {
      double pts[] = {0.0, 1.0};
      int status;
      
      status = gsl_integration_qagp (&f, pts, 2, epsabs, epsrel, limit,
                                     work, y, &abserr);
      if (status)
        {
           if (status != GSL_EROUND)           
             fprintf (stderr, "status = %d:  %s\n",
                      status,
                      gsl_strerror (status));
        }      
   }   
#endif   

   gsl_set_error_handler (gsl_error_handler);

   gsl_integration_workspace_free (work);

   return 0;
}

/*}}}*/

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

static double synchrotron_integrand (double pc, void *pt) /*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)pt;
   Particle_Type *elec = s->electrons;
   double pcomc2, gamma2, x, y, ne;

   pcomc2 = pc / (elec->mass * C_SQUARED);
   gamma2 = 1.0 + pcomc2*pcomc2;
   /* x = (photon energy) / (critical energy) */
   x = (s->photon_energy
           / (SYNCHROTRON_CRIT_ENERGY_COEF * s->B_tot * gamma2));

   (void) eval_angular_integral (x, pt, &y);
   (void) (*elec->spectrum) (elec, pc, &ne);

   return ne * y;
}

/*}}}*/

int syn_calc_synchrotron (void *vs, double photon_energy, double *emissivity)/*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)vs;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr, integral, eph;
   double pc_min, pc_max;
   size_t limit;
   int status;

   /* eV */
   s->photon_energy = photon_energy;

   f.function = &synchrotron_integrand;
   f.params = s;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   pc_min = (*s->electrons->momentum_min) (s->electrons);
   pc_max = (*s->electrons->momentum_max) (s->electrons);

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, pc_min, pc_max, epsabs, epsrel,
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

   eph = photon_energy * GSL_CONST_CGSM_ELECTRON_VOLT;
   *emissivity = (SYNCHROTRON_COEF * s->B_tot * integral) / eph;

#if 0
   fprintf (stdout, "%15.7e  %15.7e\n", photon_energy, *emissivity);
#endif

   return 0;
}

/*}}}*/

