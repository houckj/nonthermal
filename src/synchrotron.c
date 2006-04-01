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

#define DO_ANGULAR_INTEGRAL 1

#if DO_ANGULAR_INTEGRAL
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
   *y = R_factored(x);
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

static double synchrotron_integrand (double pc, void *pt) /*{{{*/
{
   Synchrotron_Type *s = (Synchrotron_Type *)pt;
   Particle_Type *elec = s->electrons;
   double pcomc2, gamma2, x, y, ne;
   
   /* change integration variable */
   pc = exp(pc);

   pcomc2 = pc / (elec->mass * C_SQUARED);
   gamma2 = 1.0 + pcomc2*pcomc2;

   /* x = (photon energy) / (critical energy) */
   x = (s->photon_energy
           / (SYNCHROTRON_CRIT_ENERGY_COEF * s->B_tot * gamma2));

   (void) eval_angular_integral (x, pt, &y);
   (void) (*elec->spectrum) (elec, pc, &ne);

   /* change integration variable */
   y *= pc;
   
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
   double gmin2, xmax, pc_min, pc_max;
   size_t limit;
   int status;

   /* eV */
   s->photon_energy = photon_energy;

   f.function = &synchrotron_integrand;
   f.params = s;
   epsabs = 0.0;
   epsrel = 1.e-12;
   limit = MAX_QAG_SUBINTERVALS;

#if 0
   pc_max = (*s->electrons->momentum_max) (s->electrons);
   pc_min = (*s->electrons->momentum_min) (s->electrons);
#else
   /* Coefficient chosen by trying several values at random 
    * to see which one satisfied the recurrence relation best.
    */
   pc_max = 1.e3 * (*s->electrons->momentum_max) (s->electrons);
   /* FIXME:  xmax value is set by lookup table coverage. */
   xmax = 100.0;    
   gmin2 = s->photon_energy / SYNCHROTRON_CRIT_ENERGY_COEF / s->B_tot / xmax;
   if (sqrt(gmin2) > GAMMA_MIN_DEFAULT)
     pc_min = ELECTRON_REST_ENERGY * sqrt(gmin2 - 1.0);
#endif

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, log(pc_min), log(pc_max), epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS31,
                                 work, &integral, &abserr);

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

