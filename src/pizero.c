/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Feb 2003
 *
 * Neutral pion production cross-sections from
 *   Blattnig et al, NASA Technical Report, NASA/TP-2000-210640
 *
 * See also
 * Dermer, C.D., 1986, A&A, 157, 223.
 *
*/

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "_nonthermal.h"
#include "pizero.h"

#define SQR_PROTON_REST_ENERGY (PROTON_REST_ENERGY * PROTON_REST_ENERGY)

#define SQR_PIZERO_REST_ENERGY   (PIZERO_REST_ENERGY * PIZERO_REST_ENERGY)
#define PIZERO_MASS_FACTOR       (SQR_PIZERO_REST_ENERGY - 4*SQR_PROTON_REST_ENERGY)

#define SQR_PIZERO_MASS_FACTOR   (PIZERO_MASS_FACTOR * PIZERO_MASS_FACTOR)

static double pizero_total_xsec (double proton_kinetic)
{
   double sigma;

   /* Aharonian and Atoyan, 2000, A&A, 362, 937 */
   if (proton_kinetic < GEV)
     return 0.0;
   sigma = 30.0 * (0.95 + 0.06 * log (proton_kinetic/GEV));

   return sigma * MILLIBARN;
}

static double pizero_differential_xsec (double proton_kinetic, double pizero_kinetic)
{
   double a, sigma;

   /* Blattnig, et al, 2000 Lab-frame differential cross-section
    * parameterization (their equation 32)
    */
   a = (-5.8
        - 1.82/pow(proton_kinetic, 0.4)
        + 13.5/pow(pizero_kinetic, 0.2)
        - 4.5/pow(pizero_kinetic, 0.4));

   /* d(sigma)/dE_pizero */
   sigma = exp(a) * (MILLIBARN / GEV);

   return sigma;
}

static double proton_integrand (double pc, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   Particle_Type *proton = p->protons;
   double r_proton, gamma, gm1, beta, v, np, xsec;
   double pizero_kinetic, proton_kinetic, result;

   (void)(*proton->spectrum) (proton, pc, &np);

   r_proton = pc / PROTON_REST_ENERGY;
   gamma = sqrt (r_proton * r_proton + 1.0);
   beta = sqrt ((gamma + 1.0)*(gamma - 1.0))/gamma;
   proton_kinetic = (gamma - 1.0) * PROTON_REST_ENERGY;

   pizero_kinetic = p->energy - PIZERO_REST_ENERGY;

   /* cross-section per unit proton energy */
   xsec = pizero_differential_xsec (pizero_kinetic, proton_kinetic);

   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;
   result = np * v * xsec;

   /* We're integrating over momentum (pc, really), not energy, so:
    * dE = d(pc) (pc/E) = d(pc) * beta */
   result *= beta;

   return result;
}
/*}}}*/

/* delta-function approximation: see Aharonian and Atoyan (2000) */
static double delta_function_approximation (Pizero_Type *p)
{
   Particle_Type *protons = p->protons;
   double eproton_delta, proton_pc, proton_kinetic;
   double np, gamma, beta, sigma, v, val;

   /* kappa = mean fraction of proton kinetic energy transferred
    *         to the secondary meson per collision. */
   double kappa = 0.17;

   eproton_delta = (p->energy /kappa + PROTON_REST_ENERGY);
   proton_pc = sqrt (eproton_delta * eproton_delta - SQR_PROTON_REST_ENERGY);
   (void)(*protons->spectrum)(protons, proton_pc, &np);

   gamma = eproton_delta /PROTON_REST_ENERGY;
   beta = sqrt ((gamma + 1.0)*(gamma - 1.0))/gamma;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   proton_kinetic = eproton_delta - PROTON_REST_ENERGY;
   sigma = pizero_total_xsec (proton_kinetic);

   val = np * v * sigma /kappa;

   return val;
}

static int integral_over_proton_momenta (Pizero_Type *p, double *val) /*{{{*/
{
   Particle_Type *protons = p->protons;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double x2, epizero2, eproton_thresh2, pc_thresh, pc_min, pc_max;
   size_t limit;
   int status;
   
   *val = 0.0;
   
   if (p->energy > 100.0 * GEV)
     {
        *val = delta_function_approximation (p);
        return 0;
     }

   pc_max = (*protons->momentum_max)(protons);
   pc_min = (*protons->momentum_min)(protons);

   /* threshold proton momentum to produce the given pi-zero */
   epizero2 = p->energy * p->energy;
   x2 = SQR_PIZERO_MASS_FACTOR / epizero2;
   eproton_thresh2 = (0.25 * SQR_PIZERO_MASS_FACTOR
                      * (1.0 + (2/x2) * (1.0 + /*+/-*/ sqrt(1.0 + x2))));
   pc_thresh = sqrt (eproton_thresh2 - SQR_PROTON_REST_ENERGY);

   if (pc_thresh < pc_min)
     {
        fprintf (stderr, "ERROR:  pizero production threshold falls below proton momentum distribution lower bound\n");
        fprintf (stderr, "        threshold pc = %g   lower bound pc = %g\n",
                 pc_thresh, pc_min);
        return 0;
     }

   f.function = &proton_integrand;
   f.params = p;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, pc_thresh, pc_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS61,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** %s\n", gsl_strerror (status));

   return 0;
}

/*}}}*/

static double pizero_integrand (double epizero, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   double pizero_pc, y;

   p->energy = epizero;

   /* FIXME!!!  far better to pre-compute pizero distribution,
    * and then interpolate on the saved table here.
    */
   if (-1 == integral_over_proton_momenta (p, &y))
     return 0.0;

   pizero_pc = sqrt (epizero * epizero - SQR_PIZERO_REST_ENERGY);

   return y / pizero_pc;
}

/*}}}*/

static int pizero_min_energy (double photon_energy, double *epizero) /*{{{*/
{
   *epizero = photon_energy + 0.25 * SQR_PIZERO_REST_ENERGY / photon_energy;

   return 0;
}

/*}}}*/

static int pizero_max_energy (Particle_Type *protons, double *epizero) /*{{{*/
{
   double pc_max = (*protons->momentum_max)(protons);
   double ep_max = sqrt (pc_max*pc_max + SQR_PROTON_REST_ENERGY);
   double root_s = 2 * ep_max;

   /* See Blattnig et al (2000), Appendix A */
   *epizero = 0.5 * (root_s + PIZERO_MASS_FACTOR / root_s);

   return 0;
}

/*}}}*/

static int integral_over_pizero_energies (Pizero_Type *p, double photon_energy, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double epi_min, epi_max;
   size_t limit;
   int status;

   *val = 0.0;

   if (-1 == pizero_min_energy (photon_energy, &epi_min))
     return -1;
   if (-1 == pizero_max_energy (p->protons, &epi_max))
     return -1;

   if (epi_min >= epi_max)
     return 0;

   /* FIXME:
    * At this point, should pre-compute pizero distribution over
    * the energy range [epi_min, epi_max]
    */

   f.function = &pizero_integrand;
   f.params = p;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, epi_min, epi_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS61,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** %s\n", gsl_strerror (status));

   return 0;
}

/*}}}*/

int pizero_decay (void *x, double photon_energy, double *emissivity)
{
   Pizero_Type *p = (Pizero_Type *)x;
   photon_energy *= GSL_CONST_CGSM_ELECTRON_VOLT;
   return integral_over_pizero_energies (p, photon_energy, emissivity);
}


