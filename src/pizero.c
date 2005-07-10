/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Feb 2003
 *
 * Neutral pion production cross-sections from
 *   Blattnig et al, NASA Technical Report, NASA/TP-2000-210640
 *   Aharonian and Atoyan, 2000, A&A, 362, 937
 *
 * See also
 *   Domingo-Santamaria & Torres, 2005, astro-ph/0506240
 *   Dermer, C.D., 1986, A&A, 157, 223.
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

#define SQR_PROTON_REST_ENERGY  (PROTON_REST_ENERGY * PROTON_REST_ENERGY)
#define SQR_PIZERO_REST_ENERGY  (PIZERO_REST_ENERGY * PIZERO_REST_ENERGY)
#define PIZERO_MASS_FACTOR      (SQR_PIZERO_REST_ENERGY - 4*SQR_PROTON_REST_ENERGY)
#define SQR_PIZERO_MASS_FACTOR  (PIZERO_MASS_FACTOR * PIZERO_MASS_FACTOR)

static double pizero_total_xsec (double proton_kinetic) /*{{{*/
{
   double sigma;

   if (proton_kinetic > GEV)
     sigma = 30.0 * (0.95 + 0.06 * log (proton_kinetic/GEV));
   else return 0.0;

   return sigma * MILLIBARN;
}

/*}}}*/

static double pizero_differential_xsec (double proton_kinetic, double pizero_kinetic) /*{{{*/
{
   double a, sigma;

   /* Blattnig, et al, 2000 Lab-frame differential cross-section
    * parameterization (their equation 32)
    */
   a = (-5.8 - 1.82/pow(proton_kinetic, 0.4)
        + 13.5/pow(pizero_kinetic, 0.2) - 4.5/pow(pizero_kinetic, 0.4));

   /* d(sigma)/dE_pizero */
   sigma = exp(a) * (MILLIBARN / GEV);

   return sigma;
}

/*}}}*/

static double proton_integrand (double pc, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   Particle_Type *proton = p->protons;
   double r_proton, gamma, beta, v, np, xsec;
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
static double delta_function_approximation (Pizero_Type *p) /*{{{*/
{
   Particle_Type *protons = p->protons;
   double eproton_delta, proton_pc, proton_kinetic;
   double np, gamma, beta, sigma, v, val;

   /* kappa = mean fraction of proton kinetic energy transferred
    *         to the secondary meson per collision. */
   double kappa = 0.17;

   eproton_delta = PROTON_REST_ENERGY + p->energy /kappa;
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

/*}}}*/

static int Threshold_Flag;

static int integral_over_proton_momenta (Pizero_Type *p, double *val) /*{{{*/
{
   Particle_Type *protons = p->protons;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double pc_min, pc_max, x0, eproton_thresh, pc_thresh;
   size_t limit;
   int status;

   *val = 0.0;

   if (p->energy > 100.0 * GEV)
     {
        *val = delta_function_approximation (p);
        return 0;
     }

   if (p->energy < 100.0 * MEV)
     return 0;

   pc_max = (*protons->momentum_max)(protons);
   pc_min = (*protons->momentum_min)(protons);

   /* threshold proton momentum to produce the given pi-zero */
   x0 = (p->energy * p->energy - 0.5*PIZERO_MASS_FACTOR) / PROTON_REST_ENERGY;
   eproton_thresh = (x0 * (1.0 + sqrt (1.0 - SQR_PIZERO_MASS_FACTOR / (x0*x0)))
                     - PROTON_REST_ENERGY);
   pc_thresh = sqrt (eproton_thresh*eproton_thresh - SQR_PROTON_REST_ENERGY);

   if (pc_thresh < pc_min)
     {
        Threshold_Flag = 1;
        return 0;
     }

   f.function = &proton_integrand;
   f.params = p;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, pc_thresh, pc_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** pizero: proton integral:  %s\n", gsl_strerror (status));

   return 0;
}

/*}}}*/

static double pizero_integrand (double epizero, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   double pizero_pc, y, s;

   p->energy = epizero;

   /* FIXME!!!  far better to pre-compute pizero distribution,
    * and then interpolate on the saved table here.
    */
   if (-1 == integral_over_proton_momenta (p, &y))
     return 0.0;

   pizero_pc = sqrt (epizero * epizero - SQR_PIZERO_REST_ENERGY);
   s = y / pizero_pc;

   return s;
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
   double pc_max, ep_max, root_s;

   pc_max = (*protons->momentum_max)(protons);
   ep_max = sqrt (pc_max*pc_max + SQR_PROTON_REST_ENERGY);
   root_s = sqrt (2*PROTON_REST_ENERGY*(PROTON_REST_ENERGY + ep_max));

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

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, epi_min, epi_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);
   /* two photons per pion */
   *val *= 2.0;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** pizero: pizero integral: %s\n", gsl_strerror (status));

   return 0;
}

/*}}}*/

int pizero_decay (void *x, double photon_energy, double *emissivity)
{
   Pizero_Type *p = (Pizero_Type *)x;
   int status;

   *emissivity = 0.0;
   photon_energy *= GSL_CONST_CGSM_ELECTRON_VOLT;

   Threshold_Flag = 0;
   status = integral_over_pizero_energies (p, photon_energy, emissivity);
   if (Threshold_Flag)
     {
        fprintf (stderr, "WARNING:  Underestimated pizero rate for %g GeV photons (proton distribution doesn't extend to low enough momenta)\n",
                 photon_energy / GEV);
     }

   return status;
}

