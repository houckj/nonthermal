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
#include "pizero_table.h"

#define SQR_PROTON_REST_ENERGY  (PROTON_REST_ENERGY * PROTON_REST_ENERGY)
#define SQR_PIZERO_REST_ENERGY  (PIZERO_REST_ENERGY * PIZERO_REST_ENERGY)
#define PIZERO_MASS_FACTOR      (SQR_PIZERO_REST_ENERGY - 4*SQR_PROTON_REST_ENERGY)
#define SQR_PIZERO_MASS_FACTOR  (PIZERO_MASS_FACTOR * PIZERO_MASS_FACTOR)

#define PIZERO_MIN_ENERGY (100.0 * MEV)
#define MIN_PHOTON_ENERGY (100.0 * MEV)

double Pizero_Approx_Min_Energy = 50.0;  /* GeV */

static double pizero_total_xsec (double proton_kinetic) /*{{{*/
{
   double sigma;

   if (proton_kinetic > GEV)
     sigma = 30.0 * (0.95 + 0.06 * log (proton_kinetic/GEV));
   else return 0.0;

   return sigma * MILLIBARN;
}

/*}}}*/

/* delta-function approximation: see Aharonian and Atoyan (2000) */
static double delta_function_approximation (Pizero_Type *p) /*{{{*/
{
   Particle_Type *protons = p->protons;
   double eproton_delta, proton_pc, proton_kinetic;
   double np, beta, sigma, v, val;

   /* kappa = mean fraction of proton kinetic energy transferred
    *         to the secondary meson per collision. */
   double kappa = 0.17;

   eproton_delta = PROTON_REST_ENERGY + p->energy /kappa;
   proton_pc = sqrt (eproton_delta * eproton_delta - SQR_PROTON_REST_ENERGY);
   (void)(*protons->spectrum)(protons, proton_pc, &np);

   beta = proton_pc / eproton_delta;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   proton_kinetic = eproton_delta - PROTON_REST_ENERGY;
   sigma = pizero_total_xsec (proton_kinetic);

   val = np * v * sigma /kappa;

   return val;
}

/*}}}*/

static double pizero_differential_xsec (double proton_kinetic, double pizero_kinetic) /*{{{*/
{
   double a, sigma, xpi;

   /* Blattnig, et al, 2000 Lab-frame differential cross-section
    * parameterization (their equation 32)
    */

   proton_kinetic /= GEV;
   pizero_kinetic /= GEV;

   xpi = pow(pizero_kinetic, 0.2);
   a = -5.8 - 1.82/pow(proton_kinetic, 0.4) + (13.5 - 4.5/xpi)/xpi;

   /* d(sigma)/dE_pizero */
   sigma = exp(a) * (MILLIBARN / GEV);

   return sigma;
}

/*}}}*/

static double proton_integrand (double pc, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   Particle_Type *proton = p->protons;
   double e_proton, beta, v, np, xsec;
   double pizero_kinetic, proton_kinetic, result;

   (void)(*proton->spectrum) (proton, pc, &np);

   e_proton = sqrt (pc*pc + SQR_PROTON_REST_ENERGY);
   proton_kinetic = e_proton - PROTON_REST_ENERGY;
   pizero_kinetic = p->energy - PIZERO_REST_ENERGY;

   /* cross-section per unit proton energy */
   xsec = pizero_differential_xsec (pizero_kinetic, proton_kinetic);

   beta = pc / e_proton;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;
   result = np * v * xsec;

   /* We're integrating over momentum (pc, really), not energy, so:
    * dE = d(pc) (pc/E) = d(pc) * beta */
   result *= beta;

   return result;
}
/*}}}*/

static double proton_momentum_threshold (double e_pizero) /*{{{*/
{
   double eproton_thresh, pc;

   /* threshold proton momentum to produce the given pi-zero */

   eproton_thresh = 2*e_pizero + PROTON_REST_ENERGY
     + 0.5 * SQR_PIZERO_REST_ENERGY / PROTON_REST_ENERGY;
   pc = sqrt (eproton_thresh*eproton_thresh - SQR_PROTON_REST_ENERGY);

   return pc;
}

/*}}}*/

static struct
{
   int flag;
   double max;
}
Threshold;

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

   if (p->energy < PIZERO_MIN_ENERGY)
     return 0;

   if (p->energy > Pizero_Approx_Min_Energy * GEV)
     {
        *val = delta_function_approximation (p);
        return 0;
     }

   pc_max = (*protons->momentum_max)(protons);
   pc_min = (*protons->momentum_min)(protons);

   pc_thresh = proton_momentum_threshold (p->energy);
   if (pc_thresh < pc_min)
     {
        if ((Threshold.flag == 0)
            || (pc_thresh > Threshold.max))
          {
             Threshold.max = pc_thresh;
          }
        Threshold.flag = 1;
        pc_thresh = pc_min;
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

static int _pizero_integrand (Pizero_Type *p, double e_pizero, double *s) /*{{{*/
{
   double pc_pizero, q;

   p->energy = e_pizero;

   if (-1 == integral_over_proton_momenta (p, &q))
     return -1;

   pc_pizero = sqrt (e_pizero * e_pizero - SQR_PIZERO_REST_ENERGY);
   *s = q / pc_pizero;

   return 0;
}

/*}}}*/

static int pizero_build_table (Pizero_Type *p, double epi_min, double epi_max) /*{{{*/
{
   double lg_min, lg_max, dlg, lg;
   double *x, *y;
   unsigned int i, n = PIZERO_TABLE_SIZE;
   int status;

   if (NULL == (x = malloc (2*n * sizeof(double))))
     return -1;
   y = x + n;

   /* FIXME -- an adaptive grid might be better ... */

   lg_min = log10(epi_min);
   lg_max = log10(epi_max);
   dlg = (lg_max - lg_min) / (n - 1);

   for (i = 0; i < n; i++)
     {
        double xx, yy;
        lg = lg_min + dlg * i;
        xx = pow(10.0, lg);
        if (-1 == _pizero_integrand (p, xx, &yy))
          {
             free(x);
             return -1;
          }
        x[i] = lg;
        y[i] = log10(yy);
     }

   status = pizero_spline_table (p->client_data, x, y, n);
   free(x);

   return status;
}

/*}}}*/

static double pizero_integrand (double e_pizero, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   double s;
   int status;

   if (p->interpolate)
     {
        status = pizero_interp_pizero_integral (p->client_data, e_pizero, &s);
     }
   else
     status = _pizero_integrand (p, e_pizero, &s);

   return status ? 0.0 : s;
}

/*}}}*/

static int pizero_min_energy (double photon_energy, double *e_pizero) /*{{{*/
{
   *e_pizero = photon_energy + 0.25 * SQR_PIZERO_REST_ENERGY / photon_energy;

   return 0;
}

/*}}}*/

static int pizero_max_energy (Particle_Type *protons, double *e_pizero) /*{{{*/
{
   double pc_max, ep_max, root_s;

   pc_max = (*protons->momentum_max)(protons);
   ep_max = sqrt (pc_max*pc_max + SQR_PROTON_REST_ENERGY);
   root_s = sqrt (2*PROTON_REST_ENERGY*(PROTON_REST_ENERGY + ep_max));

   /* See Blattnig et al (2000), Appendix A */
   *e_pizero = 0.5 * (root_s + PIZERO_MASS_FACTOR / root_s);

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

   /* FIXME(?) can re-use table as long as proton distribution is fixed */
   if (p->interpolate == 1)
     {
        double min_epi_min;
        if (-1 == pizero_min_energy (MIN_PHOTON_ENERGY, &min_epi_min))
          return -1;
        if (-1 == pizero_build_table (p, min_epi_min, epi_max))
          return -1;
        p->interpolate = 2;
     }

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

   /* two photons per pion
    * (but we only ever see one of them...)*/
#if 0
   *val *= 2.0;
#endif

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

   Threshold.flag = 0;
   status = integral_over_pizero_energies (p, photon_energy, emissivity);

   if (Threshold.flag)
     {
#if 0
        fprintf (stderr, "WARNING: pizero: underestimated photon rate < %g GeV;  needed %g GeV proton threshold\n",
                 photon_energy / GEV, Threshold.max / GEV);
#endif
     }

   return status;
}

