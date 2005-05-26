/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Feb 2003
 *
 * Neutral pion production cross-sections from
 *   Blattnig et al, NASA Technical Report, NASA/TP-2000-210640
 *      Equations 23 and 24, page 10.  See also Fig 105 on page 76.
 *      (Note that the curves in the figure are labeled
 *      in the wrong order!)
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
#include "pion.h"

#define SQR_PION_MASS_ENERGY   (PI_ZERO_REST_ENERGY * PI_ZERO_REST_ENERGY)
#define PION_MASS_FACTOR       (4*PROTON_REST_ENERGY * PROTON_REST_ENERGY \
                               - SQR_PION_MASS_ENERGY)

#if 0
#include "pion_xsec_nasa.inc"
#else
#include "pion_xsec_dermer.inc"
#endif

static double pion_integrand (double epion, void *x) /*{{{*/
{
   Pion_Type *p = (Pion_Type *)x;
   double val;

   (void) npxsec (epion, p->proton_kinetic_energy, &val);
   val /= sqrt (epion * epion - SQR_PION_MASS_ENERGY);

   return val;
}

/*}}}*/

static int pion_min_energy (double photon_energy, double *epion) /*{{{*/
{
   *epion = photon_energy + 0.25 * SQR_PION_MASS_ENERGY / photon_energy;
   return 0;
}

/*}}}*/

static int pion_max_energy (Particle_Type *proton, double *pion_kinetic_energy) /*{{{*/
{
   double pc_max = (*proton->momentum_max)(proton);
   double ep_max = sqrt (pc_max*pc_max + PROTON_REST_ENERGY*PROTON_REST_ENERGY);
   double root_s = 2 * ep_max;
   double epion;

   /* From Blattnig et al (2000), Appendix A */
   epion = 0.5 * (root_s - PION_MASS_FACTOR / root_s);
   *pion_kinetic_energy = epion - PI_ZERO_REST_ENERGY;

   return 0;
}

/*}}}*/

static int integral_over_pion_energies (Pion_Type *p, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double epi_min, epi_max;
   size_t limit;
   int status;

   *val = 0.0;

   if (-1 == pion_min_energy (p->photon_energy, &epi_min))
     return -1;
   if (-1 == pion_max_energy (p->protons, &epi_max))
     return -1;

   if (epi_min >= epi_max)
     return 0;

   f.function = &pion_integrand;
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

static double proton_integrand (double proton_kinetic_energy, void *x) /*{{{*/
{
   Pion_Type *p = (Pion_Type *)x;
   Particle_Type *protons = p->protons;
   double pion_flux, np, pc;

   p->proton_kinetic_energy = proton_kinetic_energy;
   pc = sqrt (proton_kinetic_energy
              * (proton_kinetic_energy + 2*PROTON_REST_ENERGY));

   (void)(*protons->spectrum) (protons, pc, &np);
   if (-1 == integral_over_pion_energies (p, &pion_flux))
     return 0.0;

   return np * pion_flux;
}
/*}}}*/

static double lg_proton_integrand (double t, void *x)
{
   /* change of variables yields better result */
   return exp(t) * proton_integrand (exp(t), x);
}

static int proton_min_energy (Particle_Type *proton, double *ep_min) /*{{{*/
{
   double pc_min = (*proton->momentum_min)(proton);
   *ep_min = sqrt (pc_min * pc_min + PROTON_REST_ENERGY*PROTON_REST_ENERGY);
   return 0;
}

/*}}}*/

static int proton_max_energy (Particle_Type *proton, double *ep_max) /*{{{*/
{
   double pc_max = (*proton->momentum_max)(proton);
   *ep_max = sqrt (pc_max * pc_max + PROTON_REST_ENERGY*PROTON_REST_ENERGY);
   return 0;
}

/*}}}*/

static int integral_over_proton_energies (Pion_Type *p, double *val) /*{{{*/
{
   Particle_Type *proton = p->protons;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double ep_min, ep_max;
   size_t limit;
   int status;

   if (-1 == proton_min_energy (proton, &ep_min))
     return -1;
   if (-1 == proton_max_energy (proton, &ep_max))
     return -1;

   f.function = &lg_proton_integrand;
   f.params = p;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, log(ep_min), log(ep_max), epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS61,
                                 work, val, &abserr);

   *val *= 2 * (4 * M_PI * p->density);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** %s\n", gsl_strerror (status));

   return 0;
}

/*}}}*/

int pi_calc_pion_decay (Pion_Type *p, double photon_energy, double *val)
{
   p->photon_energy = photon_energy * GSL_CONST_CGSM_ELECTRON_VOLT;
   return integral_over_proton_energies (p, val);
}

#if 1
#define MILLIBARNS_PER_GEV  (1.e-27/GEV)
int main (void) /*{{{*/
{
#if 1
   double lg_epi, lg_epi_min, lg_epi_max, dlg_epi, ep;

   ep = 1.0e2;

   lg_epi_min = log10(ep)-3.0;
   lg_epi_max = log10(ep);
   dlg_epi = 0.1;

   for (lg_epi = lg_epi_min; lg_epi < lg_epi_max; lg_epi += dlg_epi)
     {
        double xsec, epi;

        epi = pow(10.0,lg_epi);

        if (-1 == npxsec (epi*GEV, ep*GEV, &xsec))
          exit(1);

        xsec /= MILLIBARNS_PER_GEV;
        fprintf (stdout, "%12.4e %12.4e %12.4e\n", ep, epi, xsec);
     }
#else
   double xsec;
   if (-1 == npxsec (0.1*GEV, 0.5*GEV, &xsec))
     exit(1);
   xsec /= MILLIBARNS_PER_GEV;
   fprintf (stdout, "dsigma/dE = %g\n", xsec);
#endif

   return 0;
}

/*}}}*/
#endif
