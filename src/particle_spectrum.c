/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "_nonthermal.h"

#define SCALE_GAMMA(g,mass) (1.0 + ((g)-1)*(GSL_CONST_CGSM_MASS_ELECTRON/(mass)))

static double particle_gamma_min (Particle_Type *p) /*{{{*/
{
   (void) p;
   return SCALE_GAMMA (GAMMA_MIN_DEFAULT, p->mass);
}

/*}}}*/

static double particle_gamma_max (Particle_Type *p) /*{{{*/
{
   double gamma_max;

   if ((p == NULL) || (p->cutoff_energy == 0))
     return SCALE_GAMMA (GAMMA_MAX_DEFAULT, p->mass);

   gamma_max = - log(1.e-8) * (p->cutoff_energy * TEV) / (p->mass * C_SQUARED);

   return gamma_max;
}

/*}}}*/

static int particle_spectrum (Particle_Type *p, double gamma, double *ne) /*{{{*/
{
   double mc2, e_k, e0, pc, x, g, f, w;

   if (p == NULL || ne == NULL)
     return -1;

   *ne = 0.0;

   if (gamma <= 1.0)
     return 0;

   mc2 = p->mass * C_SQUARED;
   e0  = p->cutoff_energy * TEV;

   e_k = (gamma - 1.0) * mc2;
   w = sqrt ((gamma + 1.0)*(gamma - 1.0));
   pc = w * mc2;

   x = pc / GEV;
   g = p->index;

   /* no curvature below 1 GeV */
   if ((p->curvature != 0.0) && (x > 1.0))
     g += p->curvature * log10(x);

   /* dn/d(Pc) */
   f = pow (x, g) * exp ((GEV-e_k)/e0);

#if 0
   {
      double beta = w / gamma;
      f /= beta; /* Sturner uses this */
   }
#endif

   /* want dn/d\gamma = dn/d(Pc) * d(Pc)/d\gamma */
   f *= gamma / w;

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

int init_particle_spectrum (Particle_Type *p) /*{{{*/
{
   if (p == NULL)
     return -1;

   p->spectrum = &particle_spectrum;
   p->gamma_min = &particle_gamma_min;
   p->gamma_max = &particle_gamma_max;

   return 0;
}

/*}}}*/

