/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "_nonthermal.h"

static double momentum (double gamma, double mass) /*{{{*/
{
   return mass * C_SQUARED * sqrt ((gamma + 1.0) * (gamma - 1.0));
}

/*}}}*/

static double particle_momentum_min (Particle_Type *pt) /*{{{*/
{
   return momentum (GAMMA_MIN_DEFAULT, pt->mass);
}

/*}}}*/

static double particle_momentum_max (Particle_Type *pt) /*{{{*/
{
   double e_cutoff, mc2, r, f=1.e-8;

   if (pt->cutoff_energy == 0)
     return momentum (GAMMA_MAX_DEFAULT, pt->mass);

   /* Choose max momentum high enough so that the exponential cutoff
    * represents a factor of 'f' decline in particle density
    */

   mc2 = pt->mass * C_SQUARED;
   e_cutoff = pt->cutoff_energy * TEV;
   r = (GEV - e_cutoff * log(f)) / mc2;

   return mc2 * sqrt (r*r - 1.0);
}

/*}}}*/

static int particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/
{
   double mc2, e_k, e0, pcomc2, x, g, f;

   if (pt == NULL || ne == NULL)
     return -1;

   *ne = 0.0;

   mc2 = pt->mass * C_SQUARED;
   e0  = pt->cutoff_energy * TEV;

   pcomc2 = pc/mc2;
   e_k = mc2 * (sqrt (1.0 + pcomc2 * pcomc2) - 1.0);

   x = pc / GEV;
   g = pt->index;

   /* no curvature below 1 GeV */
   if ((pt->curvature != 0.0) && (x > 1.0))
     g += pt->curvature * log10(x);

   /* dn/d(Pc) (norm factored out) */
   f = pow (x, g) * exp ((GEV-e_k)/e0);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

int init_particle_spectrum (Particle_Type *pt) /*{{{*/
{
   if (pt == NULL)
     return -1;

   pt->spectrum = &particle_spectrum;
   pt->momentum_min = &particle_momentum_min;
   pt->momentum_max = &particle_momentum_max;

   return 0;
}

/*}}}*/

