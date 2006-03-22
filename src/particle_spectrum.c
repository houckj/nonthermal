/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "_nonthermal.h"

extern double Min_Curvature_Pc;
double Min_Curvature_Pc = 1.0;  /* GeV */

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
   double e_cutoff, mc2, r, p_max, f=1.e-8;

   if (pt->cutoff_energy == 0)
     return momentum (GAMMA_MAX_DEFAULT, pt->mass);

   /* Choose max momentum high enough so that the exponential cutoff
    * represents a factor of 'f' decline in particle density
    */

   mc2 = pt->mass * C_SQUARED;
   e_cutoff = pt->cutoff_energy * TEV;
   r = (GEV - e_cutoff * log(f)) / mc2;
   p_max = mc2 * sqrt ((r + 1.0)*(r - 1.0));

   return p_max;
}

/*}}}*/

static int pc_cutoff_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
{
   double x, g, e0, f;

   if (pt == NULL || ne == NULL)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->index;
   e0  = pt->cutoff_energy * TEV;

   /* no curvature below 1 GeV */
   if ((pt->curvature != 0.0) && (x > Min_Curvature_Pc))
     g += pt->curvature * log10(x);

   /* dn/d(Pc) (norm factored out) */
   f = pow (x, g) * exp ((GEV-pc)/e0);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

static int ke_cutoff_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
{
   double x, g, f;
   double mc2, e0, r, sq, T;

   if (pt == NULL || ne == NULL)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->index;
   e0  = pt->cutoff_energy * TEV;

   /* no curvature below 1 GeV */
   if ((pt->curvature != 0.0) && (x > Min_Curvature_Pc))
     g += pt->curvature * log10(x);

   mc2 = pt->mass * C_SQUARED;
   r = pc/mc2;
   sq = sqrt (1.0 + r*r);
   /* r >> 1 so no worries about precision loss here */
   T = mc2 * (sq - 1.0);

   /* dn/d(Pc) (norm factored out) */
   f = pow (x, g) * exp ((GEV-T)/e0);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

static int dermer_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
{
   double x, g, f;
   double mc2, r, E;

   if (pt == NULL || ne == NULL)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->index;

   /* no curvature below 1 GeV */
   if ((pt->curvature != 0.0) && (x > Min_Curvature_Pc))
     g += pt->curvature * log10(x);

   mc2 = pt->mass * C_SQUARED;
   r = pc/mc2;
   E = mc2 * sqrt (1.0 + r*r);

   /* dn/d(Pc) (norm factored out) */

   /* Dermer 1986 proton spectrum */
   f = pow (E/GEV, g);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

static int mori_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
{
   double f0 = 4*M_PI/GSL_CONST_CGSM_SPEED_OF_LIGHT;
   double mc2, r, E, f, beta;

   if (pt == NULL || ne == NULL)
     return -1;

   *ne = 0.0;

   mc2 = pt->mass * C_SQUARED;
   r = pc/mc2;
   E = mc2 * sqrt (1.0 + r*r);
   
   beta = E/pc;

   /* dn/d(Pc) (norm factored out) */

   /* Mori 1997 proton spectrum */
   if (E >= 100.0 * GEV)
     {
        double t = 2.5 * GEV/ pc;
        f = 1.67*pow(pc/GEV, -2.7) /sqrt(1.0 + t*t);
     }
   else
     {
        f = 6.65e-6*pow(E/(100.0*GEV), -2.75);
     }

   f *= f0 / beta;

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

int init_particle_spectrum (Particle_Type *pt)
{
   if (pt == NULL)
     return -1;

   /* silence compiler complaints about unused functions */
   
   if (0) {(void) &pc_cutoff_particle_spectrum;} /* <- default */
   if (0) {(void) &ke_cutoff_particle_spectrum;}
   if (0) {(void) &mori_particle_spectrum;}
   if (0) {(void) &dermer_particle_spectrum;}

   pt->spectrum = &pc_cutoff_particle_spectrum;
   pt->momentum_min = &particle_momentum_min;
   pt->momentum_max = &particle_momentum_max;

   return 0;
}



