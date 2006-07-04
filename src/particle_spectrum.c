/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "_nonthermal.h"

extern double Min_Curvature_Pc;
double Min_Curvature_Pc = 1.0;  /* GeV */

static double momentum (double gamma, double mass) /*{{{*/
{
   return mass * C_SQUARED * gamma * sqrt ((1.0 + 1.0/gamma) * (1.0 - 1.0/gamma));
}

/*}}}*/

static double particle_momentum_min (Particle_Type *pt) /*{{{*/
{
   return momentum (GAMMA_MIN_DEFAULT, pt->mass);
}

/*}}}*/

static double fixed_particle_momentum_max (Particle_Type *pt) /*{{{*/
{
   return momentum (GAMMA_MAX_DEFAULT, pt->mass);
}

/*}}}*/

static double particle_momentum_max (Particle_Type *pt) /*{{{*/
{
   double e_cutoff, mc2, r, p_max, f=1.e-8;
   double cutoff_energy;

   cutoff_energy = pt->params[2];

   if (cutoff_energy == 0)
     return fixed_particle_momentum_max (pt);

   /* Choose max momentum high enough so that the exponential cutoff
    * represents a factor of 'f' decline in particle density
    */

   mc2 = pt->mass * C_SQUARED;
   e_cutoff = cutoff_energy * TEV;
   r = (GEV - e_cutoff * log(f)) / mc2;
   p_max = mc2 * sqrt ((r + 1.0)*(r - 1.0));

   return p_max;
}

/*}}}*/

static int pc_cutoff_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/
{
   double x, g, e0, f;

   if (pt == NULL || ne == NULL)
     return -1;

   /* params:  [index, curvature, cutoff] */
   if (pt->num_params != 3)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->params[0];
   e0  = pt->params[2] * TEV;

   /* no curvature below 1 GeV */
   if ((pt->params[1] != 0.0) && (x > Min_Curvature_Pc))
     g -= pt->params[1] * log10(x);

   /* dn/d(Pc) (norm factored out) */
   f = pow (x, -g) * exp ((GEV-pc)/e0);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

static int cbreak_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/
{
   double x, g, e0, f;

   if (pt == NULL || ne == NULL)
     return -1;

   /* params:  [index, curvature, cutoff, cooling_break_momentum] */
   if (pt->num_params != 4)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->params[0];
   e0  = pt->params[2] * TEV;

   /* no curvature below 1 GeV */
   if ((pt->params[1] != 0.0) && (x > Min_Curvature_Pc))
     g -= pt->params[1] * log10(x);

   /* dn/d(Pc) (norm factored out) */
   f = pow (x, -g) * exp ((GEV-pc)/e0);

   /* steeper above the cooling break */
   if (pt->params[3] > 0)
     {
        double pc_break = pt->params[3] * TEV;
        if (pc > pc_break) f *= pc_break / pc;
     }

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

   /* params:  [index, curvature, cutoff] */
   if (pt->num_params != 3)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->params[0];
   e0  = pt->params[2] * TEV;

   /* no curvature below 1 GeV */
   if ((pt->params[1] != 0.0) && (x > Min_Curvature_Pc))
     g -= pt->params[1] * log10(x);

   mc2 = pt->mass * C_SQUARED;
   r = pc/mc2;
   sq = sqrt (1.0 + r*r);
   /* r >> 1 so no worries about precision loss here */
   T = mc2 * (sq - 1.0);

   /* dn/d(Pc) (norm factored out) */
   f = pow (x, -g) * exp ((GEV-T)/e0);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

static double etot_particle_momentum_max (Particle_Type *pt) /*{{{*/
{
   return momentum (1.e25, pt->mass);
}

/*}}}*/

/* use etot_particle_spectrum for testing sync analytic solution */
static int etot_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
{
   double mc2, r, e, x, g, f;

   if (pt == NULL || ne == NULL)
     return -1;

   /* params:  [index] */
   if (pt->num_params != 1)
     return -1;

   *ne = 0.0;

   mc2 = pt->mass * C_SQUARED;
   r = pc/mc2;
   e = mc2 * r * sqrt (1.0 + 1.0 /r /r);
   x = e / GEV;
   g = pt->params[0];

   f = pow (x, -g);

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

   /* params:  [index, curvature] */
   if (pt->num_params != 2)
     return -1;

   *ne = 0.0;

   x = pc / GEV;
   g = pt->params[0];

   /* no curvature below 1 GeV */
   if ((pt->params[1] != 0.0) && (x > Min_Curvature_Pc))
     g -= pt->params[1] * log10(x);

   mc2 = pt->mass * C_SQUARED;
   r = pc/mc2;
   E = mc2 * sqrt (1.0 + r*r);

   /* dn/d(Pc) (norm factored out) */

   /* Dermer 1986 proton spectrum */
   f = pow (E/GEV, -g);

   if (!finite(f))
     f = 0.0;

   *ne = f;

   return 0;
}

/*}}}*/

static double mori_particle_momentum_max (Particle_Type *pt) /*{{{*/
{
   return momentum (0.25e-3*GAMMA_MAX_DEFAULT, pt->mass);
}

/*}}}*/

static int mori_particle_spectrum (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
{
   double f0 = 4*M_PI/GSL_CONST_CGSM_SPEED_OF_LIGHT;
   double mc2, r, E, f, beta;

   if (pt == NULL || ne == NULL)
     return -1;

   /* no params:  */

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

static struct Particle_Type Particle_Methods[] =
{
   PARTICLE_METHODS("default",
                    pc_cutoff_particle_spectrum,
                    particle_momentum_min,
                    particle_momentum_max),
   PARTICLE_METHODS("etot",
                    etot_particle_spectrum,
                    particle_momentum_min,
                    etot_particle_momentum_max),
   PARTICLE_METHODS("ke_cutoff",
                    ke_cutoff_particle_spectrum,
                    particle_momentum_min,
                    particle_momentum_max),
   PARTICLE_METHODS("mori",
                    mori_particle_spectrum,
                    particle_momentum_min,
                    mori_particle_momentum_max),
   PARTICLE_METHODS("dermer",
                    dermer_particle_spectrum,
                    particle_momentum_min,
                    fixed_particle_momentum_max),
   PARTICLE_METHODS("cbreak",
                    cbreak_particle_spectrum,
                    particle_momentum_min,
                    particle_momentum_max),
   NULL_PARTICLE_TYPE
};

static int init_particle_params (Particle_Type *pt, unsigned int type, /*{{{*/
                                 double *pars, unsigned int num_pars)
{
   pt->num_params = num_pars;
   pt->params = NULL;
   if (num_pars > 0)
     {
        unsigned int size = num_pars * sizeof(double);
        if (NULL == (pt->params = malloc (size)))
          return -1;
        memcpy ((char *)pt->params, pars, size);
     }
   switch (type)
     {
      case ELECTRON:
        pt->mass = GSL_CONST_CGSM_MASS_ELECTRON;
        break;
      case PROTON:
        pt->mass = GSL_CONST_CGSM_MASS_PROTON;
        break;
      default:
        fprintf (stderr, "*** unrecognized particle type: %d\n", type);
        return -1;
     }
   return 0;
}

/*}}}*/

void free_particle_spectrum (Particle_Type *pt) /*{{{*/
{
   if (pt == NULL)
     return;
   free (pt->params);
   pt->params = NULL;
}

/*}}}*/

int init_particle_spectrum (Particle_Type *pt, unsigned int type, /*{{{*/
                            char *method, SLang_Array_Type *sl_pars)
{
   Particle_Type *t = Particle_Methods;
   Particle_Type *q;
   double *params = NULL;
   unsigned int num_pars = 0;
   
   if (pt == NULL)
     return -1;
   
   if ((sl_pars != NULL) && (sl_pars->num_elements > 0))
     {
        params = (double *)sl_pars->data;
        num_pars = sl_pars->num_elements;
     }   

   if (method == NULL)
     {
        /* struct copy */
        *pt = t[0];
        return init_particle_params (pt, type, params, num_pars);
     }

   for (q = t; q->method != NULL; q++)
     {
        if (strcmp(q->method,method) == 0)
          {
             /* struct copy */
             *pt = *q;
             return init_particle_params (pt, type, params, num_pars);
          }
     }

   return -1;
}

/*}}}*/

