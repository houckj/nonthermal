/* -*- mode: C; mode: fold -*- */
/*
  Copyright (C) 2002, 2003, 2004, 2005, 2006 John C. Houck 

  This file is part of the nonthermal module

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "config.h"
#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "_nonthermal.h"
#include "isis.h"

extern double Min_Curvature_Pc;
double Min_Curvature_Pc = 1.0;  /* GeV */

static double momentum (double gamma, double mass) /*{{{*/
{
   return mass * C_SQUARED * gamma * sqrt ((1.0 + 1.0/gamma) * (1.0 - 1.0/gamma));
}

/*}}}*/

static double min_momentum (Particle_Type *pt) /*{{{*/
{
   return momentum (GAMMA_MIN_DEFAULT, pt->mass);
}

/*}}}*/

static double fixed_max_momentum (Particle_Type *pt) /*{{{*/
{
   return momentum (GAMMA_MAX_DEFAULT, pt->mass);
}

/*}}}*/

static double max_momentum (Particle_Type *pt) /*{{{*/
{
   double e_cutoff, mc2, r, p_max, f=1.e-8;
   double cutoff_energy;

   cutoff_energy = pt->params[2];

   if (cutoff_energy == 0)
     return fixed_max_momentum (pt);

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

static int pdf_pc_cutoff (Particle_Type *pt, double pc, double *ne) /*{{{*/
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

static int pdf_cbreak (Particle_Type *pt, double pc, double *ne) /*{{{*/
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

static int pdf_ke_cutoff (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
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

static double etot_max_momentum (Particle_Type *pt) /*{{{*/
{
   return momentum (1.e25, pt->mass);
}

/*}}}*/

/* use pdf_etot for testing sync analytic solution */
static int pdf_etot (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
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

static int pdf_dermer (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
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

static double mori_max_momentum (Particle_Type *pt) /*{{{*/
{
   return momentum (0.25e-3*GAMMA_MAX_DEFAULT, pt->mass);
}

/*}}}*/

static int pdf_mori (Particle_Type *pt, double pc, double *ne) /*{{{*/ /*{{{*/
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

static int init_pdf_params1 (Particle_Type *pt, unsigned int type, /*{{{*/
                            double *pars, unsigned int num_pars)
{
   if (num_pars != pt->num_params)
     {
        fprintf (stderr, "*** Error:  expecting %d parameters, got %d\n",
                 pt->num_params, num_pars);
        return -1;
     }
   pt->params = NULL;
   if (pt->num_params > 0)
     {
        unsigned int size = pt->num_params * sizeof(double);
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

void free_pdf (Particle_Type *pt) /*{{{*/
{
   if (pt == NULL)
     return;
   free (pt->params);
   pt->params = NULL;
}

/*}}}*/

static struct Particle_Type Particle_Methods[] =
{
   PARTICLE_METHOD("default", 3, pdf_pc_cutoff, min_momentum, max_momentum),
   PARTICLE_METHOD("etot", 1, pdf_etot, min_momentum, etot_max_momentum),
   PARTICLE_METHOD("mori", 0, pdf_mori, min_momentum, mori_max_momentum),
   PARTICLE_METHOD("ke_cutoff", 3, pdf_ke_cutoff, min_momentum, max_momentum),
   PARTICLE_METHOD("dermer", 2, pdf_dermer, min_momentum, fixed_max_momentum),
   PARTICLE_METHOD("cbreak", 4, pdf_cbreak, min_momentum, max_momentum),
   NULL_PARTICLE_TYPE
};

static Particle_Type *User_Pdf_Methods = NULL;

int init_pdf_params (Particle_Type *pt, unsigned int type, /*{{{*/
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
        return init_pdf_params1 (pt, type, params, num_pars);
     }

   for (q = t; q->method != NULL; q++)
     {
        if (strcmp(q->method,method) == 0)
          {
             /* struct copy */
             *pt = *q;
             return init_pdf_params1 (pt, type, params, num_pars);
          }
     }

   if (User_Pdf_Methods)
     {
        for (q = User_Pdf_Methods; q != NULL; q = q->next)
          {
             if (strcmp(q->method,method) == 0)
               {
                  /* struct copy */
                  *pt = *q;
                  return init_pdf_params1 (pt, type, params, num_pars);
               }
          }
     }

   return -1;
}

/*}}}*/

static void handle_link_error (char *path, char *name) /*{{{*/
{
   const char * error = dlerror ();
   if (error == NULL)
     error = "(unknown)";

   fprintf (stderr, "Link error:  %s\n", error);

   if (name == NULL)
     name = "(null)";
   if (path == NULL)
     path  = "(null)";

   fprintf (stderr, "Link error:  failed loading %s from %s\n",
            name, path);
}

/*}}}*/

typedef int Pdf_Init_Type (Particle_Type *, char *);

Particle_Type *load_pdf (char *path, char *name, char *options) /*{{{*/
{
   Particle_Type *p;
   Pdf_Init_Type *init = NULL;
   void *handle = NULL;
   char *s;

   if (path == NULL || name == NULL)
     return NULL;

#ifndef RTLD_GLOBAL
# define RTLD_GLOBAL 0
#endif
#ifdef RTLD_NOW
# define DLOPEN_FLAG  (RTLD_NOW | RTLD_GLOBAL)
#else
# define DLOPEN_FLAG  (RTLD_LAZY | RTLD_GLOBAL)
#endif

   if (NULL == (handle = dlopen (path, DLOPEN_FLAG)))
     {
        handle_link_error (path, NULL);
        return NULL;
     }

   if (NULL == (s = isis_mkstrcat ("Pdf_", name, "_init", NULL)))
     {
        dlclose (handle);
        return NULL;
     }

   if (NULL == (init = (Pdf_Init_Type *) dlsym (handle, s)))
     {
        handle_link_error (path, s);
        dlclose (handle);
        free(s);
        return NULL;
     }

   free(s);
   s = NULL;

   if (NULL == (p = malloc (sizeof *p)))
     {
        fprintf (stderr, "*** Error: load_pdf: malloc failed\n");
        dlclose (handle);
        return NULL;
     }
   memset ((char *)p, 0, sizeof *p);

   if (-1 == (*init)(p, options))
     {
        free(p);
        dlclose (handle);
        return NULL;
     }

   return p;
}

/*}}}*/

int append_pdf (Particle_Type *pt) /*{{{*/
{
   Particle_Type *p = User_Pdf_Methods;
   Particle_Type *x;

   pt->next = NULL;

   for (x = p; x != NULL; x = x->next)
     {
        /* replace method if name already exists */
        if (0 == strcmp (x->method, pt->method))
          {
             Particle_Type *n = x->next;
             free (x->params);
             /* struct copy */
             *x = *pt;  
             x->next = n;
             free_pdf (pt);
             free (pt);
             return 0;
          }
        p = x;
     }

   if (p == NULL)
     User_Pdf_Methods = pt;
   else p->next = pt;

   return 0;
}

/*}}}*/

void free_user_pdf_methods (void) /*{{{*/
{
   Particle_Type *q = User_Pdf_Methods;
   while (q)
     {
        Particle_Type *n = q->next;
        free_pdf (q);
        free(q);
        q = n;
     }
}

/*}}}*/

