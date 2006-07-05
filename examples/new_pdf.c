#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nonthermal.h"

#define SPEED_OF_LIGHT 3.e10
#define GEV (1.e9 * 1.602e-12)

#define GAMMA_MIN 10
#define GAMMA_MAX 2.e9

static double min_momentum (Particle_Type *pt)
{
   return GAMMA_MIN * pt->mass * SPEED_OF_LIGHT;
}

static double max_momentum (Particle_Type *pt)
{
   return GAMMA_MAX * pt->mass * SPEED_OF_LIGHT;
}

static int pdf_simple (Particle_Type *pt, double pc, double *ne)
{
   double mc2, r, e, x, g, f;

   if (pt == NULL || ne == NULL)
     return -1;

   /* params:  [index] */
   if (pt->num_params != 1)
     return -1;

   *ne = 0.0;

   mc2 = pt->mass * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
   r = pc/mc2;
   e = mc2 * r * sqrt (1.0 + 1.0 /r /r);
   x = e / GEV;
   g = pt->params[0];

   f = pow (x, -g);

   *ne = f;

   return 0;
}

NONTHERMAL_PDF_MODULE(simple,pt,options)
{
   if (pt == NULL)
     return -1;

   (void) options;

   pt->method = "simple";
   pt->num_params = 1;
   pt->spectrum = &pdf_simple;
   pt->momentum_min = &min_momentum;
   pt->momentum_max = &max_momentum;

   return 0;
}
