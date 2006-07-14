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
#include <math.h>

extern 
int bisection (double (*func)(double, void *), double a, double b, void *cd, double *xp);  
int bisection (double (*func)(double, void *), double a, double b, void *cd, double *xp)
{
   unsigned int count = 1;
   unsigned int bisect_count = 5;
   double fa, fb;
   double x = a;

   if (a > b)
     {
        double tmp = a; a = b; b = tmp;
     }

   a = nextafter (a, b);
   b = nextafter (b, a);

   fa = (*func)(a, cd);
   fb = (*func)(b, cd);

   if (fa * fb > 0)
     {
        /* "bisection: the interval may not bracket a root: f(a)*f(b)>0" */
        *xp = a;
        return -1;
     }

   while (b > a)
     {
        double fx;

        if (fb == 0.0)
          {
             x = b;
             break;
          }
        if (fa == 0.0)
          {
             x = a;
             break;
          }

        if (count % bisect_count)
          {
             x = (a*fb - b*fa) / (fb - fa);
             if (x < a || b < x)
               x = 0.5 * (a + b);
          }
        else x = 0.5 * (a + b);

        x = nextafter (x, a);

        if (x <= a || b <= x)
          break;

        fx = (*func)(x, cd);
        count++;

        if (fx*fa < 0.0)
          {
             fb = fx;
             b = x;
          }
        else
          {
             fa = fx;
             a = x;
          }
     }

   *xp = x;
   return 0;
}

#ifdef TESTING

#include <float.h>
#include <stdio.h>

static double f (double x, void *cd)
{
   (void) cd;

   return x - cos(x);
}

int main (void)
{
   double x;

   (void) bisection (&f, 0.0, 1.0, NULL, &x);

   if (fabs(f(x,NULL)) > DBL_EPSILON)
     fprintf (stdout, "Error\n");

   fprintf (stdout, "x = %15.13f\n", x);

   return 0;
}

#endif
