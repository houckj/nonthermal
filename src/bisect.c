#include "config.h"
#include <math.h>
#include <stdio.h>

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

