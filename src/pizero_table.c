/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  July 2005
 *
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <slang.h>

#include "_nonthermal.h"
#include "pizero_table.h"

typedef struct
{
   gsl_interp_accel *accel;
   gsl_spline *spline;
   double *x;
   double *y;
   unsigned int n;
}
Table_Type;

void pizero_free_table (void *p) /*{{{*/
{
   Table_Type *t = (Table_Type *)p;

   if (t == NULL)
     return;

   gsl_interp_accel_free (t->accel);
   gsl_spline_free (t->spline);
   free (t->x);
   free (t);
}

/*}}}*/

void *pizero_alloc_table (unsigned int n) /*{{{*/
{
   Table_Type *t;

   if (NULL == (t = malloc (sizeof *t)))
     return NULL;
   memset ((char *)t, 0, sizeof *t);

   t->n = n;

   if (NULL == (t->x = malloc (2 * n * sizeof (double))))
     {
        pizero_free_table (t);
        return NULL;
     }

   t->y = t->x + n;

   t->spline = gsl_spline_alloc (gsl_interp_cspline, t->n);
   t->accel = gsl_interp_accel_alloc ();
   if (t->spline == NULL || t->accel == NULL)
     {
        pizero_free_table (t);
        return NULL;
     }

   return t;
}

/*}}}*/

int pizero_spline_table (void *pt, double *x, double *y, unsigned int n) /*{{{*/
{
   Table_Type *t = (Table_Type *)pt;

   if (x == NULL || y == NULL || t == NULL || t->n != n)
     return -1;

   memcpy ((char *)t->x, (char *)x, n * sizeof(double));
   memcpy ((char *)t->y, (char *)y, n * sizeof(double));

   gsl_spline_init (t->spline, t->x, t->y, t->n);

   return 0;
}

/*}}}*/

int pizero_interp_pizero_integral (void *p, double x, double *yy) /*{{{*/
{
   Table_Type *t = (Table_Type *)p;
   double lgx, y;

   *yy = 0.0;

   if ((x <= 0.0) || t == NULL)
     return 0;

   lgx = log10(x);

   if (lgx < t->x[0] || t->x[t->n-1] <= lgx)
     return 0;

   if (-1 == gsl_spline_eval_e (t->spline, lgx, t->accel, &y))
     return -1;

   *yy = pow (10.0, y);

   return 0;
}

/*}}}*/

