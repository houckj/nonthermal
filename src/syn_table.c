/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Oct 2002
 *
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_spline.h>

#include <slang.h>

#include "_nonthermal.h"
#include "syn_table.h"

typedef struct
{
   gsl_interp_accel *accel;
   gsl_spline *spline;
   double *x;
   double *y;
   unsigned int n;
}
Table_Type;

static double angular_integrand (double t, void *p) /*{{{*/
{
   double x = *(double *)p;
   double t2;
   gsl_sf_result f;

   if (t == 0.0 || t == 1.0)
     return 0.0;

   (void) gsl_sf_synchrotron_1_e (x/t, &f);

   t2 = t*t;

   return f.val * t2 / sqrt (1.0 - t2);
}

/*}}}*/

int syn_angular_integral (double x, double *y) /*{{{*/
{
   double epsabs, epsrel, abserr;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   size_t limit;

   *y = 0.0;

   f.function = &angular_integrand;
   f.params = &x;
   epsabs = 0.0;
   epsrel = 1.e-10;
   limit = MAX_QAG_SUBINTERVALS;

   work = gsl_integration_workspace_alloc (limit);
   if (work == NULL)
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   (void) gsl_integration_qag (&f, 0.0, 1.0, epsabs, epsrel, limit,
                               GSL_INTEG_GAUSS31,
                               work, y, &abserr);

   gsl_set_error_handler (gsl_error_handler);

   gsl_integration_workspace_free (work);

   return 0;
}

/*}}}*/

void syn_free_table (void *p) /*{{{*/
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

void *syn_alloc_table (unsigned int n) /*{{{*/
{
   Table_Type *t;

   if (NULL == (t = malloc (sizeof *t)))
     return NULL;
   memset ((char *)t, 0, sizeof *t);

   t->n = n;

   if (NULL == (t->x = malloc (2 * n * sizeof (double))))
     {
        syn_free_table (t);
        return NULL;
     }

   t->y = t->x + n;

   t->spline = gsl_spline_alloc (gsl_interp_akima, t->n);
   t->accel = gsl_interp_accel_alloc ();
   if (t->spline == NULL || t->accel == NULL)
     {
        syn_free_table (t);
        return NULL;
     }

   return t;
}

/*}}}*/

typedef struct
{
   SLang_Array_Type *x;
   SLang_Array_Type *y;
}
Table_CStruct_Type;

static SLang_CStruct_Field_Type Table_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, x, "x", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, y, "y", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

void *syn_load_table (char *file) /*{{{*/
{
   Table_CStruct_Type tc;
   Table_Type *t = NULL;

   (void) SLang_run_hooks ("fits_read_table", 1, file);

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&tc, Table_Type_Layout))
     return NULL;

   if ((tc.x == NULL || tc.y == NULL)
       || (tc.x->num_elements != tc.y->num_elements))
     return NULL;

   if (NULL == (t = syn_alloc_table (tc.x->num_elements)))
     {
        SLang_free_cstruct ((VOID_STAR)&tc, Table_Type_Layout);
        return NULL;
     }

   memcpy ((char *)t->x, (char *)tc.x->data, t->n * sizeof(double));
   memcpy ((char *)t->y, (char *)tc.y->data, t->n * sizeof(double));
   SLang_free_cstruct ((VOID_STAR)&tc, Table_Type_Layout);

   gsl_spline_init (t->spline, t->x, t->y, t->n);

   /* fprintf (stderr, "loaded %s\n", file); */

   return t;
}

/*}}}*/

int syn_push_table (void *p) /*{{{*/
{
   Table_Type *t = (Table_Type *)p;
   Table_CStruct_Type tc;
   int num = t->n;

   tc.x = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &num, 1);
   if (NULL == tc.x)
     return -1;

   tc.y = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &num, 1);
   if (NULL == tc.y)
     {
        SLang_free_array (tc.x);
        return -1;
     }

   memcpy ((char *)tc.x->data, (char *)t->x, num*sizeof(double));
   memcpy ((char *)tc.y->data, (char *)t->y, num*sizeof(double));

   if (-1 == SLang_push_cstruct ((VOID_STAR)&tc, Table_Type_Layout))
     {
        fprintf (stderr, "*** syn_write_table:  Error pushing cstruct\n");
        SLang_free_array (tc.x);
        SLang_free_array (tc.y);
        return -1;
     }

   return 0;
}

/*}}}*/

int syn_create_table (void *p) /*{{{*/
{
   Table_Type *t = (Table_Type *)p;
   unsigned int i;
   double lg_xmin = -38.0;
   double lg_xmax =   1.5;

   for (i = 0; i < t->n; i++)
     {
        double y, lgx;
        lgx = lg_xmin + (lg_xmax - lg_xmin) * i / (t->n - 1.0);
        (void) syn_angular_integral (pow(10.0,lgx), &y);
        t->x[i] = lgx;
        t->y[i] = (y > DBL_MIN) ? log10(y) : (double) DBL_MIN_10_EXP;
        fprintf (stderr, "%4d/%4d\r", i+1, t->n);
     }
   fputc ('\n', stderr);

   return 0;
}

/*}}}*/

int syn_interp_angular_integral (void *p, double x, double *y) /*{{{*/
{
   Table_Type *t = (Table_Type *)p;
   double lgx;

   *y = 0.0;

   if ((x <= 0.0) || t == NULL)
     return 0;

   lgx = log10(x);

   if (lgx < t->x[0] || t->x[t->n-1] <= lgx)
     return 0;

   if (-1 == gsl_spline_eval_e (t->spline, lgx, t->accel, y))
     return -1;

   *y = pow (10.0, *y);

   return 0;
}

/*}}}*/

