/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <slang.h>

#include "interp.h"

#include "_nonthermal.h"
#include "ntb_table.h"

typedef struct
{
   Bspline_Info_Type *q;
   double xscale_log10[2];
   double yscale_log10[2];
   double sigma0;
   double xlg_epsilon;
   double *x_bdry;
   double *y_bdry;
   unsigned int n;
}
Ntb_Table_Type;

void ntb_free_client_data (void *v) /*{{{*/
{
   Ntb_Table_Type *t = (Ntb_Table_Type *)v;
   if (t == NULL)
     return;
   bspline_free (t->q);
   free(t->x_bdry);
}

/*}}}*/

static void *pop_table (const char *file) /*{{{*/
{
   Ntb_Table_Type *t = NULL;
   SLang_Array_Type *x=NULL, *y=NULL, *k=NULL;
   double *keys;
   int status = 1;

   if (-1 == SLang_run_hooks ("_ntbrem_table2_init", 1, file))
     return NULL;
   
   if (SLANG_NULL_TYPE == SLang_peek_at_stack())
     {
        SLdo_pop();
        return NULL;
     }   

   if ((-1 == SLang_pop_array_of_type (&k, SLANG_DOUBLE_TYPE))
       || (-1 == SLang_pop_array_of_type (&y, SLANG_DOUBLE_TYPE))
       || (-1 == SLang_pop_array_of_type (&x, SLANG_DOUBLE_TYPE)))
     goto return_error;

   if ((k == NULL || x == NULL || y == NULL)
       || (x->num_elements != y->num_elements)
       || (k->num_elements != 6))
     goto return_error;

   if (NULL == (t = malloc (sizeof *t)))
     goto return_error;
   memset ((char *)t, 0, sizeof *t);

   keys = k->data;
   t->xscale_log10[0] = keys[0];
   t->xscale_log10[1] = keys[1];
   t->yscale_log10[0] = keys[2];
   t->yscale_log10[1] = keys[3];
   t->sigma0          = keys[4];
   t->xlg_epsilon     = keys[5];

   t->n = x->num_elements;
   if (NULL == (t->x_bdry = malloc (2*t->n * sizeof(double))))
     goto return_error;
   t->y_bdry = t->x_bdry + t->n;

   memcpy ((char *)t->x_bdry, (char *)x->data, t->n * sizeof(double));
   memcpy ((char *)t->y_bdry, (char *)y->data, t->n * sizeof(double));

   if (NULL == (t->q = copy_bspline_info ()))
     goto return_error;

   status = 0;
   return_error:

   SLang_free_array (k);
   SLang_free_array (y);
   SLang_free_array (x);
   if (status)
     {
        ntb_free_client_data ((void *)t);
        t = NULL;
     }

   return t;
}

/*}}}*/

void *ntb_init_client_data (const char *file) /*{{{*/
{
   return pop_table (file);
}

/*}}}*/

int ntb_push_table (int process) /*{{{*/
{
   (void) process;
   fprintf (stderr, "ntb_push_table:  not implemented\n");
   return 0;
}

/*}}}*/

static int bsearch_d (double t, double *x, int n) /*{{{*/
{
   int n0, n1, n2;
   double xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

/*}}}*/

/* taken from jdmath by John Davis <davis@space.mit.edu> */
static double interpolate_d (double x, double *xp, double *yp, unsigned int n) /*{{{*/
{
   unsigned int n0, n1;
   double x0, x1;

   n0 = bsearch_d (x, xp, n);

   x0 = xp[n0];
   n1 = n0 + 1;

   if (x == x0)
     return yp[n0];
   if (n1 == n)
     {
        if (n0 == 0)
          return yp[n0];
        n1 = n0 - 1;
     }

   x1 = xp[n1];
   if (x1 == x0) return yp[n0];

   return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * (x - x0);
}

/*}}}*/

int ntb_interp_sigma (Brems_Type *ntb, double gm1, double ephoton, double *value) /*{{{*/
{
   Ntb_Table_Type *t = (Ntb_Table_Type *)ntb->client_data;
   double *xscale = t->xscale_log10;
   double *yscale = t->yscale_log10;
   double xlg, ylg, xb, xx, f;

   *value = 0.0;

   /* [gm1] = gamma-1 = electron kinetic energy in units of m_e c^2
    * [ephoton] = photon energy in units of m_e c^2
    * [sigma] cross-section in cm^2
    */
   if (gm1 > 0)
     xlg = log10(gm1);
   else return -1;

   if ((xlg < xscale[0]) || (xscale[1] < xlg))
     return 0;

   if (ephoton > 0)
     ylg = log10 (ephoton);
   else return -1;

   if ((ylg < yscale[0]) || (yscale[1] < ylg))
     return 0;

   xb = interpolate_d (ylg, t->y_bdry, t->x_bdry, t->n);

   if (xlg - xb > t->xlg_epsilon)
     xx = log (xlg - xb);
   else return 0;

   f = bspline_eval (t->q, xx, ylg);
   
   if (f == DBL_MAX)
     return 0;

   *value = t->sigma0 * pow (10.0, f);

   return 0;
}

/*}}}*/
