/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Jan 2006
 */

#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <slang.h>

#include "interp.h"

#include "_nonthermal.h"
#include "inverse_compton.h"
#include "ic_table.h"

typedef struct IC_Table_Type IC_Table_Type;
struct IC_Table_Type
{
   IC_Table_Type *next;
   Bspline_Info_Type *q;
   double xscale_log10[2];
   double yscale_log10[2];
   double *x_bdry, *y_bdry;
   double sigma0, xlg_epsilon;
   double dilution_factor;
   unsigned int n;
};

int ic_push_table (Inverse_Compton_Type *ic) /*{{{*/
{
   (void) ic;
   fprintf (stderr, "*** ic_push_table: not implemented\n");
   return 0;
}

/*}}}*/

void ic_free_client_data (void *v) /*{{{*/
{
   IC_Table_Type *sa = (IC_Table_Type *)v;

   while (sa != NULL)
     {
        IC_Table_Type *next = sa->next;
        bspline_free (sa->q);
        free(sa->x_bdry);
        free(sa);
        sa = next;
     }
}

/*}}}*/

static void *pop_table (char *file) /*{{{*/
{
   IC_Table_Type *t = NULL;
   SLang_Array_Type *k=NULL, *x=NULL, *y=NULL;
   double *keys;
   int status;

   if (-1 == SLang_run_hooks ("invc_table_init_hook", 1, file))
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

   t->dilution_factor = 1.0;

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
        ic_free_client_data ((void *)t);
        t = NULL;
     }

   return t;
}

/*}}}*/

void ic_set_dilution_factors (void *v, double *df, unsigned int n) /*{{{*/
{
   IC_Table_Type *sa = (IC_Table_Type *)v;

   if (df == NULL)
     {
        fprintf (stderr, "dilution factors not set (null ptr)\n");
        return;
     }

   /* linked list is in reverse order from input list of tables */
   while ((sa != NULL) && (n > 0))
     {
        sa->dilution_factor = df[--n];
        sa = sa->next;
     }

   if (n != 0)
     fprintf (stderr, "*** number of dilution factors != number of tables \n");
}

/*}}}*/

static int copy_string_array (char ***a, SLang_Array_Type *s) /*{{{*/
{
   int i;

   *a = malloc ((s->num_elements + 1) * sizeof (char *));
   if (NULL == *a)
     return -1;

   for (i = 0; i < (int) s->num_elements; i++)
     {
        if (-1 == SLang_get_array_element (s, &i, (VOID_STAR) &(*a)[i]))
          {
             free(*a);
             *a=NULL;
             return -1;
          }
     }

   (*a)[s->num_elements] = NULL;    /* ensure NULL terminated */

   return 0;
}
/*}}}*/

static SLang_Array_Type *get_file_list (const char *file) /*{{{*/
{
   SLang_Array_Type *a = NULL;
   int size, pos;
   size = 1;
   pos = 0;

   if (-1 == SLang_run_hooks ("ic_get_table_names", 1, file))
     return NULL;

   if (-1 == SLang_pop_array_of_type (&a, SLANG_STRING_TYPE))
     {
        if (NULL == (a = SLang_create_array (SLANG_STRING_TYPE, 0, NULL, &size, 1)))
          return NULL;
        if (-1 == SLang_set_array_element (a, &pos, &file))
          {
             SLang_free_array (a);
             return NULL;
          }
     }

   return a;
}

/*}}}*/

void *ic_init_client_data (const char *file) /*{{{*/
{
   SLang_Array_Type *sl_fs;
   char **fs;
   char **f;
   IC_Table_Type *sa = NULL;

   if (NULL == (sl_fs = get_file_list (file)))
     return NULL;

   if (-1 == copy_string_array (&fs, sl_fs))
     {
        SLang_free_array (sl_fs);
        return NULL;
     }

   SLang_free_array (sl_fs);

   for (f = fs; *f != NULL; f++)
     {
        IC_Table_Type *s;

        if (NULL == (s = pop_table (*f)))
          goto error_return;

        s->next = sa;
        sa = s;
     }

   free (fs);

   return sa;

   error_return:
   ic_free_client_data (sa);
   free (fs);
   return NULL;
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

static int interp_photon_integral (IC_Table_Type *t, int complain, /*{{{*/
                                   double gamma, double efinal,
                                   double *value)
{
   double *xscale = t->xscale_log10;
   double *yscale = t->yscale_log10;
   double xlg, ylg, xb, xx, f;
   
   (void) complain;

   *value = 0.0;

   if (gamma > 0)
     xlg = log10(gamma);
   else return -1;

   if ((xlg < xscale[0]) || (xscale[1] < xlg))
     return 0;

   if (efinal > 0)
     ylg = log10 (efinal);
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

int ic_interp_photon_integral (Inverse_Compton_Type *ic, double *value) /*{{{*/
{
   IC_Table_Type *s = (IC_Table_Type *)ic->client_data;
   double gamma, efinal;
   int flag;

   efinal = ic->energy_final_photon;
   gamma = ic->electron_gamma;
   flag = ic->complain_on_extrapolate;

   *value = 0.0;
   while (s != NULL)
     {
        double v;

        if (-1 == interp_photon_integral (s, flag, gamma, efinal, &v))
          return -1;

        *value += v * s->dilution_factor;
        s = s->next;
     }

   return 0;
}

/*}}}*/

