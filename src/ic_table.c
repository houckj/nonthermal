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
   double sigma0, omega0, xlg_epsilon;
   double dilution_factor;
};

void ic_free_client_data (void *v) /*{{{*/
{
   IC_Table_Type *sa = (IC_Table_Type *)v;

   while (sa != NULL)
     {
        IC_Table_Type *next = sa->next;
        bspline_free (sa->q);
        free(sa);
        sa = next;
     }
}

/*}}}*/

static void *pop_table (char *file) /*{{{*/
{
   IC_Table_Type *t = NULL;
   SLang_Array_Type *k=NULL;
   double *keys;
   int status = -1;

   if (-1 == SLang_run_hooks ("invc_table_init_hook", 1, file))
     return NULL;

   if (SLANG_NULL_TYPE == SLang_peek_at_stack())
     {
        SLdo_pop();
        return NULL;
     }

   if (-1 == SLang_pop_array_of_type (&k, SLANG_DOUBLE_TYPE))
     goto return_error;

   if ((k == NULL) || (k->num_elements != 7))
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
   t->omega0          = keys[5];
   t->xlg_epsilon     = keys[6];

   if (NULL == (t->q = copy_bspline_info ()))
     goto return_error;

   status = 0;
   return_error:

   SLang_free_array (k);
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

static int interp_photon_integral (IC_Table_Type *t, int complain, /*{{{*/
                                   double gamma, double efinal,
                                   double *value)
{
   double *xscale = t->xscale_log10;
   double *yscale = t->yscale_log10;
   double xlg, ylg, xb, xx, f, gamma_min;

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

   gamma_min = 0.5 * (efinal + sqrt (efinal * (efinal + 1.0 / t->omega0)));
   xb = log10(gamma_min);

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

