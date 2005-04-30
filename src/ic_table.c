/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <slang.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "_nonthermal.h"
#include "inverse_compton.h"
#include "ic_table.h"

#define MAX_LOG_GAMMA_CHANGE  (0.025)
#define MAX_LOG_EFINAL_CHANGE (0.025)
#define MAX_LOG_FUNC_CHANGE   (0.125)

#define MIN_FUNC_VALUE (1.0e-38)

/* final photon energy, units of mc^2 */
#define EFINAL_SCALE (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY)
#define EFINAL_MIN   (1.e2 * EFINAL_SCALE)
#define EFINAL_MAX   (1.e15 * EFINAL_SCALE)

#define MIN_EXPONENT   (-737.0)
#define SAFE_LOG(x) (((x) > 0.0) ? log(x) : (double) MIN_EXPONENT)

typedef struct IC_Spline_Type IC_Spline_Type;
struct IC_Spline_Type
{
   IC_Spline_Type *next;
   double dilution_factor;
   gsl_spline **spline;
   gsl_interp_accel **accel;
   gsl_interp_accel *gamma_accel;
   double *gammas;
   unsigned int num_gammas;
   double gamma_min, gamma_max;
   double efinal_min, efinal_max;
};

typedef struct
{
   double *gammas;
   double **efinals;
   double **y;
   unsigned int *num_efinals;
   unsigned int num_gammas;
   double gamma_min, gamma_max;
   double efinal_min, efinal_max;
}
Table_Type;

#define NULL_TABLE_TYPE {NULL,NULL,NULL,NULL,0,0.0,0.0,0.0,0.0}

typedef struct
{
   double *x;
   unsigned int n;
   unsigned int size;
}
Array_Type;

#define NULL_ARRAY_TYPE {NULL,0,0}

static int append_to_array (Array_Type *a, double value) /*{{{*/
{
   if (a->n == a->size)
     {
        double *tmp = realloc (a->x, 2 * a->size * sizeof(double));
        if (tmp == NULL) return -1;
        a->x = tmp;
        a->size *= 2;
     }

   a->x[a->n++] = value;

   return 0;
}

/*}}}*/

static int alloc_array (Array_Type *a) /*{{{*/
{
   a->n = 0;
   a->size = 128;
   if (NULL == (a->x = malloc (a->size * sizeof(double))))
     return -1;

   return 0;
}

/*}}}*/

static void free_array (Array_Type *a) /*{{{*/
{
   free (a->x);
}

/*}}}*/

static double choose_next_energy (double last_e, double last_f, double e, double f, double efinal_max) /*{{{*/
{
   double factor, next_e;

   if (f > 0 && last_f > 0)
     {
        double de, dfde;

        dfde = log(f/last_f) / log(e/last_e);
        de = 0.5 * MAX_LOG_EFINAL_CHANGE * log(sqrt(f*last_f)) / dfde;
        factor = exp(de);

        if ((factor <= 1.0) || (10.0 < factor))
          factor = 2.0;
     }
   else factor = 2.0;

   if ((e*factor > efinal_max) && (f > 0))
     next_e = 0.5 * (e + efinal_max);
   else
     next_e = e * factor;

   return next_e;
}

/*}}}*/

static int map_vs_efinal (Inverse_Compton_Type *ic, double gamma, /*{{{*/
                          double efinal_min, double efinal_max,
                          unsigned int *num_efinals, double **efinals,
                          double **integrals)
{
   enum {MIN_NUM_POINTS = 5};
   Array_Type x = NULL_ARRAY_TYPE;
   Array_Type y = NULL_ARRAY_TYPE;
   unsigned int n = 0;
   double f, last_f;
   double e, last_e;

   *num_efinals = 0;
   *efinals = NULL;
   *integrals = NULL;

   if (-1 == alloc_array (&x)
       || -1 == alloc_array (&y))
     goto return_error;

   e = efinal_min;

   ic->energy_final_photon = e;
   ic->gamma_electron = gamma;
   if (-1 == ic_integral_over_incident_photons (ic, &f))
     goto return_error;

   (void) append_to_array (&x, e);
   (void) append_to_array (&y, f);

   last_e = e;
   last_f = f;

   e = 1.1 * last_e;
   do
     {
        unsigned int k;
        double next_e;

        k = 0;
        do
          {
             ic->energy_final_photon = e;
             if (-1 == ic_integral_over_incident_photons (ic, &f))
               goto return_error;

             if ((last_f == 0.0)
                 && ((f > 0.0)
                     || ((f == 0.0) && (k > 3))))
               break;

             if ((f > 0)
                 && (fabs(log(f/last_f)) < MAX_LOG_FUNC_CHANGE))
               break;

             e = 0.5 * (e + last_e);
          }
        while (k++ < 50);

        if ((n < MIN_NUM_POINTS) || (f != last_f))
          {
             if ((-1 == append_to_array (&x, e))
                 || (-1 == append_to_array (&y, f)))
               goto return_error;
             n++;

             if (n > MIN_NUM_POINTS && f < MIN_FUNC_VALUE)
               break;
          }

        next_e = choose_next_energy (last_e, last_f, e, f, efinal_max);

        last_e = e;
        last_f = f;

        e = next_e;

     }
   while (e < efinal_max);

   *num_efinals = x.n;
   *efinals = x.x;
   *integrals = y.x;
   return 0;

   return_error:
   free_array (&x);
   free_array (&y);
   return -1;
}

/*}}}*/

static int compute_table (Inverse_Compton_Type *ic, Table_Type *t) /*{{{*/
{
   unsigned int g;

   for (g = 0; g < t->num_gammas; g++)
     {
        double *ef, *y;
        unsigned int n;

        if (-1 == map_vs_efinal (ic, t->gammas[g], t->efinal_min, t->efinal_max,
                                 &n, &ef, &y))
          return -1;

        t->num_efinals[g] = n;
        t->efinals[g] = ef;
        t->y[g] = y;

        fprintf (stderr, "%3d/%3d [%3d]\r", g, t->num_gammas-1, n);
     }

   fputc ('\n',stderr);

   return 0;
}

/*}}}*/

static int alloc_table (Table_Type *t, unsigned int num) /*{{{*/
{
   t->num_gammas = num;

   if ((NULL == (t->gammas = malloc (num * sizeof (double))))
       || (NULL == (t->num_efinals = malloc (num * sizeof (unsigned int))))
       || (NULL == (t->efinals = malloc (num * sizeof (double *))))
       || (NULL == (t->y = malloc (num * sizeof (double *)))))
     return -1;

   return 0;
}

/*}}}*/

static void free_table (Table_Type *t) /*{{{*/
{
   unsigned int i;

   if (t == NULL)
     return;

   for (i = 0; i < t->num_gammas; i++)
     {
        free (t->efinals[i]);
        free (t->y[i]);
     }

   free (t->gammas);
   free (t->num_efinals);
   free (t->efinals);
   free (t->y);
}

/*}}}*/

static void compute_logs (double *x, unsigned int n) /*{{{*/
{
   unsigned int i;

   for (i = 0; i < n; i++)
     {
        double v = x[i];
        x[i] = SAFE_LOG(v);
     }
}

/*}}}*/

typedef struct
{
   SLang_Array_Type *gamma_range;
   SLang_Array_Type *efinal_range;
   SLang_Array_Type *gammas;
   SLang_Array_Type *efinals;
   SLang_Array_Type *y;
}
Table_CStruct_Type;

static SLang_CStruct_Field_Type Table_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, gamma_range, "gamma_range", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, efinal_range, "efinal_range", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, gammas, "gammas", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, efinals, "efinals", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, y, "y", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int read_table (const char *file, Table_Type *t)  /*{{{*/
{
   Table_CStruct_Type tc;
   double *gamma_range, *efinal_range;
   int i, num;
   int status = -1;

   (void) SLang_run_hooks ("ic_read_table_hook", 1, file);

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&tc, Table_Type_Layout))
     return -1;

   gamma_range = (double *)tc.gamma_range->data;
   efinal_range = (double *)tc.efinal_range->data;

   t->gamma_min = gamma_range[0];
   t->gamma_max = gamma_range[1];
   t->efinal_min = efinal_range[0];
   t->efinal_max = efinal_range[1];

   t->num_gammas = tc.gammas->num_elements;
   num = t->num_gammas;

   if (-1 == alloc_table (t, num))
     goto return_error;

   memcpy ((char *)t->gammas, (char *)tc.gammas->data, num * sizeof(double));
   compute_logs (t->gammas, num);

   for (i = 0; i < num; i++)
     {
        SLang_Array_Type *pe, *py;
        unsigned int ne, size;

        if (-1 == SLang_get_array_element (tc.efinals, &i, (VOID_STAR)&pe))
          goto return_error;
        if (-1 == SLang_get_array_element (tc.y, &i, (VOID_STAR)&py))
          goto return_error;

        ne = pe->num_elements;
        size = ne * sizeof(double);

        t->num_efinals[i] = ne;

        if ((NULL == (t->efinals[i] = malloc (size)))
            || (NULL == (t->y[i] = malloc (size))))
          goto return_error;

        memcpy ((char *)t->efinals[i], (char *)pe->data, size);
        memcpy ((char *)t->y[i], (char *)py->data, size);

        compute_logs (t->efinals[i], ne);
        compute_logs (t->y[i], ne);

        SLang_free_array (pe);
        SLang_free_array (py);
     }

   status = 0;
   return_error:

   if (status)
     free_table (t);

   SLang_free_cstruct ((VOID_STAR)&tc, Table_Type_Layout);

   return 0;
}

/*}}}*/

static int push_table (Table_Type *t)  /*{{{*/
{
   Table_CStruct_Type tc;
   double *gamma_range, *efinal_range;
   int i, num, status=-1, two=2;

   memset ((char *)&tc, 0, sizeof (Table_CStruct_Type));

   tc.gamma_range = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &two, 1);
   if (tc.gamma_range == NULL)
     goto push_result;

   tc.efinal_range = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &two, 1);
   if (tc.efinal_range == NULL)
     goto push_result;

   gamma_range = (double *)tc.gamma_range->data;
   efinal_range = (double *)tc.efinal_range->data;

   gamma_range[0] = t->gamma_min;
   gamma_range[1] = t->gamma_max;
   efinal_range[0] = t->efinal_min;
   efinal_range[1] = t->efinal_max;

   num = t->num_gammas;

   tc.gammas = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &num, 1);
   if (tc.gammas == NULL)
     goto push_result;
   memcpy ((char *)tc.gammas->data, (char *)t->gammas, num*sizeof(double));

   tc.efinals = SLang_create_array (SLANG_ARRAY_TYPE, 1, NULL, &num, 1);
   if (tc.efinals == NULL)
     goto push_result;
   tc.y = SLang_create_array (SLANG_ARRAY_TYPE, 1, NULL, &num, 1);
   if (tc.y == NULL)
     goto push_result;

   for (i = 0; i < num; i++)
     {
        SLang_Array_Type *pe;
        SLang_Array_Type *py;
        int ne = t->num_efinals[i];

        pe = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &ne, 1);
        if (pe == NULL)
          goto push_result;

        py = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &ne, 1);
        if (py == NULL)
          goto push_result;

        memcpy ((char *)pe->data, (char *)t->efinals[i], ne*sizeof(double));
        memcpy ((char *)py->data, (char *)t->y[i], ne*sizeof(double));

        ((SLang_Array_Type **)tc.efinals->data)[i] = pe;
        ((SLang_Array_Type **)tc.y->data)[i] = py;
     }

   status = 0;
   push_result:

   if (status)
     {
        SLang_free_cstruct ((VOID_STAR)&tc, Table_Type_Layout);
        memset ((char *)&tc, 0, sizeof (Table_CStruct_Type));
     }

   SLang_push_cstruct ((VOID_STAR)&tc, Table_Type_Layout);

   return status;
}

/*}}}*/

static void make_log_grid (double *x, double xmin, double xmax, unsigned int n) /*{{{*/
{
   double dlgx, min_lgx;
   unsigned int i;

   if (n < 2)
     return;

   min_lgx = log(xmin);
   dlgx = log(xmax/xmin) / (n - 1.0);

   for (i = 0; i < n; i++)
     {
        double lgx = min_lgx + i * dlgx;
        x[i] = exp(lgx);
     }
}

/*}}}*/

int ic_push_table (Inverse_Compton_Type *ic) /*{{{*/
{
   Table_Type t;
   unsigned int num;
   int status = -1;

   t.gamma_min = GAMMA_MIN_DEFAULT;
   t.gamma_max = GAMMA_MAX_DEFAULT;
   t.efinal_min = EFINAL_MIN;
   t.efinal_max = EFINAL_MAX;

   num = log(t.gamma_max/t.gamma_min) / MAX_LOG_GAMMA_CHANGE;

   if (-1 == alloc_table (&t, num))
     return -1;

   make_log_grid (t.gammas, t.gamma_min, t.gamma_max, num);

   if (-1 == compute_table (ic, &t))
     goto return_error;

   if (-1 == push_table (&t))
     goto return_error;

   status = 0;
   return_error:

   free_table (&t);

   return status;
}

/*}}}*/

static void free_splines (IC_Spline_Type *sa) /*{{{*/
{
   unsigned int i;

   if (sa == NULL)
     return;

   gsl_interp_accel_free (sa->gamma_accel);

   for (i = 0; i < sa->num_gammas; i++)
     {
        gsl_spline_free (sa->spline[i]);
        gsl_interp_accel_free (sa->accel[i]);
     }

   free (sa->spline);
   free (sa->accel);
   free (sa->gammas);
   free (sa);
}

/*}}}*/

void ic_free_client_data (void *v) /*{{{*/
{
   IC_Spline_Type *sa = (IC_Spline_Type *)v;

   while (sa != NULL)
     {
        IC_Spline_Type *next = sa->next;
        free_splines (sa);
        sa = next;
     }
}

/*}}}*/

void ic_set_dilution_factors (void *v, double *df, unsigned int n) /*{{{*/
{
   IC_Spline_Type *sa = (IC_Spline_Type *)v;

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

static IC_Spline_Type *alloc_splines (unsigned int num) /*{{{*/
{
   IC_Spline_Type *sa;
   int status = -1;

   if (NULL == (sa = malloc (sizeof *sa)))
     return NULL;
   sa->num_gammas = num;
   sa->next = NULL;
   sa->dilution_factor = 1.0;

   if (NULL == (sa->spline = malloc (num * sizeof *sa->spline)))
     goto return_error;
   memset ((char *)sa->spline, 0, num * sizeof *sa->spline);

   if (NULL == (sa->accel = malloc (num * sizeof *sa->accel)))
     goto return_error;
   memset ((char *)sa->accel, 0, num * sizeof *sa->accel);

   if (NULL == (sa->gammas = malloc (num * sizeof *sa->gammas)))
     goto return_error;
   memset ((char *)sa->gammas, 0, num * sizeof *sa->gammas);

   if (NULL == (sa->gamma_accel = gsl_interp_accel_alloc ()))
     goto return_error;

   status = 0;
   return_error:

   if (status)
     {
        ic_free_client_data (sa);
        sa = NULL;
     }

   return sa;
}

/*}}}*/

static IC_Spline_Type *compute_splines (Table_Type *t) /*{{{*/
{
   IC_Spline_Type *sa;
   unsigned int i;
   int status = -1;

   if (NULL == (sa = alloc_splines (t->num_gammas)))
     return NULL;

   memcpy ((char *)sa->gammas, (char *)t->gammas, t->num_gammas * sizeof(double));

   sa->gamma_min = t->gamma_min;
   sa->gamma_max = t->gamma_max;
   sa->efinal_min = t->efinal_min;
   sa->efinal_max = t->efinal_max;

   for (i = 0; i < t->num_gammas; i++)
     {
        gsl_spline *s;
        gsl_interp_accel *a;
        unsigned int num;
        double *y, *x;

        num = t->num_efinals[i];
        x = t->efinals[i];
        y = t->y[i];

        if ((NULL == (s = gsl_spline_alloc (gsl_interp_akima, num)))
             || (0 != gsl_spline_init (s, x, y, num)))
          goto return_error;

        if (NULL == (a = gsl_interp_accel_alloc ()))
          goto return_error;

        sa->spline[i] = s;
        sa->accel[i] = a;
     }

   status = 0;
   return_error:

   if (status)
     {
        ic_free_client_data (sa);
        sa = NULL;
     }

   return sa;
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
   IC_Spline_Type *sa = NULL;
   Table_Type t = NULL_TABLE_TYPE;

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
        IC_Spline_Type *s;

        if (-1 == read_table (*f, &t))
          goto error_return;

        /* fprintf (stderr, "loaded %s\n", *f); */

        if (NULL == (s = compute_splines (&t)))
          goto error_return;

        s->next = sa;
        sa = s;

        free_table (&t);
        memset ((char *)&t, 0, sizeof (t));
     }

   free (fs);

   return sa;

   error_return:
   ic_free_client_data (sa);
   free_table (&t);
   free (fs);
   return NULL;
}

/*}}}*/

static int interp_photon_integral (IC_Spline_Type *x, int complain, /*{{{*/
                               double gamma, double efinal,
                               double *value)
{
   enum {NUM_PTS=2};
   gsl_interp_accel *a[NUM_PTS];
   gsl_spline *s[NUM_PTS];
   double y[NUM_PTS], v, f;
   size_t g[NUM_PTS], g0;
   int i;

   *value = 0.0;

   if (x == NULL)
     return -1;

   if ((gamma < x->gamma_min || x->gamma_max <= gamma)
       || (efinal < x->efinal_min || x->efinal_max <= efinal))
     {
        if (complain == 0)
          return 0;
        fprintf (stderr, "*** extrapolation not supported:  gamma=%g [%g,%g), efinal=%g [%g,%g)\n",
                 gamma, x->gamma_min, x->gamma_max,
                 efinal, x->efinal_min, x->efinal_max);
        return -1;
     }

   /* interpolate on logs */
   gamma = log (gamma);
   efinal = log (efinal);

   g0 = gsl_interp_accel_find (x->gamma_accel, x->gammas, x->num_gammas, gamma);

   g[0] = g0;
   g[1] = g0+1;

   for (i = 0; i < 2; i++)
     {
        int status;
        a[i] = x->accel[g[i]];
        s[i] = x->spline[g[i]];

        status = gsl_spline_eval_e (s[i], efinal, a[i], &y[i]);
        if (status == GSL_EDOM)
          y[i] = MIN_EXPONENT;
        else if (status)
          {
             fprintf (stderr, "*** %s\n", gsl_strerror(status));
             return -1;
          }
     }

   f = (x->gammas[g[1]] - gamma) / (x->gammas[g[1]] - x->gammas[g[0]]);
   v = y[0] * f + y[1] * (1.0 - f);

   if (v > MIN_EXPONENT)
     v = exp(v);
   else v = 0.0;

   *value = v;

   return 0;
}

/*}}}*/

int ic_interp_photon_integral (Inverse_Compton_Type *ic, double *value) /*{{{*/
{
   IC_Spline_Type *s = (IC_Spline_Type *)ic->client_data;
   double gamma, efinal;
   int flag;

   efinal = ic->energy_final_photon;
   gamma = ic->gamma_electron;
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

