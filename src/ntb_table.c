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
#include "ntbrems.h"
#include "ntb_table.h"

#define MAX_LOG_EPHOTON_CHANGE   (0.05)
#define MAX_LOG_EKINETIC_CHANGE  (0.05)
#define MAX_LOG_FUNC_CHANGE      (0.05)

#define MIN_FUNC_VALUE (1.0e-50)

#define EPHOTON_SCALE (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY)
/* photon energy, units of mc^2 */
#define EPHOTON_MIN   (1.e2 * EPHOTON_SCALE)
/* min e- kinetic energy, units of mc^2 */
#define EKINETIC_MIN  (EPHOTON_MIN)
#define EPHOTON_MAX   ((GAMMA_MAX_DEFAULT-1.0)-EKINETIC_MIN)
#define EKINETIC_MAX  ((GAMMA_MAX_DEFAULT-1.0)-EPHOTON_MIN)

#define MIN_EXPONENT   (-737.0)
#define SAFE_LOG(x) (((x) > 0.0) ? log(x) : (double) MIN_EXPONENT)

typedef struct NTB_Spline_Type NTB_Spline_Type;
struct NTB_Spline_Type
{
   NTB_Spline_Type *next;
   double process_weight;
   gsl_spline **spline;
   gsl_interp_accel **accel;
   gsl_interp_accel *ephoton_accel;
   double *ephotons;
   unsigned int num_ephotons;
   double ephoton_min, ephoton_max;
   double ekinetic_min, ekinetic_max;
};

typedef struct
{
   double *ephotons;
   double **ekinetics;
   double **y;
   unsigned int *num_ekinetics;
   unsigned int num_ephotons;
   double ephoton_min, ephoton_max;
   double ekinetic_min, ekinetic_max;
   int process;
}
Table_Type;

#define NULL_TABLE_TYPE {NULL,NULL,NULL,NULL,0,0.0,0.0,0.0,0.0,0}

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

static int Num_Short_Steps;
static double choose_next_energy (double last_e, double last_f, double e, double f, double ekinetic_max) /*{{{*/
{
   double next_e, factor = MAX_LOG_EKINETIC_CHANGE;

   if (f > 0 && last_f > 0)
     {
        double dfde = log(f/last_f) / log(e/last_e);
        double r = MAX_LOG_FUNC_CHANGE / fabs(dfde);
        if (r < factor)
          factor = r;
     }

   factor = exp(factor);

   if ((e*factor > ekinetic_max) && (f > 0))
     {
        Num_Short_Steps++;
        if (Num_Short_Steps < 10)
          next_e = 0.5 * (e + ekinetic_max);
        else next_e = ekinetic_max;
     }
   else next_e = e * factor;

   return next_e;
}

/*}}}*/

static int map_vs_ekinetic (double (*f)(double, double), double ephoton, /*{{{*/
                           double ekinetic_min, double ekinetic_max,
                           unsigned int *num_ekinetics, double **ekinetics,
                           double **sigmas)
{
   enum {MIN_NUM_POINTS = 5};
   Array_Type x = NULL_ARRAY_TYPE;
   Array_Type y = NULL_ARRAY_TYPE;
   unsigned int n = 0;
   double s, last_s;
   double e, last_e;

   *num_ekinetics = 0;
   *ekinetics = NULL;
   *sigmas = NULL;

   Num_Short_Steps = 0;
   if (-1 == alloc_array (&x)
       || -1 == alloc_array (&y))
     goto return_error;

   e = ekinetic_min;
   s = (*f)(e, ephoton);

   (void) append_to_array (&x, e);
   (void) append_to_array (&y, s);

   last_e = e;
   last_s = s;

   e = 1.1 * last_e;
   do
     {
        unsigned int k;
        double next_e;

        k = 0;
        do
          {
             s = (*f)(e, ephoton);
             if ((last_s == 0.0)
                 && ((s > 0.0)
                     || ((s == 0.0) && (k > 1))))
               break;

             if ((s > 0)
                 && (fabs(log(s/last_s)) < MAX_LOG_FUNC_CHANGE))
               break;

             e = 0.5 * (e + last_e);
          }
        while (k++ < 10);

        if ((n < MIN_NUM_POINTS) || (s != last_s))
          {
             if ((-1 == append_to_array (&x, e))
                 || (-1 == append_to_array (&y, s)))
               goto return_error;
             n++;

             if ((n > MIN_NUM_POINTS) && (s < MIN_FUNC_VALUE))
               break;
          }

        next_e = choose_next_energy (last_e, last_s, e, s, ekinetic_max);

        last_e = e;
        last_s = s;

        e = next_e;
     }
   while (e < ekinetic_max);

   *num_ekinetics = x.n;
   *ekinetics = x.x;
   *sigmas = y.x;
   return 0;

   return_error:
   free_array (&x);
   free_array (&y);
   return -1;
}

/*}}}*/

static int compute_table (Table_Type *t) /*{{{*/
{
   double (*cross_section)(double, double);
   unsigned int g;

   switch (t->process)
     {
      case NTB_ee:
        cross_section = &_ntb_ee_sigma_haug;
        break;
      case NTB_ep:
        cross_section = &_ntb_ei_sigma;
        break;
      default:
        fprintf (stderr, "*** Invalid process type %d\n", t->process);
        return -1;
        break;
     }

   for (g = 0; g < t->num_ephotons; g++)
     {
        double *ef, *y;
        unsigned int n;

        if (-1 == map_vs_ekinetic (cross_section, t->ephotons[g],
                                  t->ekinetic_min, t->ekinetic_max, &n, &ef, &y))
          return -1;

        t->num_ekinetics[g] = n;
        t->ekinetics[g] = ef;
        t->y[g] = y;

        fprintf (stderr, "%3d/%3d [%3d]\r", g, t->num_ephotons-1, n);
     }

   fputc ('\n',stderr);

   return 0;
}

/*}}}*/

static int alloc_table (Table_Type *t, unsigned int num) /*{{{*/
{
   t->num_ephotons = num;

   if ((NULL == (t->ephotons = malloc (num * sizeof (double))))
       || (NULL == (t->num_ekinetics = malloc (num * sizeof (unsigned int))))
       || (NULL == (t->ekinetics = malloc (num * sizeof (double *))))
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

   for (i = 0; i < t->num_ephotons; i++)
     {
        free (t->ekinetics[i]);
        free (t->y[i]);
     }

   free (t->ephotons);
   free (t->num_ekinetics);
   free (t->ekinetics);
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
   SLang_Array_Type *ephoton_range;
   SLang_Array_Type *ekinetic_range;
   SLang_Array_Type *ephotons;
   SLang_Array_Type *ekinetics;
   SLang_Array_Type *y;
}
Table_CStruct_Type;

static SLang_CStruct_Field_Type Table_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, ephoton_range, "ephoton_range", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, ekinetic_range, "ekinetic_range", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, ephotons, "ephotons", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, ekinetics, "ekinetics", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Table_CStruct_Type, y, "y", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int read_table (const char *file, Table_Type *t)  /*{{{*/
{
   Table_CStruct_Type tc;
   double *ephoton_range, *ekinetic_range;
   int i, num;
   int status = -1;

   (void) SLang_run_hooks ("ntb_read_table_hook", 1, file);

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&tc, Table_Type_Layout))
     return -1;

   ephoton_range = (double *)tc.ephoton_range->data;
   ekinetic_range = (double *)tc.ekinetic_range->data;

   t->ephoton_min = ephoton_range[0];
   t->ephoton_max = ephoton_range[1];
   t->ekinetic_min = ekinetic_range[0];
   t->ekinetic_max = ekinetic_range[1];

   t->num_ephotons = tc.ephotons->num_elements;
   num = t->num_ephotons;

   if (-1 == alloc_table (t, num))
     goto return_error;

   memcpy ((char *)t->ephotons, (char *)tc.ephotons->data, num * sizeof(double));
   compute_logs (t->ephotons, num);

   for (i = 0; i < num; i++)
     {
        SLang_Array_Type *pe, *py;
        unsigned int ne, size;

        if (-1 == SLang_get_array_element (tc.ekinetics, &i, (VOID_STAR)&pe))
          goto return_error;
        if (-1 == SLang_get_array_element (tc.y, &i, (VOID_STAR)&py))
          goto return_error;

        ne = pe->num_elements;
        size = ne * sizeof(double);

        t->num_ekinetics[i] = ne;

        if ((NULL == (t->ekinetics[i] = malloc (size)))
            || (NULL == (t->y[i] = malloc (size))))
          goto return_error;

        memcpy ((char *)t->ekinetics[i], (char *)pe->data, size);
        memcpy ((char *)t->y[i], (char *)py->data, size);

        compute_logs (t->ekinetics[i], ne);
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
   double *ephoton_range, *ekinetic_range;
   int i, num, status=-1, two=2;

   memset ((char *)&tc, 0, sizeof (Table_CStruct_Type));

   tc.ephoton_range = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &two, 1);
   if (tc.ephoton_range == NULL)
     goto push_result;

   tc.ekinetic_range = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &two, 1);
   if (tc.ekinetic_range == NULL)
     goto push_result;

   ephoton_range = (double *)tc.ephoton_range->data;
   ekinetic_range = (double *)tc.ekinetic_range->data;

   ephoton_range[0] = t->ephoton_min;
   ephoton_range[1] = t->ephoton_max;
   ekinetic_range[0] = t->ekinetic_min;
   ekinetic_range[1] = t->ekinetic_max;

   num = t->num_ephotons;

   tc.ephotons = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &num, 1);
   if (tc.ephotons == NULL)
     goto push_result;
   memcpy ((char *)tc.ephotons->data, (char *)t->ephotons, num*sizeof(double));

   tc.ekinetics = SLang_create_array (SLANG_ARRAY_TYPE, 1, NULL, &num, 1);
   if (tc.ekinetics == NULL)
     goto push_result;
   tc.y = SLang_create_array (SLANG_ARRAY_TYPE, 1, NULL, &num, 1);
   if (tc.y == NULL)
     goto push_result;

   for (i = 0; i < num; i++)
     {
        SLang_Array_Type *pe;
        SLang_Array_Type *py;
        int ne = t->num_ekinetics[i];

        pe = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &ne, 1);
        if (pe == NULL)
          goto push_result;

        py = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &ne, 1);
        if (py == NULL)
          goto push_result;

        memcpy ((char *)pe->data, (char *)t->ekinetics[i], ne*sizeof(double));
        memcpy ((char *)py->data, (char *)t->y[i], ne*sizeof(double));

        ((SLang_Array_Type **)tc.ekinetics->data)[i] = pe;
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

int ntb_push_table (int process) /*{{{*/
{
   Table_Type t;
   unsigned int num;
   int status = -1;

   t.process = process;

   t.ephoton_min = EPHOTON_MIN;
   t.ephoton_max = EPHOTON_MAX;
   t.ekinetic_min = EKINETIC_MIN;
   t.ekinetic_max = EKINETIC_MAX;

   num = log(t.ephoton_max/t.ephoton_min) / MAX_LOG_EPHOTON_CHANGE;

   if (-1 == alloc_table (&t, num))
     return -1;

   make_log_grid (t.ephotons, t.ephoton_min, t.ephoton_max, num);

   if (-1 == compute_table (&t))
     goto return_error;

   if (-1 == push_table (&t))
     goto return_error;

   status = 0;
   return_error:

   free_table (&t);

   return status;
}

/*}}}*/

static void free_splines (NTB_Spline_Type *sa) /*{{{*/
{
   unsigned int i;

   if (sa == NULL)
     return;

   gsl_interp_accel_free (sa->ephoton_accel);

   for (i = 0; i < sa->num_ephotons; i++)
     {
        gsl_spline_free (sa->spline[i]);
        gsl_interp_accel_free (sa->accel[i]);
     }

   free (sa->spline);
   free (sa->accel);
   free (sa->ephotons);
   free (sa);
}

/*}}}*/

void ntb_free_client_data (void *v) /*{{{*/
{
   NTB_Spline_Type *sa = (NTB_Spline_Type *)v;

   while (sa != NULL)
     {
        NTB_Spline_Type *next = sa->next;
        free_splines (sa);
        sa = next;
     }
}

/*}}}*/

#if 0
void ntb_set_process_weights (void *v, double *pw, unsigned int n) /*{{{*/
{
   NTB_Spline_Type *sa = (NTB_Spline_Type *)v;

   if (pw == NULL)
     {
        fprintf (stderr, "process weights not set (null ptr)\n");
        return;
     }

   /* linked list is in reverse order from input list of tables */
   while ((sa != NULL) && (n > 0))
     {
        sa->process_weight = pw[--n];
        sa = sa->next;
     }

   if (n != 0)
     fprintf (stderr, "*** number of weights != number of tables \n");
}

/*}}}*/
#endif

static NTB_Spline_Type *alloc_splines (unsigned int num) /*{{{*/
{
   NTB_Spline_Type *sa;
   int status = -1;

   if (NULL == (sa = malloc (sizeof *sa)))
     return NULL;
   sa->num_ephotons = num;
   sa->next = NULL;
   sa->process_weight = 1.0;

   if (NULL == (sa->spline = malloc (num * sizeof *sa->spline)))
     goto return_error;
   memset ((char *)sa->spline, 0, num * sizeof *sa->spline);

   if (NULL == (sa->accel = malloc (num * sizeof *sa->accel)))
     goto return_error;
   memset ((char *)sa->accel, 0, num * sizeof *sa->accel);

   if (NULL == (sa->ephotons = malloc (num * sizeof *sa->ephotons)))
     goto return_error;
   memset ((char *)sa->ephotons, 0, num * sizeof *sa->ephotons);

   if (NULL == (sa->ephoton_accel = gsl_interp_accel_alloc ()))
     goto return_error;

   status = 0;
   return_error:

   if (status)
     {
        ntb_free_client_data (sa);
        sa = NULL;
     }

   return sa;
}

/*}}}*/

static NTB_Spline_Type *compute_splines (Table_Type *t) /*{{{*/
{
   NTB_Spline_Type *sa;
   unsigned int i;
   int status = -1;

   if (NULL == (sa = alloc_splines (t->num_ephotons)))
     return NULL;

   memcpy ((char *)sa->ephotons, (char *)t->ephotons, t->num_ephotons * sizeof(double));

   sa->ephoton_min = t->ephoton_min;
   sa->ephoton_max = t->ephoton_max;
   sa->ekinetic_min = t->ekinetic_min;
   sa->ekinetic_max = t->ekinetic_max;

   for (i = 0; i < t->num_ephotons; i++)
     {
        gsl_spline *s;
        gsl_interp_accel *a;
        unsigned int num;
        double *y, *x;

        num = t->num_ekinetics[i];
        x = t->ekinetics[i];
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
        ntb_free_client_data (sa);
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

   if (-1 == SLang_run_hooks ("ntb_get_table_names", 1, file))
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

void *ntb_init_client_data (const char *file) /*{{{*/
{
   SLang_Array_Type *sl_fs;
   char **fs;
   char **f;
   NTB_Spline_Type *sa = NULL;
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
        NTB_Spline_Type *s;

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
   ntb_free_client_data (sa);
   free_table (&t);
   free (fs);
   return NULL;
}

/*}}}*/

static int interp_sigma (NTB_Spline_Type *x, int complain, /*{{{*/
                         double ekinetic, double ephoton,
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

   if ((ephoton < x->ephoton_min || x->ephoton_max < ephoton)
       || (ekinetic < x->ekinetic_min || x->ekinetic_max < ekinetic))
     {
        if (complain == 0)
          return 0;
        fprintf (stderr, "*** extrapolation not supported:  ephoton=%15.9e [%15.9e,%15.9e), ekinetic=%15.9e [%15.9e,%15.9e)\n",
                 ephoton, x->ephoton_min, x->ephoton_max,
                 ekinetic, x->ekinetic_min, x->ekinetic_max);
        /* exit(1); */
        return -1;
     }

   /* interpolate on logs */
   ephoton = log (ephoton);
   ekinetic = log (ekinetic);

   g0 = gsl_interp_accel_find (x->ephoton_accel, x->ephotons, x->num_ephotons, ephoton);

   g[0] = g0;
   g[1] = g0+1;

   for (i = 0; i < 2; i++)
     {
        int status;
        a[i] = x->accel[g[i]];
        s[i] = x->spline[g[i]];

        status = gsl_spline_eval_e (s[i], ekinetic, a[i], &y[i]);
        if (status == GSL_EDOM)
          y[i] = MIN_EXPONENT;
        else if (status)
          {
             fprintf (stderr, "*** %s\n", gsl_strerror(status));
             return -1;
          }
     }

   f = (x->ephotons[g[1]] - ephoton) / (x->ephotons[g[1]] - x->ephotons[g[0]]);
   v = y[0] * f + y[1] * (1.0 - f);

   if (v > MIN_EXPONENT)
     v = exp(v);
   else v = 0.0;

   *value = v;

   return 0;
}

/*}}}*/

/* ekinetic = gamma - 1 */
int ntb_interp_sigma (Brems_Type *ntb, double ekinetic, double ephoton, double *sigma) /*{{{*/
{
   NTB_Spline_Type *s = (NTB_Spline_Type *)ntb->client_data;

   *sigma = 0.0;
   while (s != NULL)
     {
        double v; /* FIXME? returns zero rather than complain about extrapolation */
        if (-1 == interp_sigma (s, 0, ekinetic, ephoton, &v))
          return -1;
        *sigma += v * s->process_weight;
        s = s->next;
     }

   return 0;
}

/*}}}*/

