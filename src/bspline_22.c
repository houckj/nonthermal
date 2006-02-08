#include "config.h"
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <math.h>
#include <float.h>

#include <slang.h>

#include "cfortran/cfortran.h"
#include "interp.h"

#define FEWER_CFORTRAN_COMPILER_WARNINGS \
     (void)c2fstrv; (void)f2cstrv; (void)kill_trailing; (void)vkill_trailing; \
     (void)num_elem

PROTOCCALLSFSUB4(DBSNAK,dbsnak,INT,DOUBLEV,INT,DOUBLEV)
#define DBSNAK(ndata,data,kord,xknot) \
     CCALLSFSUB4(DBSNAK,dbsnak,INT,DOUBLEV,INT,DOUBLEV,\
               ndata,data,kord,xknot)
  
PROTOCCALLSFFUN9(DOUBLE,DBS2VL,dbs2vl,DOUBLE,DOUBLE,INT,INT,DOUBLEV,DOUBLEV,INT,INT,DOUBLEVV)
#define DBS2VL(x,y,kx,ky,xknot,yknot,nx,ny,bcoef) \
            CCALLSFFUN9(DBS2VL,dbs2vl,DOUBLE,DOUBLE,INT,INT,DOUBLEV,DOUBLEV,INT,INT,DOUBLEVV,\
               x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

PROTOCCALLSFSUB11(DBS2IN,dbs2in,INT,DOUBLEV,INT,DOUBLEV,DOUBLEVV,INT,INT,INT,DOUBLEV,DOUBLEV,DOUBLEVV)
#define DBS2IN(nxdata,xdata,nydata,ydata,fdata,ldf,kxord,kyord,xknot,yknot,bscoef) \
     CCALLSFSUB11(DBS2IN,dbs2in,INT,DOUBLEV,INT,DOUBLEV,DOUBLEVV,INT,INT,INT,DOUBLEV,DOUBLEV,DOUBLEVV,\
               nxdata,xdata,nydata,ydata,fdata,ldf,kxord,kyord,xknot,yknot,bscoef)

struct Bspline_Info_Type
{
   double *xknot;   /* size nx+kx */
   double *yknot;   /* size ny+ky */
   double **bcoef;  /* size nx*ny */
   int kx, ky;      /* order of B-spline in each direction */
   int nx, ny;      /* number of B-spline coefficients in each direction */
};

double bspline_eval (Bspline_Info_Type *q, double x, double y)
{
   if ((x < q->xknot[0]) || (q->xknot[q->nx-1] < x)
       || (y < q->yknot[0]) || (q->yknot[q->ny-1] < y))
     {
        return DBL_MAX;
     }   
   FEWER_CFORTRAN_COMPILER_WARNINGS;      
   return DBS2VL(x,y,q->kx,q->ky,q->xknot,q->yknot,q->nx,q->ny,q->bcoef);
}

typedef int Init_Method_Type (Bspline_Info_Type *, int, int, double *, int, double *, int, double **);
static Init_Method_Type bspline_init;

static int bspline_init (Bspline_Info_Type *q, int kx, int ky,
                         double *x, int nx, double *y, int ny, double **f)
{
   int ldf = nx;
   q->kx = kx;
   q->ky = ky;
   DBSNAK(nx,x,kx,q->xknot);
   DBSNAK(ny,y,ky,q->yknot);
   DBS2IN(nx,x,ny,y,f,ldf,kx,ky,q->xknot,q->yknot,q->bcoef);
   return 0;
}

typedef struct
{
   SLang_Array_Type *xknot;
   SLang_Array_Type *yknot;
   SLang_Array_Type *bcoef;
   int kx, ky;
}
Bspline_Init_Info_Type;

static SLang_CStruct_Field_Type Bspline_Init_Info_Layout[] =
{
   MAKE_CSTRUCT_FIELD (Bspline_Init_Info_Type, xknot, "xknot", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Bspline_Init_Info_Type, yknot, "yknot", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Bspline_Init_Info_Type, bcoef, "bcoef", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Bspline_Init_Info_Type, kx, "kx", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Bspline_Init_Info_Type, ky, "ky", SLANG_INT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int pop_bspline_init_info (Bspline_Info_Type *q,  /*{{{*/
                                  int kx, int ky, double *x, int nx, double *y, int ny, double **f)
{
   Bspline_Init_Info_Type info;
   SLang_Array_Type *xknot, *yknot, *bcoef;
   int i, bcoef_ok;

   (void) x; (void) y; (void) f;

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&info, Bspline_Init_Info_Layout))
     return -1;

   xknot = info.xknot;
   yknot = info.yknot;
   bcoef = info.bcoef;
   kx = info.kx;
   ky = info.ky;
   nx = xknot->num_elements - kx;
   ny = yknot->num_elements - ky;

   bcoef_ok = ((bcoef->num_dims == 2) && (bcoef->dims[0] == ny)
               && (bcoef->dims[1] == nx));
   if (bcoef_ok == 0)
     {
        fprintf (stderr, "*** Error:  bcoef must be Double_Type[%d,%d]\n", ny,nx);
        SLang_free_cstruct ((VOID_STAR)&info, Bspline_Init_Info_Layout);
        return -1;
     }

   q->kx = kx;
   q->ky = ky;
   q->nx = nx;
   q->ny = ny;

   memcpy ((char *)q->xknot, (char *)xknot->data, xknot->num_elements * sizeof(int));
   memcpy ((char *)q->yknot, (char *)yknot->data, yknot->num_elements * sizeof(int));

   for (i = 0; i < q->ny; i++)
     {
        double *qb = q->bcoef[i];
        int j, dims[2];

        dims[0] = i;

        for (j = 0; j < q->nx; j++)
          {
             dims[1] = j;
             if (-1 == SLang_get_array_element (bcoef, dims, &qb[j]))
               {
                  SLang_free_cstruct ((VOID_STAR)&info, Bspline_Init_Info_Layout);
                  return -1;
               }
          }
     }

   SLang_free_cstruct ((VOID_STAR)&info, Bspline_Init_Info_Layout);

   return 0;
}

/*}}}*/

static int push_bspline_init_info (Bspline_Info_Type *q) /*{{{*/
{
   Bspline_Init_Info_Type info;
   SLang_Array_Type *xknot=NULL, *yknot=NULL, *bcoef=NULL;
   int i, j, dims[2], num_xk, num_yk, status;

   dims[0] = q->ny;
   dims[1] = q->nx;
   num_xk = q->nx + q->kx;
   num_yk = q->ny + q->ky;

   if ((NULL == (bcoef = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, dims, 2)))
       || (NULL == (xknot = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &num_xk, 1)))
       || (NULL == (yknot = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &num_yk, 1))))
     {
        SLang_free_array (bcoef);
        SLang_free_array (xknot);
        SLang_free_array (yknot);
        return -1;
     }

   for (i = 0; i < q->ny; i++)
     {
        double *qb = q->bcoef[i];
        dims[0] = i;

        for (j = 0; j < q->nx; j++)
          {
             dims[1] = j;
             if (-1 == SLang_set_array_element (bcoef, dims, &qb[j]))
               {
                  SLang_free_array (bcoef);
                  SLang_free_array (xknot);
                  SLang_free_array (yknot);
                  return -1;
               }
          }
     }

   memcpy ((char *)xknot->data, (char *)q->xknot, num_xk*sizeof(double));
   memcpy ((char *)yknot->data, (char *)q->yknot, num_yk*sizeof(double));

   info.xknot = xknot;
   info.yknot = yknot;
   info.bcoef = bcoef;
   info.kx = q->kx;
   info.ky = q->ky;

   status = SLang_push_cstruct ((VOID_STAR)&info, Bspline_Init_Info_Layout);
   SLang_free_array (bcoef);
   SLang_free_array (xknot);
   SLang_free_array (yknot);

   return status;
}

/*}}}*/

void bspline_free (Bspline_Info_Type *q) /*{{{*/
{
   if (q == NULL)
     return;
   free (q->xknot);
   if (q->bcoef != NULL)
     {
        free (q->bcoef[0]);
        free (q->bcoef);
     }
   free (q);
}

/*}}}*/

static Bspline_Info_Type *alloc_bspline_info_type
  (unsigned kx, unsigned ky, unsigned int nx, unsigned ny)
{
   Bspline_Info_Type *q;
   unsigned int i, n;

   if (NULL == (q = malloc (sizeof *q)))
     return NULL;
   memset ((char *)q, 0, sizeof *q);

   q->nx = nx;
   q->ny = ny;
   q->kx = kx;
   q->ky = ky;

   n = (nx + kx) + (ny + ky);

   if (NULL == (q->xknot = malloc (n*sizeof(double))))
     {
        bspline_free (q);
        return NULL;
     }
   q->yknot = q->xknot + (nx + kx);

   if ((NULL == (q->bcoef = malloc (ny*sizeof(double *))))
       || (NULL == (q->bcoef[0] = malloc ((ny*nx)*sizeof(double)))))
     {
        bspline_free (q);
        return NULL;
     }
   memset ((char *)q->bcoef[0], 0, (ny*nx)*sizeof(double));

   for (i = 1; i < ny; i++)
     {
        double *p = q->bcoef[i-1];
        q->bcoef[i] = p + q->nx;
     }

   return q;
}

static Bspline_Info_Type *dup_bspline_info (Bspline_Info_Type *q)
{
   Bspline_Info_Type *qc;
   int i, n;

   if (q == NULL)
     return NULL;

   if (NULL == (qc = alloc_bspline_info_type (q->kx, q->ky, q->nx, q->ny)))
     return NULL;

   n = (q->nx + q->kx) + (q->ny + q->ky);
   memcpy ((char *)qc->xknot, (char *)q->xknot, n * sizeof(double));

   for (i = 0; i < q->ny; i++)
     {
        int j;
        double *q_b = q->bcoef[i];
        double *qc_b = qc->bcoef[i];
        for (j = 0; j < q->nx; j++)
          {
             qc_b[j] = q_b[j];
          }
     }

   return qc;
}

static Bspline_Info_Type *generic_bspline_open /*{{{*/
  (int kx, int ky, double *x, int nx, double *y, int ny, double **f,
   Init_Method_Type *init_method)
{
   Bspline_Info_Type *q;

   if (NULL == (q = alloc_bspline_info_type (kx, ky, nx, ny)))
     return NULL;

   if (-1 == (*init_method)(q, kx, ky, x, nx, y, ny, f))
     {
        bspline_free (q);
        return NULL;
     }

   return q;
}

/*}}}*/

Bspline_Info_Type *bspline_open /*{{{*/
  (int kx, int ky, double *x, int nx, double *y, int ny, double **f)
{
   return generic_bspline_open (kx, ky, x, nx, y, ny, f, &bspline_init);
}

int Bspline_Type_Id = -1;

static int push_bspline_info (Bspline_Info_Type *x) /*{{{*/
{
   Bspline_Type *qt;
   SLang_MMT_Type *mmt;

   if (x == NULL)
     return -1;

   if (NULL == (qt = (Bspline_Type *)SLmalloc (sizeof *qt)))
     return -1;
   memset ((char *) qt, 0, sizeof *qt);

   qt->info = x;

   mmt = SLang_create_mmt (Bspline_Type_Id, (VOID_STAR) qt);
   if (NULL == mmt)
     {
        SLfree ((char *)qt);
        bspline_free (x);
        return -1;
     }

   if (-1 == SLang_push_mmt (mmt))
     {
        SLang_free_mmt (mmt);
        return -1;
     }

   return 0;
}

/*}}}*/

Bspline_Info_Type *copy_bspline_info (void) /*{{{*/
{
   SLang_MMT_Type *mmt;
   Bspline_Type *qt;
   Bspline_Info_Type *info;

   if (Bspline_Type_Id != SLang_peek_at_stack ())
     {
        SLdo_pop ();
        return NULL;
     }

   if ((NULL == (mmt = SLang_pop_mmt (Bspline_Type_Id)))
       || (NULL == (qt = SLang_object_from_mmt (mmt))))
     {
        SLang_free_mmt (mmt);
        return NULL;
     }

   info = dup_bspline_info (qt->info);
   SLang_free_mmt (mmt);

   return info;
}

/*}}}*/

typedef struct
{
   SLang_Array_Type *x;
   SLang_Array_Type *y;
   SLang_Array_Type *f;
}
Grid_Type;

static SLang_CStruct_Field_Type Grid_Layout[] =
{
   MAKE_CSTRUCT_FIELD (Grid_Type, x, "x", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Grid_Type, y, "y", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Grid_Type, f, "f", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int pop_grid (Grid_Type *g) /*{{{*/
{
   unsigned int nxny;

   if (-1 == SLang_pop_cstruct ((VOID_STAR)g, Grid_Layout))
     return -1;

   if ((g->x == NULL) || (g->y == NULL) || (g->f == NULL))
     return -1;

   nxny = g->x->num_elements * g->y->num_elements;
   if (nxny != g->f->num_elements)
     return -1;

   return 0;
}

/*}}}*/

void bspline_open_intrin (int *kx, int *ky) /*{{{*/
{
   Grid_Type g;
   Bspline_Info_Type *q = NULL;
   Init_Method_Type *init_method;
   double *x, *y, *pf, **f=NULL;
   int i, nx, ny;
   int status = -1;

   if (-1 == pop_grid (&g))
     return;

   x = (double *)g.x->data;
   y = (double *)g.y->data;
   nx = g.x->num_elements;
   ny = g.y->num_elements;

   if ((NULL == (f = malloc (ny*sizeof(double *))))
       || (NULL == (f[0] = malloc ((ny*nx)*sizeof(double)))))
     goto return_error;
   
   for (i = 1; i < ny; i++)
     {
        pf = f[i-1];
        f[i] = pf + nx;
     }   

   for (i = 0; i < ny; i++)
     {
        int j, dims[2];
        pf = f[i];
        dims[0] = i;
        for (j = 0; j < nx; j++)
          {
             dims[1] = j;
             if (-1 == SLang_get_array_element (g.f, dims, &pf[j]))
               goto return_error;
          }
     }

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLdo_pop ();
        init_method = &bspline_init;
     }
   else init_method = &pop_bspline_init_info;

   if (NULL == (q = generic_bspline_open (*kx, *ky, x, nx, y, ny, f, init_method)))
     goto return_error;

   if (-1 == push_bspline_info (q))
     goto return_error;

   status = 0;

   return_error:
   if (status) SLang_set_error (SL_INTRINSIC_ERROR);
   SLang_free_cstruct ((VOID_STAR)&g, Grid_Layout);
   if (f != NULL) free(f[0]);
   free(f);

   return;
}

/*}}}*/

double bspline_eval_intrin (Bspline_Type *qt, double *x, double *y) /*{{{*/
{
   return bspline_eval (qt->info, *x, *y);
}

/*}}}*/

void bspline_init_info_intrin (Bspline_Type *qt) /*{{{*/
{
   if (-1 == push_bspline_init_info (qt->info))
     SLang_set_error (SL_INTRINSIC_ERROR);
}

/*}}}*/

void destroy_bspline_type (SLtype type, VOID_STAR f) /*{{{*/
{
   Bspline_Type *qt = (Bspline_Type *) f;
   (void) type;
   if (qt == NULL)
     return;
   bspline_free (qt->info);
   qt->info = NULL;
   SLfree ((char *) qt);
}

/*}}}*/
