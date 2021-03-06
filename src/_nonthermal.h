/*
  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008 John C. Houck 

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
#ifndef U_NONTHERMAL_H
#define U_NONTHERMAL_H

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>

#include <math.h>

#include "slang.h"

#ifdef NAN
#define NT_NAN   NAN
#elif defined(INFINITY)
#define NT_NAN  (INFINITY/INFINITY)
#else
#define NT_NAN  (0.0/0.0)
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390  /* 2/sqrt(pi) */
#endif

#ifndef isnan
#define isnan(a)  ((a)!=(a))
#endif

#define MAX_QAG_SUBINTERVALS  (1<<18)

#define CBR_TEMPERATURE 2.725         /* Kelvin */
#define ELECTRON_CHARGE 4.803250e-10  /* esu */

#define KEV  (1.0e3 * GSL_CONST_CGSM_ELECTRON_VOLT)
#define MEV  (1.0e6 * GSL_CONST_CGSM_ELECTRON_VOLT)
#define GEV  (1.0e9 * GSL_CONST_CGSM_ELECTRON_VOLT)
#define TEV  (1.e12 * GSL_CONST_CGSM_ELECTRON_VOLT)

#define BARN        (1.e-24)
#define MILLIBARN   (1.e-3 * BARN)

#define C_SQUARED  \
    (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT)

#define ELECTRON_REST_ENERGY   \
    (GSL_CONST_CGSM_MASS_ELECTRON * C_SQUARED)
#define PROTON_REST_ENERGY   \
    (GSL_CONST_CGSM_MASS_PROTON * C_SQUARED)
#define PIZERO_REST_ENERGY  (134.9766 * MEV)

#define ELECTRON_RADIUS \
    (ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_REST_ENERGY)

/* Emissivities derived assuming gamma_electron >> 1
 * and a 1000 TeV seems like a reasonable upper limit.
 *
 * Note that, for high B fields (B > 10 mG) electrons
 * with gamma < 100 contribute significant flux below
 * 1 GHz.  If you're interested in that part of the
 * radio spectrum more care may be required.
 * Anyway, that's the reason I'm using gamma_min=10
 * instead of gamma_min=100.
 */
#define  GAMMA_MIN_DEFAULT   10.0
#define  GAMMA_MAX_DEFAULT  (1.e15 \
                        * (GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY))

enum
{
   ELECTRON,
   PROTON
};

typedef struct Nonthermal_Error_Type Nonthermal_Error_Type;
struct Nonthermal_Error_Type
{
   const char *error_msg;
   double value;
   double estimated_abserr;
   double allowed_abserr;
   double allowed_epsrel;
};

#include "nonthermal.h"

extern int nonthermal_error_hook (Nonthermal_Error_Type *e, Particle_Type *p,
                                  char *file, int line);

extern Particle_Type *load_pdf (char *path, char *name, char *options);
extern int append_pdf (Particle_Type *pt);
extern void free_user_pdf_methods (void);

extern void free_pdf (Particle_Type *p);
extern int init_pdf (Particle_Type *p, unsigned int type);
extern int init_pdf_params (Particle_Type *p, unsigned int type, char *method, SLang_Array_Type *sl_pars);

extern int bisection (double (*func)(double, void *), double a, double b, void *cd, double *xp);

/* S-Lang-2 compatibility */

#if (SLANG_VERSION<20000)
extern int SLang_get_error (void);
extern void SLang_set_error (int);
#endif

#if (SLANG_VERSION<20000)
# define SLang_pop_double(a)  SLang_pop_double ((a),NULL,NULL)
#endif

#endif

