#ifndef U_NONTHERMAL_H
#define U_NONTHERMAL_H

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>

#include <math.h>

#include "slang.h"

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

#include "nonthermal.h"

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

