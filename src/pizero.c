/* -*- mode: C; mode: fold -*- */
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

/* Neutral pion production cross-sections from
 *   Dermer, C.D., 1986, A&A, 157, 223.
 *   Blattnig et al, NASA Technical Report, NASA/TP-2000-210640
 *   Aharonian and Atoyan, 2000, A&A, 362, 937
 *
 * See also
 *   Domingo-Santamaria & Torres, 2005, astro-ph/0506240
*/

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "_nonthermal.h"
#include "pizero.h"
#include "pizero_table.h"

int Pizero_Method = 0;

#define TWO_PROTON_REST_ENERGY   (2 * PROTON_REST_ENERGY)
#define SQR_PROTON_REST_ENERGY  (PROTON_REST_ENERGY * PROTON_REST_ENERGY)
#define SQR_PIZERO_REST_ENERGY  (PIZERO_REST_ENERGY * PIZERO_REST_ENERGY)
#define FOUR_SQR_PROTON_REST_ENERGY  (4 * SQR_PIZERO_REST_ENERGY)
#define MASS_RATIO              (PIZERO_REST_ENERGY / TWO_PROTON_REST_ENERGY)
#define PIZERO_MASS_FACTOR      (FOUR_SQR_PROTON_REST_ENERGY * (1.0 + MASS_RATIO) * (1.0 - MASS_RATIO))
#define SQR_PIZERO_MASS_FACTOR  (PIZERO_MASS_FACTOR * PIZERO_MASS_FACTOR)

/* Placing MIN_PHOTON_ENERGY just above the threshold
 * eliminates the sharp peak in the photon spectrum near 70 MeV(?)
 */
#define MIN_PHOTON_ENERGY (0.5*PIZERO_REST_ENERGY *2)

/* \Delta_{3/2}(1232) resonance width [GeV]
 * Particle Data Group (1984)
 */
#define RESONANCE_WIDTH         ((0.115 * GEV)*0.5)
#define RESONANCE_REST_ENERGY   (1233.0 * MEV)

/* a^2 - b^2 */
static double diff_sqr (double a, double b)
{
   double r, f = 0.0;

   if (a > b)
     {
        r = b/a;
        f = a * a * (1.0 + r) * (1.0 - r);
     }
   else if (b > a)
     {
        r = a/b;
        f = b * b * (r + 1.0) * (r - 1.0);
     }

   return f;
}

/* sqrt (a^2 - b^2) */
static double root_diff_sqr (double a, double b)
{
   return sqrt(diff_sqr (a,b));
}

typedef struct
{
   double T_pi;
   double T_p;
   double s;
}
Collision_Info_Type;

static void init_collision_info (Collision_Info_Type *info, double T_p, double T_pi) /*{{{*/
{
   info->T_p = T_p;
   info->T_pi = T_pi;
   info->s = TWO_PROTON_REST_ENERGY * (T_p + TWO_PROTON_REST_ENERGY);
}

/*}}}*/

static double pizero_dermer_total_xsec (double T_p) /*{{{*/
{
   double ep = T_p + PROTON_REST_ENERGY;
   double pc = root_diff_sqr(ep, PROTON_REST_ENERGY);
   double s, sigma = 0.0;

   pc /= GEV;

   /* Dermer, C.D., 1986, A&A, 157, 223. */
   if (pc < 0.78)
     sigma = 0.0;
   else if (0.78 <= pc && pc < 0.96)
     {
        double x, eta2, eta4, tem;
        tem = 2*PIZERO_REST_ENERGY;
        s = (T_p + TWO_PROTON_REST_ENERGY) / TWO_PROTON_REST_ENERGY;
        x = T_p - SQR_PIZERO_REST_ENERGY / TWO_PROTON_REST_ENERGY;
        eta2 = diff_sqr(x,tem)/(4*SQR_PIZERO_REST_ENERGY*s);
        eta4 = eta2 * eta2;
        sigma = eta2 * (0.032 + eta4 * (0.040 + eta2 * 0.047));
     }
   else if (0.96 <= pc && pc < 1.27)
     sigma = 32.6 * pow(pc - 0.8, 3.21);
   else if (1.27 <= pc && pc < 8.0)
     sigma = 5.40 * pow(pc - 0.8, 0.81);
   else if (8.0 <= pc && pc < 1.e3)
     sigma = 32.0 * log(pc) + 48.5/sqrt(pc) - 59.5;
   else if (1.e3 <= pc)
     {
        /* Mori (1997), ApJ, 478, 225 */
        s = TWO_PROTON_REST_ENERGY * (T_p + TWO_PROTON_REST_ENERGY);
        /* s /= GEV; */
        sigma = 163.0 * pow(s/1876.0, 0.21);
     }

   return sigma * MILLIBARN;
}

/*}}}*/

static double bta (double g) /*{{{*/
{
   if (g > 1.0)
     return root_diff_sqr(g ,1.0)/g;
   else return 0.0;
}

/*}}}*/

static double gamma_fun (double g1, double g2, int sgn) /*{{{*/
{
   return g1 * g2 * (1.0 + sgn * bta(g1) * bta(g2));
}

/*}}}*/

#define COMMON_FACTORS \
   Collision_Info_Type *info = (Collision_Info_Type *)x; \
   double T_pi = info->T_pi; \
   double s = info->s; \
   double root_s = sqrt(s); \
   double g_pi = 1.0 + T_pi / PIZERO_REST_ENERGY; \
   double g_c = root_s / TWO_PROTON_REST_ENERGY; \
   double g_d_cm = (s + diff_sqr(m_d,PIZERO_REST_ENERGY)) / (2 * root_s * m_d); \
   double g_pi_i = ((m_d * m_d + diff_sqr(PIZERO_REST_ENERGY,PROTON_REST_ENERGY)) \
                       / (2 * PIZERO_REST_ENERGY * m_d))

#define PLUS_COMMON \
   COMMON_FACTORS; \
   double g_d_p = gamma_fun (g_c, g_d_cm, +1);

#define MINUS_COMMON \
   COMMON_FACTORS; \
   double g_d_m = gamma_fun (g_c, g_d_cm, -1);

static double h_plus_upper (double m_d, void *x) /*{{{*/
{
   PLUS_COMMON;
   return g_pi - gamma_fun (g_d_p, g_pi_i, -1);
}

/*}}}*/

static double h_plus_lower (double m_d, void *x) /*{{{*/
{
   PLUS_COMMON;
   return g_pi - gamma_fun (g_d_p, g_pi_i, +1);
}

/*}}}*/

static double h_minus_upper (double m_d, void *x) /*{{{*/
{
   MINUS_COMMON;
   return g_pi - gamma_fun (g_d_m, g_pi_i, -1);
}

/*}}}*/

static double h_minus_lower (double m_d, void *x) /*{{{*/
{
   MINUS_COMMON;
   return g_pi - gamma_fun (g_d_m, g_pi_i, +1);
}

/*}}}*/

static double breit_wigner (double m_d) /*{{{*/
{
   double dm = m_d - RESONANCE_REST_ENERGY;
   return 1.0 / (dm*dm + RESONANCE_WIDTH*RESONANCE_WIDTH);
}

/*}}}*/

static double minus_isobar_integrand (double m_d, void *x) /*{{{*/
{
   double d;
   MINUS_COMMON;
   (void) g_pi;
   d = 2.0 * bta(g_d_m) * g_d_m * bta(g_pi_i) * g_pi_i;
   if (d > 0.0)
     return breit_wigner (m_d) / d;
   else return 0.0;
}

/*}}}*/

static double plus_isobar_integrand (double m_d, void *x) /*{{{*/
{
   double d;
   PLUS_COMMON;
   (void) g_pi;
   d = 2.0 * bta(g_d_p) * g_d_p * bta(g_pi_i) * g_pi_i;
   if (d > 0.0)
     return breit_wigner (m_d) / d;
   else return 0.0;
}

/*}}}*/

static double bw_norm (double T_p) /*{{{*/
{
   double root_s, t1, t2, w;

   root_s = sqrt(TWO_PROTON_REST_ENERGY * (T_p + TWO_PROTON_REST_ENERGY));
   t1 = root_s - PROTON_REST_ENERGY - RESONANCE_REST_ENERGY;
   t2 = PROTON_REST_ENERGY + PIZERO_REST_ENERGY - RESONANCE_REST_ENERGY;
   w = RESONANCE_WIDTH;

   return w / (atan(t1/w) - atan(t2/w));
}

/*}}}*/

static int do_isobar_integral (gsl_function *f, double min, double max, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   double epsabs = 0.0;
   double epsrel = 1.e-9;
   size_t limit = MAX_QAG_SUBINTERVALS;
   double abserr;
   int status;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (f, min, max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);
   if (isnan(*val))
     *val = 0.0;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status > 0 && status != GSL_EROUND)
     fprintf (stderr, "*** isobar integral: %s\n", gsl_strerror (status));

   return status;
}

/*}}}*/

static int isobar_integral (double T_p, double T_pi, double *val) /*{{{*/
{
   gsl_function f;
   Collision_Info_Type info;
   double m_min, m_max;
   double minus_lower, minus_upper, minus_val;
   double plus_lower, plus_upper, plus_val;

   *val = 0.0;

   init_collision_info (&info, T_p, T_pi);
   f.params = &info;

   m_min = PROTON_REST_ENERGY + PIZERO_REST_ENERGY;
   m_max = sqrt(info.s) - PROTON_REST_ENERGY;

   if ((0 == bisection (&h_minus_lower, m_min, m_max, &info, &minus_lower))
       && (0 == bisection (&h_minus_upper, m_min, m_max, &info, &minus_upper)))
     {
        f.function = &minus_isobar_integrand;
        if (-1 == do_isobar_integral (&f, minus_lower, minus_upper, &minus_val))
          return -1;
     }
   else minus_val = 0.0;

   if ((0 == bisection (&h_plus_lower, m_min, m_max, &info, &plus_lower))
       && (0 == bisection (&h_plus_upper, m_min, m_max, &info, &plus_upper)))
     {
        f.function = &plus_isobar_integrand;
        if (-1 == do_isobar_integral (&f, plus_lower, plus_upper, &plus_val))
          return -1;
     }
   else plus_val = 0.0;

   *val = (minus_val + plus_val) * bw_norm (T_p) /(2 * PIZERO_REST_ENERGY);

   return 0;
}

/*}}}*/

/* beta = 1 - delta_beta(gamma)
 * for gamma > GAMMA_THRESH this should give delta_beta to a precision of
 * better than DBL_EPSILON
 */
#define GAMMA_THRESH 1.e3
static double delta_beta (double gamma) /*{{{*/
{
   double x = 1.0 / gamma / gamma;
   /* Derived from the series for sqrt(1 - x) */
   return (x/2.0)
     *(1.0 + (x/4.0)
       *(1.0 + (3*x/6.0)
         *(1.0 + (5*x/8.0)
           *(1.0 + (7*x/10.0)
             *(1.0 + (9*x/12.0))))));
}

/*}}}*/

static double lidcs_stephens_badhwar (double T_p, double T_pi, double mu, double delta_mu) /*{{{*/
{
   /* lidcs = Lorentz invariant differential cross-section */

   double m_p2, m_pi2, s, root_s, yp, ym, f, q, xx_cm, xbar;
   double E_p, E_pi, gamma_pi, p_pi, lidcs;
   double gamma_c, beta_c, gamma_p, p_p, beta_p, beta_pi;
   double /* E_p_cm, */ p_perp, p_parallel;
   double e_max_cm, p_max_cm, p_perp_gev;
   double f_p_cm, f_parallel, d_c, d_p, d_pi;

   if (fabs(mu) > 1.0 || delta_mu < 0.0 || 2.0 < delta_mu)
     return 0.0;

   m_p2 = SQR_PROTON_REST_ENERGY;
   m_pi2 = SQR_PIZERO_REST_ENERGY;

   /* collision invariant */
   s = TWO_PROTON_REST_ENERGY * (T_p + TWO_PROTON_REST_ENERGY);
   root_s = sqrt(s);

   /* Lorentz factor for center of momentum system */
   gamma_c = root_s / TWO_PROTON_REST_ENERGY;
   beta_c = bta(gamma_c);

   /* Lab frame proton energy/momentum */
   E_p = T_p + PROTON_REST_ENERGY;
   gamma_p = E_p / PROTON_REST_ENERGY;
   beta_p = bta(gamma_p);
   p_p = E_p * beta_p;

   /* Lab frame pion energy/momentum */
   E_pi = T_pi + PIZERO_REST_ENERGY;
   gamma_pi = E_pi / PIZERO_REST_ENERGY;
   beta_pi = bta(gamma_pi);
   p_pi = E_pi * beta_pi;

#if 0
   /* CM frame proton energy */
   if (gamma_c < GAMMA_THRESH && gamma_p < GAMMA_THRESH)
     f_p_cm = (1.0 - beta_c * beta_p);
   else
     {
        d_c = delta_beta (gamma_c);
        d_p = delta_beta (gamma_p);
        f_p_cm = d_c + d_p - d_c * d_p;
     }
   E_p_cm = gamma_c * E_p * f_p_cm;
#else
   (void) f_p_cm;  (void) d_p;
#endif

   /* CM frame pion momentum components */
   p_perp = p_pi * sqrt ((1.0 + mu)*delta_mu);
   if (gamma_pi < GAMMA_THRESH && gamma_c < GAMMA_THRESH)
     f_parallel = bta(gamma_pi) * mu - bta(gamma_c);
   else
     {
        d_c = delta_beta (gamma_c);
        d_pi = delta_beta (gamma_pi);
        f_parallel = - delta_mu - mu * d_pi + d_c;
     }
   p_parallel = gamma_c * E_pi * f_parallel;

   /* Max CM frame pion momentum */
   e_max_cm = 0.5*(root_s - PIZERO_MASS_FACTOR / root_s);
   p_max_cm = root_diff_sqr (e_max_cm, PIZERO_REST_ENERGY);

   /* S&B's cross-section parameterization */
   yp = 1.0 + 4*m_p2/s;
   ym = 1.0 - 4*m_p2/s;
   f = (1.0 + 23.0 / pow(T_p/GEV, 2.6)) * (ym*ym);
   p_perp_gev = p_perp / GEV;
   q = (6.1 + p_perp_gev * (-3.3 + 0.6 * p_perp_gev)) / sqrt(yp);

   xx_cm = p_parallel / p_max_cm;
   xbar = sqrt (xx_cm * xx_cm + (4.0/s)*(p_perp*p_perp + m_pi2));

   lidcs = 140.0 * f * pow (1.0 - xbar, q) * exp (-5.43 * p_perp_gev/yp);
   lidcs *= MILLIBARN / GEV / GEV;

   if (isnan(lidcs))
     lidcs = 0.0;

   return lidcs;
}

/*}}}*/

double pizero_lidcs (double T_p, double T_pi, double mu) /*{{{*/
{
   return lidcs_stephens_badhwar (T_p, T_pi, mu, 1.0 - mu);
}

/*}}}*/

static double angular_integrand (double mu, void *x) /*{{{*/
{
   Collision_Info_Type *info = (Collision_Info_Type *)x;
   return lidcs_stephens_badhwar (info->T_p, info->T_pi, mu, 1.0 - mu);
}

/*}}}*/

static double delta_mu_angular_integrand (double delta_mu, void *x) /*{{{*/
{
   Collision_Info_Type *info = (Collision_Info_Type *)x;
   return lidcs_stephens_badhwar (info->T_p, info->T_pi, 1.0-delta_mu, delta_mu);
}

/*}}}*/

static int angular_integral (double T_p, double T_pi, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   Collision_Info_Type info;
   double epsabs, epsrel, abserr;
   double gamma_c, E_pi, gamma_pi, e_max_cm, root_s;
   double mu_min, mu_max, v, v_eps, delta_mu_min;
   size_t limit;
   int status;

   *val = 0.0;

   init_collision_info (&info, T_p, T_pi);
   root_s = sqrt(info.s);

   /* Lorentz factor for center of momentum system */
   gamma_c = root_s /TWO_PROTON_REST_ENERGY;

   /* Lab frame pion energy/momentum */
   E_pi = T_pi + PIZERO_REST_ENERGY;
   gamma_pi = E_pi / PIZERO_REST_ENERGY;

   /* Max CM frame pion momentum */
   e_max_cm = 0.5 * (root_s - PIZERO_MASS_FACTOR / root_s);

   mu_max = 1.0;
   mu_min = (1.0 - e_max_cm/(gamma_c * E_pi)) / (bta(gamma_c) * bta(gamma_pi));
   if (mu_min < -1)
     mu_min = -1.0;
   else if (1 <= mu_min)
     return 0;
   delta_mu_min = 1.0 - mu_min;

   f.function = &angular_integrand;
   f.params = &info;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

#define EPSILON  1.e-3
   if (mu_min < mu_max-EPSILON)
     {
        status = gsl_integration_qag (&f, mu_min, mu_max-EPSILON, epsabs, epsrel, limit,
                                      GSL_INTEG_GAUSS15,
                                      work, &v, &abserr);

        if (isnan(v))
          v = 0.0;
        delta_mu_min = EPSILON;
     }
   else v = 0.0;

   f.function = &delta_mu_angular_integrand;
   status = gsl_integration_qag (&f, 0.0, delta_mu_min, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, &v_eps, &abserr);

   if (isnan(v_eps))
     v_eps = 0.0;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status > 0 && status != GSL_EROUND)
     fprintf (stderr, "*** angular integral: %s\n", gsl_strerror (status));

   *val = (v + v_eps) * (2 * M_PI * E_pi * bta (gamma_pi));

   return 0;
}

/*}}}*/

static double pizero_dermer_differential_xsec (double T_p, double T_pi) /*{{{*/
{
   double xsec_isobar, xsec_scaling, xsec;
   double e1=3.0, e2=7.0;
   double tp = T_p / GEV;

   if (tp < e2)
     {
        if (-1 == isobar_integral (T_p, T_pi, &xsec_isobar))
          return 0.0;
        xsec_isobar *= pizero_dermer_total_xsec (T_p);
     }
   else xsec_isobar = 0.0;

   if (e1 < tp)
     {
        if (-1 == angular_integral (T_p, T_pi, &xsec_scaling))
          return 0.0;
     }
   else xsec_scaling = 0.0;

   if (tp < e1)
     xsec = xsec_isobar;
   else if (e2 <= tp)
     xsec = xsec_scaling;
   else
     {
        double f = (tp - e1) / (e2 - e1);
        xsec = (1.0 - f) * xsec_isobar + f * xsec_scaling;
     }

   return xsec;
}

/*}}}*/

static double pizero_blattnig_differential_xsec (double T_p, double T_pi) /*{{{*/
{
   double a, sigma, xpi;

   /* Blattnig, et al, 2000 Lab-frame differential cross-section
    * parameterization (their equation 32)
    */

   T_p /= GEV;
   T_pi /= GEV;

#if 1
   /* Cross-sections apply only for 0.3 <= T_p <= 50 GeV */
   if (T_p < 0.3 || T_p > 50.0)
     return 0.0;
#endif

   xpi = pow(T_pi, 0.2);
   a = -5.8 - 1.82/pow(T_p, 0.4) + (13.5 - 4.5/xpi)/xpi;

   /* d(sigma)/dE_pizero */
   sigma = exp(a) * (MILLIBARN / GEV);

   return sigma;
}

/*}}}*/

double pizero_differential_xsec (double T_p, double T_pi) /*{{{*/
{
   if (Pizero_Method == 2)
     return pizero_dermer_differential_xsec (T_p, T_pi);

   return pizero_blattnig_differential_xsec (T_p, T_pi);
}

/*}}}*/

static double proton_integrand (double pc, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   Particle_Type *proton = p->protons;
   double e_proton, beta, v, np, xsec;
   double T_pi, T_p, result;

   (void)(*proton->spectrum) (proton, pc, &np);

   e_proton = hypot (pc, PROTON_REST_ENERGY);
   T_p = e_proton - PROTON_REST_ENERGY;
   T_pi = p->energy - PIZERO_REST_ENERGY;

   /* cross-section per unit proton energy */
   xsec = pizero_differential_xsec (T_p, T_pi);

   beta = pc / e_proton;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;
   result = np * v * xsec;

   /* We're integrating over momentum (pc, really), not energy, so:
    * dE = d(pc) (pc/E) = d(pc) * beta */
   result *= beta;

   return result;
}
/*}}}*/

static double proton_momentum_threshold (double e_pizero) /*{{{*/
{
   double r = 0.5 * SQR_PIZERO_REST_ENERGY / SQR_PROTON_REST_ENERGY;
   double pc, x;

   /* threshold proton momentum to produce the given pi-zero */

   x = e_pizero / PROTON_REST_ENERGY;
   pc = root_diff_sqr (2*x + 1.0 + r, 1.0);

   return pc * PROTON_REST_ENERGY;
}

/*}}}*/

static struct
{
   int flag;
   double max;
}
Threshold;

static int integral_over_proton_momenta (Pizero_Type *p, double *val) /*{{{*/
{
   Particle_Type *protons = p->protons;
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double pc_min, pc_max, pc_thresh;
   size_t limit;
   int status;

   *val = 0.0;

   if (p->energy <= PIZERO_REST_ENERGY)
     return 0;

   pc_max = (*protons->momentum_max)(protons);
   pc_min = (*protons->momentum_min)(protons);

   pc_thresh = proton_momentum_threshold (p->energy);
   if (pc_thresh < pc_min)
     {
        if ((Threshold.flag == 0)
            || (pc_thresh > Threshold.max))
          {
             Threshold.max = pc_thresh;
          }
        Threshold.flag = 1;
        pc_thresh = pc_min;
     }

   f.function = &proton_integrand;
   f.params = p;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, pc_thresh, pc_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS31,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status > 0 && status != GSL_EROUND)
     fprintf (stderr, "*** pizero: proton integral:  %s\n", gsl_strerror (status));

   if (isnan(*val))
     *val = 0.0;

   return 0;
}

/*}}}*/

static int pizero_build_table (Pizero_Type *p, double epi_min, double epi_max) /*{{{*/
{
   double lg_min, lg_max, dlg, lg;
   double *x, *y;
   unsigned int i, n = PIZERO_TABLE_SIZE;
   int status;

   if (NULL == (x = malloc (2*n * sizeof(double))))
     return -1;
   y = x + n;

   /* FIXME -- an adaptive grid might be better ... */

   lg_min = log10(epi_min);
   lg_max = log10(epi_max);
   dlg = (lg_max - lg_min) / (n - 1);

   for (i = 0; i < n; i++)
     {
        double xx, yy;
        lg = lg_min + dlg * i;
        xx = pow(10.0, lg);
        p->energy = xx;
        if (-1 == integral_over_proton_momenta(p, &yy))
          {
             free(x);
             return -1;
          }
        x[i] = lg;
        y[i] = (yy > 0.0) ? log10(yy) : DBL_MIN_10_EXP;
     }

   status = pizero_spline_table (p->client_data, x, y, n);
   free(x);

   return status;
}

/*}}}*/

static double pizero_aa_total_xsec (double T_p) /*{{{*/
{
   double sigma = 0.0;

   /* Aharonian and Atoyan (2000) */
   if (T_p > GEV)
     sigma = MILLIBARN * 30.0 * (0.95 + 0.06 * log (T_p/GEV));

   return sigma;
}

/*}}}*/

static int delta_function_approximation (Pizero_Type *p, double *val) /*{{{*/
{
   /* delta-function approximation: see Aharonian and Atoyan (2000) */
   Particle_Type *protons = p->protons;
   double proton_pc, T_p, kappa;
   double np, beta, sigma, v, x, pc;

   *val = 0.0;

   /* kappa = mean fraction of proton kinetic energy transferred
    *         to the secondary meson per collision. */
   kappa = 0.17;
   T_p = p->energy / kappa;

#if 1
   sigma = pizero_aa_total_xsec (T_p);
#else
   sigma = pizero_dermer_total_xsec (T_p);
#endif

   if (sigma == 0.0)
     return 0;

   x = T_p / PROTON_REST_ENERGY;
   proton_pc = x * sqrt (1.0 + 2.0/x);
   beta = proton_pc / (1.0 + x);
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   pc = proton_pc * PROTON_REST_ENERGY;
   (void)(*protons->spectrum)(protons, pc, &np);

   *val = np * v * sigma /kappa;

   return 0;
}

/*}}}*/

int pizero_distribution (Pizero_Type *p, double *val) /*{{{*/
{
   if (Pizero_Method != 0)
     return integral_over_proton_momenta (p, val);

   return delta_function_approximation (p, val);
}

/*}}}*/

static double pizero_integrand (double w, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   double q, s, e_pizero;
   int status;

   /* use log-variable in integral */
   w = exp(w);

   /* changed integration variable to remove
    * singularity at pizero rest energy
    */

   e_pizero = (w*w + 1.0) * PIZERO_REST_ENERGY;

   if (p->interpolate == 2)
     {
        status = pizero_interp_pizero_distribution (p->client_data, e_pizero, &q);
     }
   else
     {
        p->energy = e_pizero;
        status = pizero_distribution (p, &q);
     }

   s = q / w / sqrt (1.0 + 2.0/w/w);

   /* from change of variables removing singularity */
   s *= 2.0;

   /* from change to log integration interval */
   s *= w;

   return status ? 0.0 : s;
}

/*}}}*/

static int pizero_min_energy (double photon_energy, double *e_pizero) /*{{{*/
{
   double x;

   if (photon_energy <= 0.0)
     return -1;

   x = 0.5 * PIZERO_REST_ENERGY / photon_energy;
   *e_pizero = photon_energy * (1.0 + x*x);

   return 0;
}

/*}}}*/

static int pizero_max_energy (Particle_Type *protons, double *e_pizero) /*{{{*/
{
   double gp_max, x;
   /* double e_pizero_cm, gamma_cm, pc_max, s; */

   /* If we use a constant here, the spectrum model
    * satisfies the recurrence relation more accurately
    * near photon_energy = 0.5*pizero_rest_mass
    */
#if 0
   pc_max = (*protons->momentum_max)(protons) / PROTON_REST_ENERGY;
   gp_max = hypot (pc_max, 1.0);
#else
   (void) protons;
   gp_max = 1.e12;
#endif

#if 0
   /* lab frame max energy */
   gamma_cm = sqrt (0.5*(gp_max + 1));

   /* CM frame energy:  See Blattnig et al (2000), Appendix A */
   s = 2.0*SQR_PROTON_REST_ENERGY * (1.0 + gp_max);
   e_pizero_cm = 0.5 * sqrt(s) * (1.0 - PIZERO_MASS_FACTOR / s);
   *e_pizero = gamma_cm * e_pizero_cm;
#else
   x = 0.5 * (gp_max - 1.0) + MASS_RATIO * MASS_RATIO;
   *e_pizero = PROTON_REST_ENERGY * x;
#endif

   return 0;
}

/*}}}*/

static int integral_over_pizero_energies (Pizero_Type *p, double photon_energy, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double epi_min, epi_max, w_min, w_max;
   size_t limit;
   int status;

   *val = 0.0;

   if (-1 == pizero_min_energy (photon_energy, &epi_min))
     return -1;
   if (-1 == pizero_max_energy (p->protons, &epi_max))
     return -1;

   if (epi_min >= epi_max)
     return 0;

   /* FIXME(?) can re-use table as long as proton distribution is fixed */
   if (p->interpolate == 1)
     {
        double min_epi_min;
        if (-1 == pizero_min_energy (MIN_PHOTON_ENERGY, &min_epi_min))
          return -1;
        if (-1 == pizero_build_table (p, min_epi_min, epi_max))
          return -1;
        p->interpolate = 2;
     }

   f.function = &pizero_integrand;
   f.params = p;
   epsabs = 0.0;
   epsrel = 1.e-10;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   /* change variables to avoid singularity at lower integration limit */
   w_min = sqrt (epi_min/PIZERO_REST_ENERGY - 1.0);
   w_max = sqrt (epi_max/PIZERO_REST_ENERGY - 1.0);

   status = gsl_integration_qag (&f, log(w_min), log(w_max), epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (isnan(*val))
     *val = 0.0;

   if (status > 0)
     {
        /* if (status != GSL_EROUND) */
        if ((fabs(*val) > 0) && (abserr < Pizero_Epsrel * fabs(*val)))
          {             
             Nonthermal_Error_Type e;
             e.error_msg = gsl_strerror (status);
             e.value = *val;
             e.estimated_abserr = abserr;
             e.allowed_abserr = Pizero_Epsrel * fabs (*val);             
             e.allowed_epsrel = Pizero_Epsrel;
             nonthermal_error_hook (&e, p->protons, __FILE__, __LINE__);
          }
     }

   /* two photons per pion */
   *val *= 2.0;

   return 0;
}

/*}}}*/

int pizero_decay (void *x, double photon_energy, double *emissivity)
{
   Pizero_Type *p = (Pizero_Type *)x;
   int status;

   *emissivity = 0.0;
   photon_energy *= GSL_CONST_CGSM_ELECTRON_VOLT;

   Threshold.flag = 0;
   status = integral_over_pizero_energies (p, photon_energy, emissivity);

   if (Threshold.flag)
     {
#if 0
        fprintf (stderr, "WARNING: pizero: underestimated photon rate < %g GeV;  needed %g GeV proton threshold\n",
                 photon_energy / GEV, Threshold.max / GEV);
#endif
     }

   return status;
}

