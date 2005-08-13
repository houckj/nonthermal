/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Feb 2003
 *
 * Neutral pion production cross-sections from
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
#define PIZERO_MASS_FACTOR      (4*SQR_PROTON_REST_ENERGY - SQR_PIZERO_REST_ENERGY)
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
   double pc = sqrt (ep*ep - SQR_PROTON_REST_ENERGY);
   double s, sigma = 0.0;

   pc /= GEV;

   /* Dermer, C.D., 1986, A&A, 157, 223. */
   if (pc < 0.78)
     sigma = 0.0;
   else if (0.78 <= pc && pc < 0.96)
     {
        double mx2, x, eta2, eta4;
        s = 2*PROTON_REST_ENERGY * (T_p + 2*PROTON_REST_ENERGY);
        mx2 = 4*SQR_PROTON_REST_ENERGY;
        x = s - SQR_PIZERO_REST_ENERGY - mx2;
        eta2 = ((x*x - 4*SQR_PIZERO_REST_ENERGY*mx2)
                /(4*SQR_PIZERO_REST_ENERGY*s));
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
        s = 2*PROTON_REST_ENERGY * (T_p + 2*PROTON_REST_ENERGY);
        /* s /= GEV; */
        sigma = 163.0 * pow(s/1876.0, 0.21);
     }

   return sigma * MILLIBARN;
}

/*}}}*/

static double bta (double g) /*{{{*/
{
   if (g > 1.0)
     return sqrt((g+1.0)*(g-1.0))/g;
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
   double g_c, g_d_cm, g_pi_i, g_pi, root_s; \
   double m_d2, m_p2, m_pi2, T_pi, s; \
   T_pi = info->T_pi; \
   s = info->s; \
   root_s = sqrt(s); \
   m_d2 = m_d * m_d;  \
   m_p2 = SQR_PROTON_REST_ENERGY; \
   m_pi2 = SQR_PIZERO_REST_ENERGY; \
   g_pi = 1.0 + T_pi / PIZERO_REST_ENERGY; \
   g_c = root_s / TWO_PROTON_REST_ENERGY; \
   g_d_cm = (s + m_d2 - m_pi2) / (2 * root_s * m_d); \
   g_pi_i = (m_d2 + m_pi2 - m_p2) / (2 * PIZERO_REST_ENERGY * m_d);

#define PLUS_COMMON \
   double g_d_p; \
   COMMON_FACTORS \
   g_d_p = gamma_fun (g_c, g_d_cm, +1);

#define MINUS_COMMON \
   double g_d_m; \
   COMMON_FACTORS; \
   g_d_m = gamma_fun (g_c, g_d_cm, -1);

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

#if 0
   status = gsl_integration_qag (f, min, max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);
#else
     {
        size_t neval;
        status = gsl_integration_qng (f, min, max, epsabs, epsrel, val, &abserr,
                                      &neval);
     }
#endif

   if (isnan(*val))
     *val = 0.0;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
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
   double E_p_cm, p_perp, p_parallel;
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
   p_max_cm = sqrt (e_max_cm * e_max_cm - m_pi2);

   /* S&B's cross-section parameterization */
   yp = 1.0 + 4*m_p2/s;
   ym = 1.0 - 4*m_p2/s;
   f = (1.0 + 23.0 / pow(E_p_cm/GEV, 2.6)) * (ym*ym);
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

   if (status)
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
   double eproton_delta, proton_pc, T_p;
   double np, beta, sigma, v;

   /* kappa = mean fraction of proton kinetic energy transferred
    *         to the secondary meson per collision. */
   double kappa = 0.17;

   eproton_delta = PROTON_REST_ENERGY + p->energy /kappa;
   proton_pc = sqrt (eproton_delta * eproton_delta - SQR_PROTON_REST_ENERGY);
   (void)(*protons->spectrum)(protons, proton_pc, &np);

   beta = proton_pc / eproton_delta;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   T_p = eproton_delta - PROTON_REST_ENERGY;
#if 1
   sigma = pizero_aa_total_xsec (T_p);
#else
   sigma = pizero_dermer_total_xsec (T_p);
#endif

   *val = np * v * sigma /kappa;

   return 0;
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

#if 0
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

   e_proton = sqrt (pc*pc + SQR_PROTON_REST_ENERGY);
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
   double eproton_thresh, pc;

   /* threshold proton momentum to produce the given pi-zero */

   eproton_thresh = 2*e_pizero + PROTON_REST_ENERGY
     + 0.5 * SQR_PIZERO_REST_ENERGY / PROTON_REST_ENERGY;
   pc = sqrt (eproton_thresh*eproton_thresh - SQR_PROTON_REST_ENERGY);

   return pc;
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
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** pizero: proton integral:  %s\n", gsl_strerror (status));

   if (isnan(*val))
     *val = 0.0;

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

static double pizero_integrand (double e_pizero, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   double q, s, pc_pizero;
   int status;

   if (e_pizero < PIZERO_REST_ENERGY)
     return 0.0;

   if (p->interpolate == 2)
     {
        status = pizero_interp_pizero_distribution (p->client_data, e_pizero, &q);
     }
   else
     {
        p->energy = e_pizero;
        status = pizero_distribution (p, &q);
     }

   pc_pizero = sqrt (e_pizero * e_pizero - SQR_PIZERO_REST_ENERGY);
   s = q / pc_pizero;

   return status ? 0.0 : s;
}

/*}}}*/

static int pizero_min_energy (double photon_energy, double *e_pizero) /*{{{*/
{
   if (photon_energy > 0)
     *e_pizero = photon_energy + 0.25 * SQR_PIZERO_REST_ENERGY / photon_energy;
   else return -1;

   return 0;
}

/*}}}*/

static int pizero_max_energy (Particle_Type *protons, double *e_pizero) /*{{{*/
{
   double pc_max, ep_max, root_s, e_pizero_cm, gamma_cm;

   pc_max = (*protons->momentum_max)(protons);
   ep_max = sqrt (pc_max*pc_max + SQR_PROTON_REST_ENERGY);
   root_s = sqrt (2*PROTON_REST_ENERGY*(PROTON_REST_ENERGY + ep_max));

   /* CM frame energy:  See Blattnig et al (2000), Appendix A */
   e_pizero_cm = 0.5 * (root_s - PIZERO_MASS_FACTOR / root_s);

   /* lab frame max energy */
   gamma_cm = sqrt (0.5*(ep_max/PROTON_REST_ENERGY + 1));
   *e_pizero = gamma_cm * e_pizero_cm;

   return 0;
}

/*}}}*/

static int integral_over_pizero_energies (Pizero_Type *p, double photon_energy, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double epi_min, epi_max;
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
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, epi_min, epi_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** pizero: pizero integral: %s\n", gsl_strerror (status));

   if (isnan(*val))
     *val = 0.0;

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

