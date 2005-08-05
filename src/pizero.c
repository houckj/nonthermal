/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Feb 2003
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

#define SQR_PROTON_REST_ENERGY  (PROTON_REST_ENERGY * PROTON_REST_ENERGY)
#define SQR_PIZERO_REST_ENERGY  (PIZERO_REST_ENERGY * PIZERO_REST_ENERGY)
#define PIZERO_MASS_FACTOR      (SQR_PIZERO_REST_ENERGY - 4*SQR_PROTON_REST_ENERGY)
#define SQR_PIZERO_MASS_FACTOR  (PIZERO_MASS_FACTOR * PIZERO_MASS_FACTOR)

/* Placing MIN_PHOTON_ENERGY just above the threshold
 * eliminates the sharp peak in the photon spectrum near 70 MeV(?)
 */
#define MIN_PHOTON_ENERGY (0.5*PIZERO_REST_ENERGY *2)

int Pizero_Use_Dermer_Xsec = 0;
double Pizero_Approx_Min_Energy = 0.0; /* 50.0; */  /* GeV */

/* -*- mode: C; mode: fold -*- */
/* Dermer 1986, A&A 157, 223 */

#define H(x,a,b)  (((a) <= (x)) && ((x) < (b)))

/* \Delta_{3/2}(1232) resonance width [GeV]
 * Particle Data Group (1984)
 */
#define RESONANCE_WIDTH         ((0.115 * GEV)*0.5)
#define RESONANCE_REST_ENERGY   (1233.0 * MEV)

#define TWO_PROTON_REST_ENERGY   (2 * PROTON_REST_ENERGY);

static double pizero_dermer_total_xsec (double proton_kinetic) /*{{{*/
{
   double ep = proton_kinetic + PROTON_REST_ENERGY;
   double pc = sqrt (ep*ep - SQR_PROTON_REST_ENERGY);
   double s, sigma = 0.0;

   pc /= GEV;

   /* Dermer, C.D., 1986, A&A, 157, 223. */
   if (pc < 0.78)
     sigma = 0.0;
   else if (0.78 <= pc && pc < 0.96)
     {
        double mx2, x, eta2, eta4;
        s = 2*PROTON_REST_ENERGY * (proton_kinetic + 2*PROTON_REST_ENERGY);
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
        s = 2*PROTON_REST_ENERGY * (proton_kinetic + 2*PROTON_REST_ENERGY);
        sigma = 163.0 * pow(s/1876.0, 0.21);
     }

   return sigma * MILLIBARN;
}

/*}}}*/

static double beta_fun (double g) /*{{{*/
{
   return sqrt ((g + 1.0)*(g - 1.0)) / g;
}

/*}}}*/

static double gfun (double g1, double g2, double sgn) /*{{{*/
{
   double b1 = beta_fun (g1);
   double b2 = beta_fun (g2);

   return g1 * g2 * (1.0 + sgn * b1 * b2);
}

/*}}}*/

static double pfun (double g) /*{{{*/
{
   return g * beta_fun (g);
}

/*}}}*/

static int isobar_spectrum (double T_p, double T_pi, double m_d) /*{{{*/
{
   double two_mp, two_mpi, m_d2, m_pi2, m_p2;
   double s, root_s, gamma_c, gamma_d_cm;
   double a_p, b_p, a_m, b_m;
   double gamma_d_p, gamma_d_m, gamma_pi_i, gamma_pi, f;

   m_d2 = m_d * m_d;
   m_p2 = SQR_PROTON_REST_ENERGY;
   m_pi2 = SQR_PIZERO_REST_ENERGY;

   two_mpi = 2 * PIZERO_REST_ENERGY;
   two_mp = TWO_PROTON_REST_ENERGY;
   s = two_mp * (T_p + two_mp);
   root_s = sqrt(s);

   gamma_c = root_s / two_mp;
   gamma_d_cm = (s + m_d2 - m_pi2) / (2*root_s * m_d);
   gamma_pi_i = (m_d2 + m_pi2 - m_p2) / (two_mpi * m_d);
   gamma_d_p = gfun (gamma_c, gamma_d_cm, +1.0);
   gamma_d_m = gfun (gamma_c, gamma_d_cm, -1.0);

   a_p = gfun (gamma_d_p, gamma_pi_i, -1.0);
   b_p = gfun (gamma_d_p, gamma_pi_i, +1.0);

   a_m = gfun (gamma_d_m, gamma_pi_i, -1.0);
   b_m = gfun (gamma_d_m, gamma_pi_i, +1.0);

   gamma_pi = 1.0 + T_pi / PIZERO_REST_ENERGY;

   f = 0.0;

   if (H(gamma_pi, a_p, b_p))
     f += 1.0 / (2.0 * pfun (gamma_d_p) * pfun (gamma_pi_i));

   if (H(gamma_pi, a_m, b_m))
     f += 1.0 / (2.0 * pfun (gamma_d_m) * pfun (gamma_pi_i));

   f /= two_mpi;

   return f;
}

/*}}}*/

static double breit_wigner (double m_d) /*{{{*/
{
   double w = RESONANCE_WIDTH;
   double dm = m_d - RESONANCE_REST_ENERGY;
   return (w/M_PI) / (dm*dm + w*w);
}

/*}}}*/

static double bw_norm (double T_p) /*{{{*/
{
   double root_s, t1, t2, two_mp, x;

   two_mp = TWO_PROTON_REST_ENERGY;
   x = RESONANCE_WIDTH;

   root_s = sqrt(two_mp * (T_p + two_mp));
   t1 = root_s - PROTON_REST_ENERGY - RESONANCE_REST_ENERGY;

   t2 = PROTON_REST_ENERGY + PIZERO_REST_ENERGY - RESONANCE_REST_ENERGY;

   return M_PI / (atan2 (t1,x) - atan2 (t2,x));
}

/*}}}*/

static double lidcs_stephens_badhwar (double T_p, double T_pi, double mu) /*{{{*/
{
   double m_p2, m_pi2, two_mp, s, root_s, y, f, q, xx_cm, xbar;
   double E_p, E_pi, gamma_pi, p_pi, lidcs;
   double gamma_c, beta_c, gamma_p, p_p, beta_pi;
   double E_p_cm, p_perp, p_parallel;
   double e_max_cm, p_max_cm;

   m_p2 = PROTON_REST_ENERGY * PROTON_REST_ENERGY;
   m_pi2 = PIZERO_REST_ENERGY * PIZERO_REST_ENERGY;
   two_mp = TWO_PROTON_REST_ENERGY;

   /* collision invariant */
   s = two_mp * (T_p + two_mp);
   root_s = sqrt(s);

   /* Lorentz factor for center of momentum system */
   gamma_c = root_s/two_mp;
   beta_c = beta_fun(gamma_c);

   /* Lab frame proton energy/momentum */
   E_p = T_p + PROTON_REST_ENERGY;
   gamma_p = E_p / PROTON_REST_ENERGY;
   p_p = E_p * beta_fun (gamma_p);

   /* Lab frame pion energy/momentum */
   E_pi = T_pi + PIZERO_REST_ENERGY;
   gamma_pi = E_pi / PIZERO_REST_ENERGY;
   beta_pi = beta_fun (gamma_pi);
   p_pi = E_pi * beta_pi;

   /* CM frame proton energy */
   E_p_cm = gamma_c * (E_p - beta_c * p_p);

   /* CM frame pion momentum components */
   p_perp = p_pi * sqrt (1.0 - mu*mu);
   p_parallel = gamma_c * (p_pi * mu - beta_c * E_pi);

   /* Max CM frame pion momentum */
   e_max_cm = (s + PIZERO_MASS_FACTOR) / (2*root_s);
   p_max_cm = sqrt (e_max_cm * e_max_cm - m_pi2);

   /* S&B cross-section parameterization */
   y = 1.0 + 4*m_p2/s;
   f = (1.0 + 23.0 / pow(E_p_cm, 2.6)) / (y*y);
   q = (6.1 + p_perp * (-3.3 + 0.6 * p_perp)) / sqrt(y);

   xx_cm = p_parallel / p_max_cm;
   xbar = sqrt (xx_cm * xx_cm + (4.0/s)*(p_perp*p_perp + m_pi2));

   lidcs = 140.0 * f * pow (1.0 - xbar, q) * exp (-5.43 * p_perp/y);

   return lidcs;
}

/*}}}*/

typedef struct
{
   double T_pi;
   double T_p;
   double s;
}
Collision_Info_Type;

static void init_collision_info (Collision_Info_Type *info, double T_p, double T_pi) /*{{{*/
{
   double two_mp = TWO_PROTON_REST_ENERGY;
   info->s = two_mp * (T_p + two_mp);
   info->T_pi = T_pi;
   info->T_p = T_p;
}

/*}}}*/

static double isobar_integrand (double m, void *x) /*{{{*/
{
   Collision_Info_Type *info = (Collision_Info_Type *)x;
   double f, b;

   f = isobar_spectrum (info->T_p, info->T_pi, m);
   b = breit_wigner (m);

   return b * f;
}

/*}}}*/

static int isobar_integral (double T_p, double T_pi, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   Collision_Info_Type info;
   double epsabs, epsrel, abserr;
   double m_min, m_max;
   size_t limit;
   int status;

   *val = 0.0;

   init_collision_info (&info, T_p, T_pi);

   m_min = PROTON_REST_ENERGY - PIZERO_REST_ENERGY;
   m_max = sqrt(info.s) - PROTON_REST_ENERGY;

   f.function = &isobar_integrand;
   f.params = &info;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

#if 0
   status = gsl_integration_qag (&f, m_min, m_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);
#else
     {
        size_t neval;
        status = gsl_integration_qng (&f, m_min, m_max, epsabs, epsrel, val, &abserr,
                                      &neval);
     }
#endif

   /* Check for NaN */
   if (*val != *val)
     *val = 0.0;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     fprintf (stderr, "*** isobar integral: %s\n", gsl_strerror (status));

   *val *= bw_norm (T_p);

   return 0;
}

/*}}}*/

static double angular_integrand (double mu, void *x) /*{{{*/
{
   Collision_Info_Type *info = (Collision_Info_Type *)x;
   return lidcs_stephens_badhwar (info->T_p, info->T_pi, mu);
}

/*}}}*/

static int angular_integral (double T_p, double T_pi, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   Collision_Info_Type info;
   double epsabs, epsrel, abserr;
   double gamma_c, beta_c;
   double E_pi, gamma_pi, beta_pi, p_pi, e_max_cm;
   double mu_min, mu_max;
   size_t limit;
   int status;

   *val = 0.0;

   init_collision_info (&info, T_p, T_pi);

   /* Lorentz factor for center of momentum system */
   gamma_c = sqrt(info.s) /TWO_PROTON_REST_ENERGY;
   beta_c = beta_fun(gamma_c);

   /* Lab frame pion energy/momentum */
   E_pi = T_pi + PIZERO_REST_ENERGY;
   gamma_pi = E_pi / PIZERO_REST_ENERGY;
   beta_pi = beta_fun (gamma_pi);
   p_pi = E_pi * beta_pi;

   /* Max CM frame pion momentum */
   e_max_cm = (info.s + PIZERO_MASS_FACTOR) / (2*sqrt(info.s));

   mu_max = 1.0;
   mu_min = (gamma_c * E_pi - e_max_cm) / (beta_c * gamma_c * p_pi);
   if (mu_min < -1 || 1 <= mu_min)
     return 0;

   f.function = &angular_integrand;
   f.params = &info;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

#if 0
   status = gsl_integration_qag (&f, mu_min, mu_max, epsabs, epsrel, limit,
                                 GSL_INTEG_GAUSS15,
                                 work, val, &abserr);
#else
     {
        size_t neval;
        status = gsl_integration_qng (&f, mu_min, mu_max, epsabs, epsrel, val, &abserr,
                                      &neval);
     }
#endif

   /* Check for NaN */
   if (*val != *val)
     *val = 0.0;

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

#if 0 /* FIXME */
   if (status)
     fprintf (stderr, "*** angular integral: %s\n", gsl_strerror (status));
#endif

   *val *= (2*M_PI * p_pi) * (MILLIBARN/GEV);

   return 0;
}

/*}}}*/

static double pizero_dermer_differential_xsec (double T_p, double T_pi) /*{{{*/
{
   double xsec_isobar, xsec_scaling, xsec, e_pizero, Tp_thresh;
   double e1=3.0, e2=7.0;
   double tp = T_p / GEV;

#if 1
   if (T_pi <= 0.0)
     return 0.0;
   e_pizero = T_pi + PIZERO_REST_ENERGY;
   Tp_thresh = 2*e_pizero + 0.5*SQR_PIZERO_REST_ENERGY / PROTON_REST_ENERGY;
   if (T_p < Tp_thresh)
     return 0.0;
#endif

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

/*
 * Neutral pion production cross-sections from
 *   Blattnig et al, NASA Technical Report, NASA/TP-2000-210640
 *   Aharonian and Atoyan, 2000, A&A, 362, 937
 *
 * See also
 *   Domingo-Santamaria & Torres, 2005, astro-ph/0506240
 *   Dermer, C.D., 1986, A&A, 157, 223.
 */

static double pizero_aa_total_xsec (double proton_kinetic) /*{{{*/
{
   double sigma;

   /* Aharonian and Atoyan (2000) */
   if (proton_kinetic > GEV)
     sigma = 30.0 * (0.95 + 0.06 * log (proton_kinetic/GEV));
   else return 0.0;

   return sigma * MILLIBARN;
}

/*}}}*/

/* delta-function approximation: see Aharonian and Atoyan (2000) */
static double delta_function_approximation (Pizero_Type *p) /*{{{*/
{
   Particle_Type *protons = p->protons;
   double eproton_delta, proton_pc, proton_kinetic;
   double np, beta, sigma, v, val;

   /* kappa = mean fraction of proton kinetic energy transferred
    *         to the secondary meson per collision. */
   double kappa = 0.17;

   eproton_delta = PROTON_REST_ENERGY + p->energy /kappa;
   proton_pc = sqrt (eproton_delta * eproton_delta - SQR_PROTON_REST_ENERGY);
   (void)(*protons->spectrum)(protons, proton_pc, &np);

   beta = proton_pc / eproton_delta;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   proton_kinetic = eproton_delta - PROTON_REST_ENERGY;
   sigma = pizero_aa_total_xsec (proton_kinetic);

   val = np * v * sigma /kappa;

   return val;
}

/*}}}*/

static double pizero_blattnig_differential_xsec (double proton_kinetic, double pizero_kinetic) /*{{{*/
{
   double a, sigma, xpi;

   /* Blattnig, et al, 2000 Lab-frame differential cross-section
    * parameterization (their equation 32)
    */

   proton_kinetic /= GEV;
   pizero_kinetic /= GEV;

   xpi = pow(pizero_kinetic, 0.2);
   a = -5.8 - 1.82/pow(proton_kinetic, 0.4) + (13.5 - 4.5/xpi)/xpi;

   /* d(sigma)/dE_pizero */
   sigma = exp(a) * (MILLIBARN / GEV);

   return sigma;
}

/*}}}*/

static double pizero_differential_xsec (double T_p, double T_pi)
{
   if (Pizero_Use_Dermer_Xsec != 0)
     return pizero_dermer_differential_xsec (T_p, T_pi);

   return pizero_blattnig_differential_xsec (T_p, T_pi);
}

double pizero_diff_xsec (double T_p, double T_pi)
{
   return pizero_differential_xsec (T_p, T_pi);
}

static double proton_integrand (double pc, void *x) /*{{{*/
{
   Pizero_Type *p = (Pizero_Type *)x;
   Particle_Type *proton = p->protons;
   double e_proton, beta, v, np, xsec;
   double pizero_kinetic, proton_kinetic, result;

   (void)(*proton->spectrum) (proton, pc, &np);

   e_proton = sqrt (pc*pc + SQR_PROTON_REST_ENERGY);
   proton_kinetic = e_proton - PROTON_REST_ENERGY;
   pizero_kinetic = p->energy - PIZERO_REST_ENERGY;

   /* cross-section per unit proton energy */
   xsec = pizero_differential_xsec (proton_kinetic, pizero_kinetic);

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

   if ((Pizero_Use_Dermer_Xsec != 0)
       && (p->energy > Pizero_Approx_Min_Energy * GEV))
     {
        *val = delta_function_approximation (p);
        return 0;
     }

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

   /* check for NaN */
   if (*val != *val)
     *val = 0.0;

   return 0;
}

/*}}}*/

int pizero_distribution (Pizero_Type *p, double *val)
{
   return integral_over_proton_momenta (p, val);
}

static int _pizero_integrand (Pizero_Type *p, double e_pizero, double *s) /*{{{*/
{
   double pc_pizero, q;

   *s = 0.0;
   if (e_pizero < PIZERO_REST_ENERGY)
     return 0;

   p->energy = e_pizero;

   if (-1 == integral_over_proton_momenta (p, &q))
     return -1;

   pc_pizero = sqrt (e_pizero * e_pizero - SQR_PIZERO_REST_ENERGY);
   *s = q / pc_pizero;

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
        if (-1 == _pizero_integrand (p, xx, &yy))
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
   double s;
   int status;

   if (p->interpolate)
     {
        status = pizero_interp_pizero_integral (p->client_data, e_pizero, &s);
     }
   else
     {
        status = _pizero_integrand (p, e_pizero, &s);
     }

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
   e_pizero_cm = 0.5 * (root_s + PIZERO_MASS_FACTOR / root_s);

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

   /* Check for NaN */
   if (*val != *val)
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

