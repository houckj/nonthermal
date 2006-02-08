/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Oct 2004
 *
 * Implementation based on
 *   Haug, E., 1975, Z. Naturforsch. 30a, 1099
 *   Dermer, C.D., 1984, ApJ, 280, 328
 *   Alexanian, M., 1968, Phys. Rev., 165, 253.
 *   Heitler, W., 1953, "The Quantum Theory of Radiation",
 *         3rd edition (Dover), p 245
 *   Baring, M.G., et al, 1999, ApJ, 513, 311
*/

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "_nonthermal.h"
#include "ntbrems.h"
#include "ntb_table.h"

/* electron-ion bremsstrahlung */
double _ntb_ei_sigma (double electron_kinetic_energy, double photon_energy) /*{{{*/
{
   double k, g, g0, b, b0, sigma, Z=1.0;
   double r0, L, a, a0, phi_bar;
   double g0g, b0b, g02, g03, g2, g3, b02, b03, b2, b3, t1, t2, x;
   double xi, xi0, elwert;

   /* g0 = initial e- gamma   b0 = initial e- beta
    * g  = final e- gamma     b  = final e- beta
    * k  = photon energy
    * units me*c^2 = 1
    */

   k  = photon_energy;
   g0 = electron_kinetic_energy + 1.0;
   /* energy conservation */
   g = g0 - k;

   if ((g0 <= 1.0) || (g <= 1.0))
     return 0.0;

   g02 = g0*g0;
    g2 =  g*g;

   b  = sqrt (1.0 - 1.0/g2);
   b0 = sqrt (1.0 - 1.0/g02);

   r0 = ELECTRON_RADIUS;
   phi_bar = (Z*Z) * (r0*r0) * GSL_CONST_NUM_FINE_STRUCTURE;

   g0g = g0*g;
   b0b = b0*b;

   g03 = g02*g0; b02 = b0*b0; b03 = b02*b0;
    g3 =  g2*g;   b2 =  b*b;   b3 =  b2*b;

   a0 = 2.0 * log (g0 * (1.0 + b0));
   a  = 2.0 * log (g  * (1.0 + b));
   L  = 2.0 * log ((g0g + g0g*b0b - 1.0) / k);

   t1 = 2.0*g0g*(g02*b02 + g2*b2)/(g02*b02*g2*b2);
   t2 = a0*g/(g03*b03) + a*g0/(g3*b3) - a0*a/(g0g*b0b);
   x  = ((8.0/3.0)/b0b + k*k*(1.0 + b02*b2)/(g0g*b03*b3)
         + (0.5*k/(g0g*b0b)) * ((g + g0*b02)*a0/(g02*b03)
                              - (g0 + g*b2)*a/(g2*b3) + 2.0*k/(g0g*b02*b2)));
   sigma = (phi_bar/k) * ((g*b)/(g0*b0)) * (4.0/3.0 - t1 + t2 + L*x);

   /* Elwert correction to low-energy cross-section: */
#define TWO_PI_ALPHA (2.0 * M_PI * GSL_CONST_NUM_FINE_STRUCTURE)
   xi  = TWO_PI_ALPHA * Z / b;
   xi0 = TWO_PI_ALPHA * Z / b0;
   elwert = (xi/xi0) * (1.0 - exp(-xi0))/(1.0 - exp(-xi));
   sigma *= elwert;

   return sigma;
}

/*}}}*/

/* Electron-electron bremsstrahlung, diff. cross section in the cms
 * [barn / (keV * sr)]
 *
 *   een = incident electron kinetic energy [keV]
 *   pen = photon energy [keV]
 *    ct = cos (theta), theta = photon angle
 */

static double eebrems_diff_cms (double een, double pen, double ct) /*{{{*/
{
    /* Initialized data */

   double pi = M_PI;
   double alpha = GSL_CONST_NUM_FINE_STRUCTURE;
   double mcq = ELECTRON_REST_ENERGY / KEV;
   double ar0 = GSL_CONST_NUM_FINE_STRUCTURE*ELECTRON_RADIUS*ELECTRON_RADIUS / BARN;

   int m;
   double kmax, wiqc, zpia, f, k, l, s, w, x;
   double a1, a2, e1, g0, g1, g2, g3, g4, g5, h1, h2, h3, h4, h5,
     h6, h7, h8, l1, l2, l3, l4, p1, r1, r2, s1[2], w2, w4, x1, x2, x3,
     e11, /* ct,*/ r44, x12, kq, pm, rk, w44, ro;
   double rq, wq, wr, xq, xr, ww, p1q, rq2, rq4, wq2, rx1, rx2,
     rx3, wq4, wr1, x1q, x2q, x3q, pct, wiq, wwq, wq2q, wq4q;
   
   /* JCH added this to enforce energy conservation */
   if (pen > een)
     return 0.0;

   e11 = een / mcq;
   e1 = e11 + 1.0;
   p1q = e11 * (e11 + 2.0);
   kmax = p1q / e1;
   k = pen / mcq;

   if (k > kmax)
     return 0.0;

   p1 = sqrt(p1q);
   w = e1 + e1;
   wq = w * w;
   wq2 = wq - 2.0;
   wq2q = wq2 * wq2;
   wq4 = wq - 4.0;
   wq4q = wq4 * wq4;
   w44 = p1q;
   ww = p1 + p1;
   wwq = wq * wq4;

   zpia = pi * 2.0 * alpha;
   a1 = zpia * wq2 / (w * ww);

   pm = kmax * mcq;

   kq = k * k;
   rk = 1.0 / k;
   rq = wq - e1 * 4.0 * k;
   ro = sqrt(rq);
   rq2 = rq - 2.0;
   rq4 = rq - 4.0;
   r44 = rq4 * 0.25;
   wr = sqrt(rq4);
   l1 = log((ro + wr) * 0.5);
   h2 = w * wr + ro * ww;
   h5 = 2.0 / wq4q * (wq * wq2 * rq4 - rq2 - rq2 + rq * 4.0 / wq);
   a2 = zpia * rq2 / (ro * wr);
   f = a2 * (exp(a1) - 1.0) / (a1 * (exp(a2) - 1.0));

   /* ct = cos(th); */
   m = 1;
   pct = p1 * ct;
   x1 = k * (e1 - pct);
   x2 = k * (e1 + pct);
   x1q = x1 * x1;
   rx1 = 1.0 / x1;
   x3 = x1;
   x3q = x1q;
   rx3 = rx1;
   x2q = x2 * x2;
   rx2 = 1.0 / x2;
   x12 = x1 * x2;
   x = x1 + x2;
   xq = x * x;
   h1 = (x1 - x2) / x;
   xr = x12 / rq;
   w4 = ww * sqrt(w44 * rq4 + xr * 4.0);
   l3 = log(h2 * 0.125 * h2 / x);
   l4 = log((wr * 0.25 * w4 + w44 * 2.0 * r44) / xr + 1.0);
   h3 = rq2 * 0.125 / x * ((rq4 + wq) * (rq4 + wq));
   h4 = wq2 * 8.0 / wq4q + 4.0 + rq * 0.5 * rq2 / x + x * 0.25 / x12 *
     (wq2q + rq2 * rq4 - rq2 * 8.0 / wq4);
   h6 = rq2 * 0.125 / x12
     * (wq2q + rq2 * rq2 - (wq4 + rq) * 6.0
        + x * 16.0 / wq4) + 1.0 + 2.0 / x12 + rq4 / wq4 - 8.0 / wq4q;
   h7 = (rq + wq) * 0.25 / x12 * h1 * h1 - (rx1 - rx2) * (rx1 - rx2) * 0.25
     - rq * 0.5 / xq + 2.0 / (wq4q * xr) * (wq4 + 1.0);
   h8 = rq / wwq * (wq2q * 12.0 / wq4 * x12 - xq * (48.0 / wwq + 8.0));

   next_term:

   r1 = rq4 + (x1 + x1q / rq) * 4.0;
   wr1 = sqrt(r1);
   r2 = rq4 + x1 + x1;
   w2 = sqrt(r44 * x2q + x * 2.0 * xr);
   l = log(ro * 0.25 * rx1 * (r2 + wr * wr1));
   l2 = log(rq * rx1 / x * (r44 * x2 + wr * 0.5 * w2) + 1.0);
   g0 = h7 + h8 / (x1q * x1q) - rq * x / (w44 * x1q * x1) + rq / x1q
     * (4.0 / wq - 1.5 - (8.0 / wq4 - wq2 * x) / wwq) - rq * rx1
     / wq4 + ((wq - x2 * 4.0) / rq - wwq * 0.25 / x12 - rx1 * 4.0
              + (wq - wq2 * 0.5 * rq) / x1q) / r1;
   g1 = (8.0 - p1q - p1q) * rx2 + 4.0 - rx1 * 2.0 * (wq + wq - x2 + 4.)
     + wwq * 0.75 / x12 + rx1 * 4.0 / r2
     * (x2 - x1 * 3.0 + 4.0 + (rq + rq - 6.) * rx1)
       + rx1 / r1 * (wq4 - xr * 4.) * (wq * 0.25 * r2 * rx2 - rq2 - x1 - x1);
   g2 = ro * ((rq + 2.) / xq + 8.0 / x1q);
   g3 = wq2 * 2.0 * rx2 - (rq2 - x2) * rx1 + rq2 * 0.5 * (rq + x1) / x
     + wq4 / r2 - 2.0 + h3 * rx1
     + rx1 / r2 * (rq2 * ((wq2 - x) * x2 - wq2) - x2q - x2q - x2 * 4.)
       + (rq2 * (rq4 * 1.5 - wq * 0.5 * (rq - 5.)) + x1q + x1q - x1 * 6.0
          + wq2 * x2) / (r2 * x);
   g4 = h4 - wq2 * 1.5 * rx2 + rx1 * (x2q - wq - wq - wq * x2 + x2 + wq2q * 0.5)
     - x2 * 0.5 * wq2 / x - rx2 * 0.25 * (wq4 + rq) / x * wq2q + xq * 4.0 * rx1 / wwq
     - (rq2 * 2.0 * x2 + wq4 * x - wq2 * 4.0 + rx1 * 0.5 * (8.0 - rq) * wq2) / r2
     + (rq2 * 0.5 * (rq4 * 3.0 - wq * (rq - 5.)) - x1 * (wq - x1 - x1 + 4.))
       / (x * r2) + h5 * rx1 + wq2 * x / (w44 * x1q) * (x * 12.0 / wwq - rq2 * 0.5)
         * (1.0 - x * rx1 / wq);
   g5 = h6 - (x1 + 2.0 + x1) * rx2;
   s = wr * g0 + g1 * l / wr1 + g2 * l1 + g3 * l2 / w2 + g4 * 2.0
     * ro * l3 / (w * ww * x1) - g5 * l4 / w4;
   s1[m - 1] = s;

   if (m == 2)
     goto finish;

   m = 2;
   x1 = x2;
   x1q = x2q;
   rx1 = rx2;
   x2 = x3;
   x2q = x3q;
   rx2 = rx3;
   goto next_term;

   finish:

   wiq = ar0 * k * (s1[0] + s1[1]) / (pi * ro * w * ww * mcq);
   wiqc = wiq * f;
   /*     WQ: cross section */
   /*     WQC: cross section including Coulomb correction */
   return wiqc;
}

/*}}}*/

struct EE_Type
{
   double een, pen;
};

static double eebrems_diff_lab (double een, double pen, double mu) /*{{{*/
{
   double s, een_cm, pen_cm, cos_theta_cm;
   double gamma, beta, gamma_c, gc2, beta_c;
   double mcq = ELECTRON_REST_ENERGY / KEV;

   /* transform variables to cms */
   gamma = 1.0 + een/mcq;
   gc2 = 0.5*(1.0 + gamma);
   gamma_c = sqrt (gc2);
   beta_c = sqrt (1.0 - 1.0/gc2);

   /* photon (emitted energy=pen along angle=theta in lab frame) */
   pen_cm = pen * gamma_c * (1.0 - beta_c * mu);
   cos_theta_cm = (mu - beta_c)/(1.0 - beta_c * mu);

   /* electron (incident energy=gamma along angle=theta=0 in lab frame) */
   beta = sqrt (1.0 - 1.0/(gamma*gamma));
   een_cm = mcq * (gamma_c * gamma * (1.0 - beta_c * beta) - 1.0);

   s = eebrems_diff_cms (een_cm, pen_cm, cos_theta_cm);
   if (!finite(s))
     {
        fprintf (stderr, "*** eebrems_diff_lab: eebrems_diff_cms => s(%g,%g,%g)=%g\n",
                 een_cm, pen_cm, cos_theta_cm, s);
        exit(1);
     }

   s /= gamma_c * (1.0 - beta_c * mu);

   return s;
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

/* Use this integrand to handle round-off errors in the limit \mu->1, \beta->1.
 * These expressions handle cancellation analytically, computing
 * everything in terms of numerically small values. */
static double eebrems_diff_lab_deltamu (double een, double pen, double delta_mu) /*{{{*/
{
   double s, een_cm, pen_cm, cos_theta_cm;
   double gamma, beta, gamma_c, gc2, beta_c;
   double mcq = ELECTRON_REST_ENERGY / KEV;
   double dbeta, dbeta_c, ombm, mmb, x;

   /* transform variables to cms */
   gamma = 1.0 + een/mcq;
   gc2 = 0.5*(1.0 + gamma);
   gamma_c = sqrt (gc2);
   beta_c = sqrt (1.0 - 1.0/gc2);
   dbeta_c = (gamma_c < GAMMA_THRESH) ? (1.0 - beta_c) : delta_beta(gamma_c);

   /* ombm = 1 - beta_c*mu */
   ombm = delta_mu + dbeta_c - delta_mu * dbeta_c;

   /* mmb = mu - beta_c */
   mmb = dbeta_c - delta_mu;

   /* photon (emitted energy=pen along angle=theta in lab frame) */
   pen_cm = pen * gamma_c * ombm;
   cos_theta_cm = mmb / ombm;

   /* electron (incident energy=gamma along angle=theta=0 in lab frame) */
   beta = sqrt (1.0 - 1.0/(gamma*gamma));
   dbeta = (gamma < GAMMA_THRESH) ? (1.0 - beta) : delta_beta(gamma);

   /* x = gamma_c * gamma * (1.0 - beta_c * beta) */
   x = 1.0 / sqrt ( dbeta * dbeta_c * (2.0 - dbeta_c) * (2.0 - dbeta));
   x *= (dbeta_c + dbeta - dbeta * dbeta_c);
   een_cm = mcq * (x - 1.0);

   s = eebrems_diff_cms (een_cm, pen_cm, cos_theta_cm);
   if (!finite(s))
     {
        fprintf (stderr, "*** eebrems_diff_lab_deltamu: eebrems_diff_cms => s(%g,%g,%g)=%g\n",
                 een_cm, pen_cm, cos_theta_cm, s);
        exit(1);
     }

   s /= gamma_c * ombm;

   return s;
}

/*}}}*/

static unsigned int Num_Evaluations;
static double mu_integrand (double mu, void *p) /*{{{*/
{
   struct EE_Type *info = (struct EE_Type *)p;
   double s;

   /* fprintf (stderr, "%d\r", Num_Evaluations++); */

   /* mu = cos(theta) */
   s = eebrems_diff_lab (info->een, info->pen, mu);

   return s;
}

/*}}}*/

static double delta_mu_integrand (double deltamu, void *p) /*{{{*/
{
   struct EE_Type *info = (struct EE_Type *)p;
   double s;

   /* fprintf (stderr, "%d\r", Num_Evaluations++); */

   /* mu = cos(theta) */
   s = eebrems_diff_lab_deltamu (info->een, info->pen, deltamu);

   return s;
}

/*}}}*/

static void handle_gsl_status (int status) /*{{{*/
{
   if (status == 0)
     return;

   /* FIXME.  It would be nice to avoid complaints about
    * singularities.  I haven't tracked this down yet. */
   if ((status != GSL_EROUND)
       && (status != GSL_EMAXITER)
       && (status != GSL_ESING))
     fprintf (stderr, "*** %s\n", gsl_strerror (status));
}

/*}}}*/

static int angular_integral (double een, double pen, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   struct EE_Type s;
   double epsabs, epsrel, abserr;
   double k, kmax, gamma, beta, mu_min, eps;
   double mcq = ELECTRON_REST_ENERGY / KEV;
   size_t limit;
   int status;

   if (val == NULL)
     return -1;

   *val = 0.0;

   s.een = een;  /* incident electron kinetic energy (keV) */
   s.pen = pen;  /* photon energy (keV) */

   if ((een <= 0.0) || (pen <= 0.0))
     return -1;

   k = pen / mcq;
   gamma = 1.0 + een / mcq;
   beta = sqrt (1.0 - 1.0/(gamma*gamma));

   /* quick return if photon energy exceeds absolute maximum attainable */
   kmax = (gamma - 1.0)/(gamma *(1.0 - beta) + 1.0);
   if (k > kmax)
     return 0;

   /* beaming limits emission to narrow cone in the forward direction
    * mu = cos(theta)
    */
   mu_min = ((gamma + 1.0)*k - (gamma - 1.0)) / (k * gamma * beta);
   if (mu_min > 1.0)
     return 0;
   else if (mu_min < -1.0)
     mu_min = -1.0;

   f.params = &s;
   epsabs = 0.0;
   epsrel = 1.e-13;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   Num_Evaluations = 0;

   /* Break the integral into two pieces to avoid
    * loss of numerical precision.
    *
    * What is the optimal location of the break?
    * The problem is in terms like (\mu - \beta).
    * \beta = sqrt (1 - 1/\gamma^2) so when \gamma=1.e3,
    * 1-\beta = 5.e-7 which is really starting to lose
    * precision as \mu-> 1.  Empirically, breaking
    * the integral at 1-eps where eps=1.e-3 seems
    * to work well, but I don't claim its optimal.
    */
   eps = (een < GAMMA_THRESH) ? 0.0 : 1.e-3;

   if (mu_min < 1.0 - eps)
     {
        /* This piece is straightforward. */
        f.function = &mu_integrand;
        status = gsl_integration_qag (&f, mu_min, 1.0-eps, epsabs, epsrel,
                                      limit, GSL_INTEG_GAUSS61,
                                      work, val, &abserr);
        handle_gsl_status (status);
     }

   if (eps > 0.0)
     {
        double eps_val, tmax;
        
        /* This piece requires special handling.
         * A change of variables and a Taylor series does the trick. 
         */
        
        if (mu_min < (1.0 - eps))
          tmax = eps;
        else tmax = 1.0 - mu_min;

        f.function = &delta_mu_integrand;
        status = gsl_integration_qag (&f, 0.0, tmax, epsabs, epsrel,
                                      limit, GSL_INTEG_GAUSS61,
                                      work, &eps_val, &abserr);
        handle_gsl_status (status);
        *val += eps_val;
     }

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   /* Cross-section is per unit solid angle
    *        d\Omega = d(phi) d(cos(theta))
    * which, with azimuthal symmetry, is
    *        d\Omega = 2*PI* d(cos(theta))
    */

   *val *= 2 * M_PI;

   return 0;
}

/*}}}*/

static int eebrems_lab (double gm1, double e_ph, double *s) /*{{{*/
{
   double mcq = ELECTRON_REST_ENERGY / KEV;
   double een, pen;

   /* een = electron kinetic energy (keV) */
   een = gm1 * mcq;

   /* pen = photon energy (keV) */
   pen = e_ph * mcq;

   return angular_integral (een, pen, s);
}

/*}}}*/

double _ntb_ee_sigma_haug (double gm1, double e_ph) /*{{{*/
{
   double s;

   if (-1 == eebrems_lab (gm1, e_ph, &s))
     {
        fprintf (stderr, "*** _ntb_ee_sigma_haug:  angular integration failed\n");
        return -1.0;
     }

   s *= BARN * ELECTRON_REST_ENERGY / KEV;

   return s;
}

/*}}}*/

static double ntb_integrand (double pc, void *pt) /*{{{*/
{
   Brems_Type *b = (Brems_Type *)pt;
   Particle_Type *elec = b->electrons;
   double ne, pcomc2, gamma, gm1, beta, v, s;

   s = 0.0;

   if (pc <= 0.0)
     return s;

   pcomc2 = pc / ELECTRON_REST_ENERGY;
   gamma = sqrt (1.0 + pcomc2 * pcomc2);
   gm1 = gamma - 1.0;

   if (b->ee_weight != 0.0)
     {
        double s_ee;
        if (b->interpolate)
          {
             if (-1 == ntb_interp_sigma (b, gm1, b->photon_energy, &s_ee))
               s_ee = 0.0;
          }
        else
          {
             s_ee = _ntb_ee_sigma_haug (gm1, b->photon_energy);
          }

        s += b->ee_weight * s_ee;
     }

   if (b->ep_weight != 0.0)
     {
        double s_ep = _ntb_ei_sigma (gm1, b->photon_energy);
        s += b->ep_weight * s_ep;
     }

   (void)(*elec->spectrum) (elec, pc, &ne);

   beta = sqrt ((gamma + 1.0)*(gamma-1.0))/gamma;
   v = beta * GSL_CONST_CGSM_SPEED_OF_LIGHT;

   return v * ne * s;
}

/*}}}*/

static int integral_over_electrons (Brems_Type *b, double *val) /*{{{*/
{
   gsl_error_handler_t *gsl_error_handler;
   gsl_integration_workspace *work;
   gsl_function f;
   double epsabs, epsrel, abserr;
   double pc_min, pc_max;
   size_t limit;
   int status;

   pc_min = b->photon_energy * ELECTRON_REST_ENERGY;
   pc_max = GAMMA_MAX_DEFAULT * ELECTRON_REST_ENERGY;

   f.function = &ntb_integrand;
   f.params = b;
   epsabs = 0.0;
   epsrel = 1.e-9;
   limit = MAX_QAG_SUBINTERVALS;

   if (NULL == (work = gsl_integration_workspace_alloc (limit)))
     return -1;

   gsl_error_handler = gsl_set_error_handler_off ();

   status = gsl_integration_qag (&f, pc_min, pc_max, epsabs, epsrel,
                                 limit, GSL_INTEG_GAUSS31, work,
                                 val, &abserr);

   gsl_set_error_handler (gsl_error_handler);
   gsl_integration_workspace_free (work);

   if (status)
     {
        if (status != GSL_EROUND)
          fprintf (stderr, "*** %s\n", gsl_strerror (status));
     }

   return 0;
}

/*}}}*/

int ntb_brems (void *vb, double photon_energy, double *emissivity)
{
   Brems_Type *b = (Brems_Type *)vb;

   /* photon energy eV -> units of electron mc^2 */
   b->photon_energy = photon_energy
     * GSL_CONST_CGSM_ELECTRON_VOLT / ELECTRON_REST_ENERGY;

   if (-1 == integral_over_electrons (b, emissivity))
     return -1;
   
   *emissivity /= ELECTRON_REST_ENERGY;
   
   return 0;
}

#undef TEST
#ifdef TEST
int main (void) /*{{{*/
{
   double emin, ee, eph;
   unsigned int n;

   eph = 1.1;

   if (-1 == brems_emin (eph, &emin))
     {
        fprintf (stderr, "couldnt find emin for eph = %g\n", eph);
        return 1;
     }

   ee = emin;
   n = 1024;
   while (n-- > 0)
     {
        fprintf (stdout, "%12.4e  %12.4e\n",
                 log10(ee),
                 _ntb_ee_sigma_haug (ee,  eph));
        ee *= 1.02;
     }

   return 0;
}

/*}}}*/
#endif
