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

#include <string.h>

#include "pizero.h"
#include "cparamlib/cparamlib.h"

static int init_pizero (double *par, unsigned int npar, Pizero_Type *p, Particle_Type *proton) /*{{{*/
{
   (void) par; (void) npar;

   if (-1 == init_pdf (proton, PROTON))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   p->protons = proton;

   return 0;
}

/*}}}*/

static int pizero_init_client_data (void) /*{{{*/
{
   return 0;
}

/*}}}*/

static int binned_pizero (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   Pizero_Type p = NULL_PIZERO_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   PARAMSET params;
   double dlg_pc, lg_pc_min, lg_pc_max, pc_lo, np_lo;
   double *bin_lo, *bin_hi;
   double norm = par[0];
   int i, j, num_pc_bins;
   int n_notice, *notice_list;

   (void) npar;

   for (i = 0; i < g->n_notice; i++)
     {
        val[i] = 0.0;
     }

   if (norm == 0.0)
     return 0;

   if (-1 == init_pizero (par, npar, &p, &proton))
     return -1;

   bin_lo = g->bin_lo;
   bin_hi = g->bin_hi;
   n_notice = g->n_notice;
   notice_list = g->notice_list;

   pc_lo = PROTON_REST_ENERGY;  /* FIXME */
   (void)(*proton.spectrum)(&proton, pc_lo, &np_lo);

   lg_pc_min = log(pc_lo);
   lg_pc_max = log(proton.momentum_max (&proton));
   dlg_pc = 0.025;
   num_pc_bins = 1 + (lg_pc_max - lg_pc_min)/dlg_pc;
   
   memset (&params, 0, sizeof(PARAMSET));

   for (j = 1; j < num_pc_bins; j++)
     {
        double pc, v, pc_hi, np_hi, np_mid, num_p, x, Tp, t;

        pc_hi = exp(lg_pc_min + dlg_pc * j);
        (void)(*proton.spectrum)(&proton, pc_hi, &np_hi);

        pc = 0.5 * (pc_lo + pc_hi);
        (void)(*proton.spectrum)(&proton, pc, &np_mid);

        /* Simpson's rule */
        num_p = (pc_hi - pc_lo) * (np_lo + 4*np_mid + np_hi) / 6.0 ;

        x = pc / PROTON_REST_ENERGY;
        Tp = PROTON_REST_ENERGY * (sqrt (1.0 + x*x) - 1.0);

        v = (GSL_CONST_CGSM_SPEED_OF_LIGHT
             * pc / (Tp + PROTON_REST_ENERGY));

        t = Tp / GEV;

        gamma_param (t, &params);

        for (i=0; i < n_notice; i++)
          {
             double s_nd, s_diff, s_delta, s_res, s;
             double e, el, eh;
             int n;

             n = notice_list[i];

             /* Angstrom -> GeV */
             el = 1.e-9 * EV_ANGSTROM / bin_hi[n];
             eh = 1.e-9 * EV_ANGSTROM / bin_lo[n];
             e = 0.5 * (el + eh);

             s_nd = sigma_incl_nd (ID_GAMMA, e, t, &params);
             s_diff = sigma_incl_diff (ID_GAMMA, e, t, &params);
             s_delta = sigma_incl_delta (ID_GAMMA, e, t, &params);
             s_res = sigma_incl_res (ID_GAMMA, e, t, &params);

             s = (s_nd + s_diff + s_delta + s_res) * (eh - el) / e;
             val[i] += num_p * v * s;
          }

        pc_lo = pc_hi;
        np_lo = np_hi;
     }

   free_pdf (&proton);

   /* norm = (V/4*pi*D^2) * n0 (cm^-2 GeV^-1)
    * s(E)dE = s(y)dy = photons cm^-2 s^-1
    */
   norm *= MILLIBARN / GEV;
   for (i = 0; i < n_notice; i++)
     {
        val[i] *= norm;
     }

   return 0;
}

/*}}}*/

static int unbinned_pizero (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/ /*{{{*/
{
   Isis_Hist_t gg;
   double *x, *hi, *lo;
   double yl, yh;
   int i, status, n;

   /* g->x[i] in eV
    * val in photons /s /cm^2 /GeV
    */

   n = g->npts;
   x = g->x;

   if (-1 == Isis_Hist_allocate (n, &gg))
     return -1;

   hi = gg.bin_hi;
   lo = gg.bin_lo;

   /* make interior bins first */
   for (i = 1; i < n-1; i++)
     {
        int j;
        j = n-i-1;
        yl = EV_ANGSTROM * 0.5 * (1.0 /x[j+1] + 1.0 /x[j]);
        yh = EV_ANGSTROM * 0.5 * (1.0 /x[j] + 1.0 /x[j-1]);
        lo[i] = yl;
        hi[i] = yh;
     }
   /* first bin */
   yl = EV_ANGSTROM /x[n-1];
   lo[0] = yl - (hi[0] - yl);
   hi[0] = lo[1];
   /* last bin */
   yh = EV_ANGSTROM / x[0];
   lo[n-1] = hi[n-2];
   hi[n-1] = yh + (yh - lo[n-1]);

   /* update notice list */
   gg.n_notice = n;
   for (i = 0; i < n; i++)
     {
        gg.notice[i] = 1;
        gg.notice_list[i] = i;
     }

   status = binned_pizero (val, &gg, par, npar);

   for (i = 0; i < n; i++)
     {
        double el, eh;
        int k = gg.notice_list[i];
        el = 1.e-9 * EV_ANGSTROM / hi[k];
        eh = 1.e-9 * EV_ANGSTROM / lo[k];
        val[i] /= (eh - el);
     }

   /* reverse order for output */
   for (i = 0; i < n/2; i++)
     {
        double tmp = val[i];
        val[i] = val[n-i-1];
        val[n-i-1] = tmp;
     }

   Isis_Hist_free (&gg);

   return status;
}

/*}}}*/

