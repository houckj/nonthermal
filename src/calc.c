/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "photons.h"
#include "synchrotron.h"
#include "syn_table.h"
#include "inverse_compton.h"
#include "ic_table.h"
#include "ntbrems.h"
#include "pion.h"

#undef BENCHMARK

int main (void)
{
   Particle_Type electron = NULL_PARTICLE_TYPE;
   Particle_Type proton = NULL_PARTICLE_TYPE;
   Synchrotron_Type s = NULL_SYNCHROTRON_TYPE;
   Inverse_Compton_Type ic = NULL_INVERSE_COMPTON_TYPE;
   Pion_Type p = NULL_PION_TYPE;
#if 0   
   Brems_Type b = NULL_BREMS_TYPE;
#endif   
   double lg_e;
   
   if (-1 == SLang_init_all())
     return 1;

   /* Electron spectrum */
   (void) init_particle_spectrum (&electron);
   electron.index = -2.0;
   electron.cutoff_energy = 10.0;  /* TeV */
   electron.mass = GSL_CONST_CGSM_MASS_ELECTRON;

   /* Proton spectrum */
   (void) init_particle_spectrum (&proton);
   proton.index = -2.0;
   proton.cutoff_energy = 10.0;  /* TeV */
   proton.mass = GSL_CONST_CGSM_MASS_PROTON;
   
   /* Pion setup */
   p.protons = &proton;
   p.density = 1.0;
   
   /* Synchrotron setup */
   s.electrons = &electron;
   s.B_tot = 1.e-6;
   s.interpolate = 1;
   s.client_data = syn_load_table (SYN_TABLE_FILE_DEFAULT);
   if (s.client_data == NULL)
     goto return_error;

   /* Inverse Compton setup */
   ic.electrons = &electron;
   ic.incident_photons = &incident_photon_spectrum;
   ic.incident_photon_max_energy = incident_photon_max_energy();
   ic.interpolate = 1;
   ic.complain_on_extrapolate = 0;

   ic.client_data = ic_init_client_data (IC_TABLE_FILE_DEFAULT);
   if (ic.client_data == NULL)
     goto return_error;
   
#if 0
   /* Nonthermal Bremsstrahlung setup */
   b.electrons = &electron;

   for (lg_e = -9.0; lg_e < 15.0; lg_e += 0.25)
     {
        double e, f;

        /* photon energy in eV */
        e = pow (10.0, lg_e);
        
        (void) nt_ee_brems (&b, e, &f);

        fprintf (stderr, "%12.4e  %12.4e\n", e, f);
     }
#else

#ifdef BENCHMARK
   {int i;
   for (i = 1; i < 10; i++){
#endif

   ic.interpolate = 1;
   for (lg_e = -9.0; lg_e < 15.0; lg_e += 0.25)
     {
        double e, synch, invc, pi_gamma, e_MeV;

        /* photon energy in eV */
        e = pow (10.0, lg_e);

        if (e < 1.e5)
          (void) syn_calc_synchrotron (&s, e, &synch);
        else
          synch = 0.0;

        if (e > 1.e3)
          (void) ic_calc_inverse_compton (&ic, e, &invc);
        else
          invc = 0.0;

        if (e > 1.e3)
          (void) pi_calc_pion_decay (&p, e, &pi_gamma);
        else
          pi_gamma = 0.0;
        
        /* E^2 dn/dtdE = MeV cm^(-3) s^(-1) */
        e_MeV = e /1.e6;
#ifndef BENCHMARK
        fprintf (stdout, "%12.8e  %12.8e  %12.8e   %12.8e\n",
                 e_MeV,
                 e_MeV * (1.e3 * e_MeV) * synch,
                 e_MeV * (1.e3 * e_MeV) * invc,
                 e_MeV * (1.e3 * e_MeV) * pi_gamma
                 );
        fflush (stdout);
#endif
     }

#ifdef BENCHMARK
   }
   }
#endif

#endif

   return_error:
   if (ic.interpolate) ic_free_client_data (ic.client_data);
   if (s.interpolate) syn_free_table (s.client_data);

   return 0;
}

/*}}}*/

