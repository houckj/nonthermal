/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "photons.h"

/*{{{ physical constants */

#define COMPTON_WAVELENGTH \
   (GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR \
       / (GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT))

#define BBODY_COEF \
   (1.0 / ((M_PI * M_PI) \
              * (COMPTON_WAVELENGTH \
                    * COMPTON_WAVELENGTH * COMPTON_WAVELENGTH)))

#define BBODY_TSCALE \
   (ELECTRON_REST_ENERGY / GSL_CONST_CGSM_BOLTZMANN)

/*}}}*/

static int blackbody_photons (double e, double t, double *y) /*{{{*/
{
   double max_exp = 500.0;
   double x, f;

   *y = 0.0;

   if (t <= 0.0)
     return 0;

   /* e: in units of electron rest-energy */
   /* t: Kelvin -> dimensionless */
   t /= BBODY_TSCALE;

   x = e / t;

   if (x > max_exp)
     {
        f = 0.0;
     }   
   else if (x > 1.e-3)
     {
        f = e * e /(exp (x) - 1.0);
     }
   else
     {
        double den =
            1.0 + (x/2)
            *(1.0 + (x/3)
              *(1.0 + (x/4)
                *(1.0 + (x/5)
                  *(1.0 + (x/6)
                    *(1.0 + (x/7)
                      *(1.0 + (x/8)))))));

        f = e * t / den;
     }

   *y = f * BBODY_COEF;

   return 0;
}

/*}}}*/

static double Photon_Temperature = CBR_TEMPERATURE;

double incident_photon_max_energy (void)
{
   double peak_energy = (Photon_Temperature 
                         * GSL_CONST_CGSM_BOLTZMANN / GSL_CONST_CGSM_ELECTRON_VOLT);
   return 50.0 * peak_energy;
}

double incident_photon_min_energy (void)
{
   double peak_energy = (Photon_Temperature 
                         * GSL_CONST_CGSM_BOLTZMANN / GSL_CONST_CGSM_ELECTRON_VOLT);
   /* min energy that can scatter to an energy of interest */
   return 1.e-11 * peak_energy;
}

int incident_photon_spectrum (double e, double *num)
{
   double t = Photon_Temperature;
   return blackbody_photons (e, t, num);
}

void set_incident_photon_kelvin_temperature (double t)
{
   Photon_Temperature = t;
}
