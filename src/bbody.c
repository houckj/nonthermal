/* -*- mode: C; mode: fold -*- */

/*  author:  John Houck <houck@space.mit.edu>
 * written:  Aug 2002
 */

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
   double min_exp = 1.e-13;
   double max_exp = 500.0;
   double x;

   *y = 0.0;

   if (t <= 0.0)
     return 0.0;

   /* e: in units of electron rest-energy */
   /* t: Kelvin -> dimensionless */
   t /= BBODY_TSCALE;

   x = e / t;

   if (x < min_exp)
     *y = e * t;
   else if (x > max_exp)
     *y = 0.0;
   else
     *y = e * e /(exp (x) - 1.0);

   *y *= BBODY_COEF;

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

int incident_photon_spectrum (double e, double *num)
{
   double t = Photon_Temperature;
   return blackbody_photons (e, t, num);
}

void set_incident_photon_kelvin_temperature (double t)
{
   Photon_Temperature = t;
}
