#include "cfortran.h"
PROTOCCALLSFFUN2(DOUBLE,EECMS,eecms,DOUBLE,DOUBLE)
#define EECMS(e,p)  CCALLSFFUN2(EECMS,eecms,DOUBLE,DOUBLE,e,p)
double _ntb_ee_sigma_haug (double electron_kinetic_energy, double photon_energy) /*{{{*/
{
   double ee_kev, pe_kev, sigma;

   /* electron_kinetic_energy = gamma-1 */
   ee_kev = electron_kinetic_energy * ELECTRON_REST_ENERGY / KEV;
   pe_kev = photon_energy * ELECTRON_REST_ENERGY / KEV;

   /* returns differential cross-section in CM system in barns/keV */
   sigma = EECMS(ee_kev, pe_kev);

   if (!finite(sigma))
     sigma = 0.0;

#ifndef TEST
   sigma *= BARN * ELECTRON_REST_ENERGY / KEV;
#endif

   return sigma;
}

/*}}}*/
