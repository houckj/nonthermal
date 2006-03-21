#ifndef NTBREMS_H
#define NTBREMS_H

typedef struct
{
   Particle_Type *electrons;
   double photon_energy;
   double ee_weight;
   double ep_weight;
   void *client_data;
   int interpolate;
}
Brems_Type;

#define NULL_BREMS_TYPE  {NULL,0.0,0.0,0.0,NULL,0}

enum
{
   NTB_ee = 0,
   NTB_ep = 1   
};

/* electron_kinetic_energy = (\gamma-1)
 * photon_energy in units of electron-rest-energy
 * differential cross-section in CM system in barns/keV 
 */

extern double _ntb_ee_sigma_haug (double electron_kinetic_energy, double photon_energy);
/* _ntb_ee_sigma_haug:  (preferred method for cross-sections)
 *    sigma computed in the CM frame, then a Lorentz transform
 *    is applied to get the lab-frame cross-section, which
 *    which is then integrated over angles.
 */

extern double _ntb_ee_sigma_haug_lab (double electron_kinetic_energy, double photon_energy);
/* _ntb_ee_sigma_haug_lab:  (method used for checking)
 *    sigma is computed in the lab-frame, then integrated over
 *    angles.  The lab-frame code is written in such a way to
 *    make it prone to lose significant figures so it doesnt
 *    work very well in the ultra-relativistic regime.  
 *    The CM frame code suffers fromt he same problem, but
 *    the gamma values used in the CM frame are the sqrt()
 *    of the lab-frame gamma, so the precision loss is much less
 *    important for electron energies of interest.
 */ 


extern double _ntb_ei_sigma (double electron_kinetic_energy, double photon_energy);

extern int ntb_brems (void *vs, double photon_energy, double *emissivity);

#endif
