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
 * differential cross-section in CM system in barns/keV */
extern double _ntb_ee_sigma_haug (double electron_kinetic_energy, double photon_energy);
extern double _ntb_ei_sigma (double electron_kinetic_energy, double photon_energy);

extern int ntb_brems (void *vs, double photon_energy, double *emissivity);

#endif
