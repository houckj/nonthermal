#ifndef PION_H
#define PION_H

typedef struct
{
   Particle_Type *protons;
   double proton_kinetic_energy;
   double photon_energy;
   double density;
}
Pion_Type;

#define NULL_PION_TYPE {NULL,0.0,0.0,0.0}
   
extern int pi_calc_pion_decay (Pion_Type *p, double photon_energy, double *val);

#endif
