#ifndef PION_H
#define PION_H

typedef struct
{
   Particle_Type *protons;
   double energy;
}
Pion_Type;

#define NULL_PION_TYPE {NULL,0.0}
   
extern int pizero_decay (Pion_Type *p, double photon_energy, double *emissivity);

#endif
