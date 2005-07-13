#ifndef PIZERO_H
#define PIZERO_H

typedef struct
{
   Particle_Type *protons;
   double energy;
}
Pizero_Type;

#define NULL_PIZERO_TYPE {NULL,0.0}

extern double Pizero_Approx_Min_Energy;
extern int pizero_decay (void *p, double photon_energy, double *emissivity);

#endif
