#ifndef PIZERO_H
#define PIZERO_H

typedef struct
{
   Particle_Type *protons;
   double energy;
   int interpolate;
   void *client_data;
}
Pizero_Type;

#define NULL_PIZERO_TYPE {NULL,0.0, 0,NULL}

extern double Pizero_Approx_Min_Energy;
extern int pizero_decay (void *p, double photon_energy, double *emissivity);

#endif
