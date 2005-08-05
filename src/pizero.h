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

extern int Pizero_Use_Dermer_Xsec;
extern double Pizero_Approx_Min_Energy;
extern int pizero_decay (void *p, double photon_energy, double *emissivity);
extern int pizero_distribution (Pizero_Type *p, double *val);
extern double pizero_diff_xsec (double T_p, double T_pi);

#endif
