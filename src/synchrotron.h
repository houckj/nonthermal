
#ifndef SYN_SYNCHRO_H
#define SYN_SYNCHRO_H

typedef struct 
{
   Particle_Type *electrons;
   double B_tot;
   double photon_energy;
   int interpolate;
   void *client_data;
}
Synchrotron_Type;

#define NULL_SYNCHROTRON_TYPE  {NULL,0.0,0.0,0,NULL}

extern int syn_calc_synchrotron (void *vs, double photon_energy, double *emissivity);
extern int syn_angular_integral (double x, double *y);

#endif
