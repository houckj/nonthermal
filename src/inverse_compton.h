#ifndef IC_INVERSE_COMPTON_H
#define IC_INVERSE_COMPTON_H

typedef struct
{
   Particle_Type *electrons;
   int (*incident_photons) (double, double *);
   double electron_gamma;
   double energy_final_photon;            /* eV */
   double incident_photon_max_energy;     /* eV */
   double incident_photon_min_energy;     /* eV */
   int interpolate;
   int complain_on_extrapolate;
   void *client_data;
}
Inverse_Compton_Type;

#define NULL_INVERSE_COMPTON_TYPE  {NULL,NULL,0.0,0.0,0.0,0.0,0,0,NULL}

extern double ic_knlimit_constant (double p);
extern int ic_calc_inverse_compton
  (void *vic, double energy_final_photon, double *emissivity);

extern int
  ic_integral_over_incident_photons (Inverse_Compton_Type *ic, double *val);

#endif
