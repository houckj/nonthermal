#ifndef PHOTONS_H
#define PHOTONS_H

#include "_nonthermal.h"

extern double incident_photon_max_energy (void);
extern double incident_photon_min_energy (void);
extern int incident_photon_spectrum (double e, double *num);
extern void set_incident_photon_kelvin_temperature (double t);

#endif

