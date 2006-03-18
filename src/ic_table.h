#ifndef IC_TABLE_H
#define IC_TABLE_H

#include "inverse_compton.h"

#define IC_TABLE_ENV  "IC_TABLE_DIR"

extern void ic_free_client_data (void *);
extern void *ic_init_client_data (const char *file);
extern int ic_interp_photon_integral (Inverse_Compton_Type *ic, double *value);
extern void ic_set_dilution_factors (void *v, double *df, unsigned int n);

extern double ic_omega0 (Inverse_Compton_Type *ic);

#endif
