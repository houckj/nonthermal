#ifndef NTB_TABLE_H
#define NTB_TABLE_H

#include "ntbrems.h"

#define NTB_TABLE_FILE_DEFAULT "ntb_table.dat"
#define NTB_TABLE_ENV  "NTB_TABLE_DIR"

extern void ntb_free_client_data (void *);
extern void *ntb_init_client_data (const char *file);
extern int ntb_make_table (const char *file, int process);
/* extern void ntb_set_process_weights (void *v, double *pw, unsigned int n); */

extern int ntb_interp_sigma (Brems_Type *ntb, double gamma, double ephoton, double *value);

#endif
