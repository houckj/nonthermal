#ifndef SYN_TABLE_H
#define SYN_TABLE_H

#include "synchrotron.h"

extern int syn_angular_integral (double x, double *y);
extern int syn_interp_angular_integral (void *p, double x, double *y);

extern void *syn_alloc_table (unsigned int size);
extern int syn_create_table (void *p);
extern int syn_push_table (void *p);
extern void syn_free_table (void *p);
extern void *syn_load_table (char *file);

#endif
