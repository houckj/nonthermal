#ifndef PIZERO_TABLE_H
#define PIZERO_TABLE_H

#include "pizero.h"

#define PIZERO_TABLE_SIZE 2048

extern void *pizero_alloc_table (unsigned int size);
extern void pizero_free_table (void *p);
extern int pizero_spline_table (void *p, double *x, double *y, unsigned int n);
extern int pizero_interp_pizero_integral (void *p, double x, double *y);

#endif
