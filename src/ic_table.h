/*
  Copyright (C) 2002, 2003, 2004, 2005, 2006 John C. Houck 

  This file is part of the nonthermal module

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
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
