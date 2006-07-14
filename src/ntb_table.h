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
#ifndef NTB_TABLE_H
#define NTB_TABLE_H

#include "ntbrems.h"

#define NTB_TABLE_ENV  "NTB_TABLE_DIR"

extern void ntb_free_client_data (void *);
extern void *ntb_init_client_data (const char *file);
/* extern void ntb_set_process_weights (void *v, double *pw, unsigned int n); */

extern int ntb_interp_sigma (Brems_Type *ntb, double gm1, double ephoton, double *value);

#endif
