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
#ifndef SYN_TABLE_H
#define SYN_TABLE_H

#include "synchrotron.h"

extern int syn_interp_angular_integral (void *p, double x, double *y);

extern void *syn_alloc_table (unsigned int size);
extern int syn_push_table (void *p);
extern void syn_free_table (void *p);
extern void *syn_load_table (char *file);

#endif
