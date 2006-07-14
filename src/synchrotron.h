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
