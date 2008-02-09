/*
  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008 John C. Houck 

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
