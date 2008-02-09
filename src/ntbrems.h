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
#ifndef NTBREMS_H
#define NTBREMS_H

typedef struct Brems_Type Brems_Type;
struct Brems_Type
{
   Particle_Type *electrons;  /* relativistic electrons */
   Particle_Type *e_target;   /* non-relativistic */
   Particle_Type *i_target;   /* non-relativistic */
   double photon_energy;
   double ee_weight;
   double ep_weight;
   double pc_t, pc_r;         /* momenta */
   int interpolate;
   void *client_data;
};

#define NULL_BREMS_TYPE  {NULL,NULL,NULL,0.0,0.0,0.0,0.0,0.0,0,NULL}

enum
{
   NTB_ee = 0,
   NTB_ep = 1   
};

/* electron_kinetic_energy = (\gamma-1)
 * photon_energy in units of electron-rest-energy
 * differential cross-section in CM system in barns/keV 
 */

extern double _ntb_ee_sigma_haug (double electron_kinetic_energy, double photon_energy);
/* _ntb_ee_sigma_haug:  (preferred method for cross-sections)
 *    sigma computed in the CM frame, then a Lorentz transform
 *    is applied to get the lab-frame cross-section, which
 *    which is then integrated over angles.
 */

extern double _ntb_ee_sigma_haug_lab (double electron_kinetic_energy, double photon_energy);
/* _ntb_ee_sigma_haug_lab:  (method used for checking)
 *    sigma is computed in the lab-frame, then integrated over
 *    angles.  The lab-frame code is written in such a way to
 *    make it prone to lose significant figures so it doesnt
 *    work very well in the ultra-relativistic regime.  
 *    The CM frame code suffers fromt he same problem, but
 *    the gamma values used in the CM frame are the sqrt()
 *    of the lab-frame gamma, so the precision loss is much less
 *    important for electron energies of interest.
 */ 


extern double _ntb_ei_sigma (double electron_kinetic_energy, double photon_energy);
extern int ntb_brems_stationary (void *vs, double photon_energy, double *emissivity);
extern int ntb_brems_non_stationary (void *vs, double photon_energy, double *emissivity);

#endif
