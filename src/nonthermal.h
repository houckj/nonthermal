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
#ifndef USER_NONTHERMAL_H
#define USER_NONTHERMAL_H

typedef struct Particle_Type Particle_Type;
struct Particle_Type
{
   Particle_Type *next;
   char *method;
   int (*spectrum)(Particle_Type *, double, double *);
   int (*spectrum_nl)(Particle_Type *, double, double *);
   double (*momentum_min)(Particle_Type *);
   double (*momentum_max)(Particle_Type *);
   double *params;
   unsigned int num_params;
   double mass;
};
#define NULL_PARTICLE_TYPE {NULL,NULL,NULL,NULL,NULL,NULL,0,0.0}
#define PARTICLE_METHOD(n,num,m,mn,mx,mnl) {NULL,(n),m,mnl,mn,mx,NULL,num,0.0}

#define NONTHERMAL_PDF_MODULE(name,p,o) \
   extern int Pdf_##name##_init(Particle_Type *p, char *o); \
          int Pdf_##name##_init(Particle_Type *p, char *o)

#endif
