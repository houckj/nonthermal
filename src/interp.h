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
#ifndef INTERP_H
#define INTERP_H

extern int Bspline_Type_Id;

typedef struct Bspline_Info_Type Bspline_Info_Type;
/* This type is for the S-Lang model class */
typedef struct
{
   Bspline_Info_Type *info;
}
Bspline_Type;

extern void bspline_free (Bspline_Info_Type *q);
extern Bspline_Info_Type *bspline_open
  (int kx, int ky, double *x, int nx, double *y, int ny, double **f);
extern Bspline_Info_Type *copy_bspline_info (void);
extern double bspline_eval (Bspline_Info_Type *q, double x, double y);

extern void bspline_init_info_intrin (Bspline_Type *qt);
extern void destroy_bspline_type (SLtype type, VOID_STAR f);
extern double bspline_eval_intrin (Bspline_Type *qt, double *x, double *y);
extern void bspline_open_intrin (int *kx, int *ky);

#endif
