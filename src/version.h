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
#define MODULE_MAJOR_VERSION	1
#define MODULE_MINOR_VERSION	0
#define MODULE_PATCH_LEVEL	7
#define MKSTR1(x) #x
#define MKSTR(x) MKSTR1(x)
static char *Module_Version_String = MKSTR(MODULE_MAJOR_VERSION) "." \
   MKSTR(MODULE_MINOR_VERSION) "." MKSTR(MODULE_PATCH_LEVEL);

#define MODULE_VERSION_NUMBER	\
   (MODULE_MAJOR_VERSION*10000+MODULE_MINOR_VERSION*100+MODULE_PATCH_LEVEL)
