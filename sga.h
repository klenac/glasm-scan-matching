/*  
 *  Copyright (C) 2014  Kristijan Lenac
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "simulation.h"

//---------------------------------------------------------------------------
#ifndef sgaH
#define sgaH
//---------------------------------------------------------------------------


struct individual
{   unsigned chrom;		// 6/3/2011 carefull! max 32 bit per chromosome in this version
	double   fitness;
	int      xsite;
	int      parent[2];
};

struct bestever
{	unsigned chrom;
	double   fitness;
	int      generation;
};


void sga_parameters(int,int,int,int,int,int,double,double);
void set_sga_objfun(void (*objfun)(struct individual *critter));
#ifdef DRAW_PNG
	void set_sga_drawfun(void (*drawfun)(int run, int gen, int popsize, struct individual *critter, struct bestever *bestfit));
#endif
bestever sga(void);






#endif
