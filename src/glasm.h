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

#ifndef GLASMH
#define GLASMH


#include <math.h>
#include "utils.h"

namespace GLASM
{

void setSearchArea(double,double,double,double,double,double);
void setCenterOfSearchArea(position);
void setRefScan(const scan&);
void setNewScan(const scan&);
void setLookupTableParameters(int,int,double,double,double,double,double);
void setGeneticParameters(int,int,int,int,int,int,double,double);
void setDebugParameters(bool udraw_lookup,
		bool udraw_individual,
		bool udraw_generations,
		unsigned uverbose_level=1,
		double umap_size_x=15.0,
		double umap_size_y=15.0,
		std::string umap_filename="./cfg/bitmaps/empty_1000x1000.png",
		std::string ulookup_bitmap="map.png",
		std::string uvalid_bitmap="map.png");

void initialize_binary_lookup_table();
void initialize_gradient_lookup_table();
void initialize_lookup_table_from_bitmap(std::string,double,double,double);
void initialize_valid_table_from_bitmap(std::string);
void initialize_radial_gradient(double,double,double,double,double,double);
void initialize_invGray_table();
void delete_lookup_table();

position glasm(double* outfitness);

}

#endif
