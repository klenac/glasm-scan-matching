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

#ifndef GCPH
#define GCPH

#include "simulation.h"
#include <math.h>
#include "utils.h"

namespace GCP
{

void setGCPSearchArea(double,double,double,double,double,double);
void setGCPRefScan(const scan&);
void setGCPNewScan(const scan&);
void setGCPCorrGridParameters(int,int,double,double,double,double,double);
void setGCPGeneticParameters(int,int,int,int,int,int,double,double);

void initializeGCPCorrGrid(void);
void deleteGCPCorrGrid(void);

void initializeGCPCorrGrid_ORIGINAL(void);
void initializeGCPCorrGrid_SUGGESTED(void);

position gcp();

} // namespace

#endif
