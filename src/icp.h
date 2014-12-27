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

#ifndef ICPH
#define ICPH

//---------------------------------------------------------------------------

#ifndef utilsH
    #include "utils.h"
#endif

namespace ICP
{


void setCorrSearchMethod(char i);
void setParameters(double ucorr_angle_treshold,
        double usensor_range,
        double ustepangle,
        double umaxfloat,
        double ucorr_distance_treshold,
        unsigned umax_iteration,
        double umin_ddisp,
        double umin_rdisp);
void setDebugParameters(bool udraw_iterations,
		unsigned uverbose_level,
		double umap_size_x,
		double umap_size_y,
		std::string map_filename);



position icp(scan&, scan&, unsigned&);

std::list<corrPoints> &  findCorrPoints(scan& S1, scan& S2, double corr_distance_treshold);

} // namespace

#endif
