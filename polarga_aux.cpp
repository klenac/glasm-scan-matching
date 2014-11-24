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

#include "polarga_aux.h"
// polarga_aux.cpp functions help interfacing to sm_cli simulator


#include <iostream>
#include <fstream>
//#include <string>
//#include <stdlib.h>
//#include <stdio.h>
//#include <locale.h>
//#include <math.h>
//#include <sys/time.h>
#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "simulation.h"

#ifdef DRAW_PNG
		#include "mydraw.h"
#endif


namespace POLARGA
{

double          outfitness;

void init(const char* filename)
{
	bool draw_individual;
	bool draw_generations;
	unsigned verbose_level;
	std::string map_filename;
	double map_size_x;
	double map_size_y;

	// PolarGA - Martinez
	int             nbitx;
	int             nbity;
	int             nbitrot;
	int             pcross;
	double          pmutation;
	int             popsize;
	int             maxruns;
	int             maxgen;

	//int             RisulArea;
	double          MinX;
	double          MaxX;
	double          MinY;
	double          MaxY;
	double          MinRot;
	double          MaxRot;
	double          SearchAreaDist;
	double          SearchAreaRot;

	int             NREADINGS;
	double          CORR_DISTANCE_TRESHOLD;
	double          FOV;
	bool			RELAXED_CORR_SEARCH;


	try {

	std::ifstream config_stream(filename);
	namespace po = boost::program_options;
	po::options_description my_options("Options");

	my_options.add_options()
		("draw_individual", po::value<bool>(&draw_individual)->default_value("false"), "set true to draw_individual")
		("draw_generations", po::value<bool>(&draw_generations)->default_value("false"), "set true to draw_generations ")
		("verbose_level", po::value(&verbose_level)->default_value(0), "verbose_level")
		("nbitx", po::value(&nbitx)->default_value(8), "nbitx")
		("nbity", po::value(&nbity)->default_value(8), "nbity")
		("nbitrot", po::value(&nbitrot)->default_value(8), "nbitrot")
		("pcross", po::value(&pcross)->default_value(8), "pcross")
		("pmutation", po::value(&pmutation)->default_value(8), "pmutation")
		("popsize", po::value(&popsize)->default_value(8), "popsize")
		("maxruns", po::value(&maxruns)->default_value(8), "maxruns")
		("maxgen", po::value(&maxgen)->default_value(8), "maxgen")
		// search area (specify either exact range or by distance and rotation. the former takes precedence if both specified)
		("MinX", po::value(&MinX)->default_value(0), "SearchArea MinX")
		("MaxX", po::value(&MaxX)->default_value(1), "SearchArea MaxX")
		("MinY", po::value(&MinY)->default_value(0), "SearchArea MinY")
		("MaxY", po::value(&MaxY)->default_value(1), "SearchArea MaxY")
		("MinRot", po::value(&MinRot)->default_value(0), "SearchArea MinRot")
		("MaxRot", po::value(&MaxRot)->default_value(1), "SearchArea MaxRot")

		("SearchAreaDist", po::value(&SearchAreaDist)->default_value(0), "SearchArea Distance")
		("SearchAreaRot", po::value(&SearchAreaRot)->default_value(1), "SearchArea Rotation")

		("CORR_DISTANCE_TRESHOLD", po::value(&CORR_DISTANCE_TRESHOLD)->default_value(8), "CORR_DISTANCE_TRESHOLD")
		("FOV", po::value(&FOV)->default_value(8), "Laser Field of view")
		("NREADINGS", po::value(&NREADINGS)->default_value(8), "NREADINGS")
		("RELAXED_CORR_SEARCH", po::value<bool>(&RELAXED_CORR_SEARCH)->default_value("false"), "set true to get more correspondence points between scans by relaxing conditions")
	;

	po::variables_map vm;
	po::store(po::parse_config_file(config_stream, my_options, true), vm);
	po::notify(vm);
}
catch(std::exception& e)
{   std::cerr << "error: " << e.what() << "\n";
	return;
}
catch(...)
{   std::cerr << "Exception of unknown type!\n";
}


setPolarGAGeneticParameters(
		nbitx,
		nbity,
		nbitrot,
		popsize,
		maxruns,
		maxgen,
		pcross,
		pmutation);
setPolarGAParameters(NREADINGS, FOV, CORR_DISTANCE_TRESHOLD, RELAXED_CORR_SEARCH);
setPolarGASearchArea(MinX,MaxX,MinY,MaxY,MinRot,MaxRot);
#ifdef DRAW_PNG
setDebugParameters(draw_individual, draw_generations, verbose_level);
#endif

} // init


void pre_match(const scan &sref, const scan &snew)
{

	setPolarGARefScan(sref);
	setPolarGANewScan(snew);

	//setPolarGASearchArea( ws2.pos.x-SearchAreaDist/2.0,ws2.pos.x+SearchAreaDist/2.0,
	//                            ws2.pos.y-SearchAreaDist/2.0,ws2.pos.y+SearchAreaDist/2.0,
	//                            ws2.pos.rot-SearchAreaRot/2.0,ws2.pos.rot+SearchAreaRot/2.0);
}


matching_result match()
{
	matching_result mres;
	mres.status=0;
	mres.numiter=0;
	mres.p=PolarGA(&outfitness);
	return mres;
}

void post_match()
{
}

void deinit()
{
}





} // namespace
