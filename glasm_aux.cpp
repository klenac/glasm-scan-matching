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


// glasm_aux.cpp functions help interfacing to sm_cli simulator


#include <iostream>
#include <fstream>
#include <string>
//#include <stdlib.h>
//#include <stdio.h>
//#include <locale.h>
//#include <math.h>
//#include <sys/time.h>
#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "glasm_aux.h"

#ifdef DRAW_PNG
		#include "mydraw.h"
#endif

namespace GLASM
{

	double          outfitness;

	// todo: these parameters should be local to init(). Callbacks in glasm
	// should exist to set apropriate functions once during initialization
	// based on chosen type of matching, ie based on these ...
	// if's and namespace pollution would be avoided
	bool        global_search;
	bool        static_map;
	bool        gradient;                           // gradient based lookup generation

	std::string static_map_file;
	std::string static_map_valid_pos_file;


void init(const char* filename)
{
// reads and sets initial parameters
//#ifdef DRAW_PNG
	bool draw_lookup;
	bool draw_individual;
	bool draw_generations;
	unsigned verbose_level;
	std::string map_filename;
	double map_size_x;
	double map_size_y;
//#endif

	int         nbitx;
	int         nbity;
	int         nbitrot;
	double      pcross;
	double      pmutation;
	int         popsize;
	int         maxruns;
	int         maxgen;
	double      SearchAreaDist;
	double      SearchAreaRot;

	double      MinX;
	double      MaxX;
	double      MinY;
	double      MaxY;
	double      MinRot;
	double      MaxRot;

	unsigned    lookup_rows;                        // number of rows= 3000;
	unsigned    lookup_columns;                     // number of columns= 3000;
	double      lookup_size_x;                      // dimensione lato x area in metri =10.0
	double      lookup_size_y;                      // dimensione lato y area in metri =10.0
	double      corr_distance;                      // distanza entro la quale cercare i corrispondenti
	double      lookup_ox;                          // coord x dell'angolo inf sin dell'area a ScanRef.pos.x-lookup_size_x/2.0
	double      lookup_oy;                          // ScanRef.pos.y-lookup_size_y/2.0
	double      border_size;                        // border_size in meters added around static lookup table


	double      RGmean;
	double      RGvariance;
	double      RGextension;
	double      RGweight;
	double      RGscale;
	double      RGbitmapsize;


	try {

		std::ifstream config_stream(filename);
		namespace po = boost::program_options;
		po::options_description my_options("Options");

		my_options.add_options()
//#ifdef DRAW_PNG
		("draw_lookup", po::value<bool>(&draw_lookup)->default_value("false"), "set true to draw_lookup table ")
		("draw_individual", po::value<bool>(&draw_individual)->default_value("false"), "set true to draw_individual")
		("draw_generations", po::value<bool>(&draw_generations)->default_value("false"), "set true to draw_generations ")
		("map_filename", po::value<std::string>(&map_filename)->default_value("none"), "map filename")
	("map_size_x", po::value(&map_size_x)->default_value(15.0), "map_size_x")
	("map_size_y", po::value(&map_size_y)->default_value(15.0), "map_size_y")
		("verbose_level", po::value(&verbose_level)->default_value(0), "verbose_level")
//#endif
		("global_search", po::value<bool>(&global_search)->default_value("false"), "set true to enable global scan matching, ie without initial position estimate to start search in")
		("static_map", po::value<bool>(&static_map)->default_value("false"), "set true to create a lookup table from a-priori known map (bitmap image). Lookup is created only once during initialization.")
		("static_map_file", po::value<std::string>(&static_map_file)->default_value("map.png"), "static_map filename - full path")
		("static_map_valid_pos_file", po::value<std::string>(&static_map_valid_pos_file)->default_value("map.png"), "static_map_valid_pos_file filename - full path")

			// search area (specify either exact range or by distance and rotation. the former takes precedence if both specified)
		("MinX", po::value(&MinX)->default_value(0), "SearchArea MinX")
		("MaxX", po::value(&MaxX)->default_value(1), "SearchArea MaxX")
		("MinY", po::value(&MinY)->default_value(0), "SearchArea MinY")
		("MaxY", po::value(&MaxY)->default_value(1), "SearchArea MaxY")
		("MinRot", po::value(&MinRot)->default_value(0), "SearchArea MinRot")
		("MaxRot", po::value(&MaxRot)->default_value(1), "SearchArea MaxRot")

		("SearchAreaDist", po::value(&SearchAreaDist)->default_value(0), "SearchArea Distance")
		("SearchAreaRot", po::value(&SearchAreaRot)->default_value(1), "SearchArea Rotation")

		// genetic parameters
		("nbitx", po::value(&nbitx)->default_value(8), "nbitx")
		("nbity", po::value(&nbity)->default_value(8), "nbity")
		("nbitrot", po::value(&nbitrot)->default_value(8), "nbitrot")
		("pcross", po::value(&pcross)->default_value(0.9), "pcross")
		("pmutation", po::value(&pmutation)->default_value(0.012), "pmutation")
		("popsize", po::value(&popsize)->default_value(100), "popsize")
		("maxruns", po::value(&maxruns)->default_value(1), "maxruns")
		("maxgen", po::value(&maxgen)->default_value(10), "maxgen")

		// lookup table parameters
		("gradient", po::value<bool>(&gradient)->default_value("false"), "set true to enable gradient based lookup")
		("lookup_rows", po::value(&lookup_rows)->default_value(3000), "lookup_rows")
		("lookup_columns", po::value(&lookup_columns)->default_value(3000), "lookup_columns")
		("lookup_size_x", po::value(&lookup_size_x)->default_value(40.0), "lookup_size_x")
		("lookup_size_y", po::value(&lookup_size_y)->default_value(40.0), "lookup_size_y")
		("corr_distance", po::value(&corr_distance)->default_value(0.20), "corr_distance")
		("lookup_ox", po::value(&lookup_ox)->default_value(-20.0), "lookup_ox")
		("lookup_oy", po::value(&lookup_oy)->default_value(-20.0), "lookup_oy")
		("border_size", po::value(&border_size)->default_value(2.0), "border_size in meters added around static lookup table")


		// parameters for radial gradient method
		("RGmean", po::value(&RGmean)->default_value(8), "Radial Gradient normal curve mean")
		("RGvariance", po::value(&RGvariance)->default_value(8), "Radial Gradient normal curve variance")
		("RGextension", po::value(&RGextension)->default_value(8), "Radial Gradient  extension of normal curve over the image - good choice is RGbitmapsize*RGscale*1.2")
		("RGweight", po::value(&RGweight)->default_value(8), "Weight given to normalized normal curve when drawing Radial Gradient")
		("RGscale", po::value(&RGscale)->default_value(8), "Scale, ie, width of the normalized normal curve to map to Radial Gradient bitmap")
		("RGbitmapsize", po::value(&RGbitmapsize)->default_value(8), "RGBitmapsize is equivalent to GLASM_CORR_DISTANCE_TRESHOLD in the old binary GLASM, ie area around ref points in lookup")
		;

		po::variables_map vm;
		po::store(po::parse_config_file(config_stream, my_options, true), vm);
		po::notify(vm);
		// option values are now in their variables, also accessible by vm["SectionName.option"].as<int>()
	}
	catch(std::exception& e)
	{   std::cerr << "error: " << e.what() << "\n";
	return;
	}
	catch(...)
	{   std::cerr << "Exception of unknown type!\n";
	}

	//std::cout << filename << "glasm aux map:" <<map_filename<< " "<< map_size_x << " " << map_size_y << std::endl;

//#ifdef DRAW_PNG
	setDebugParameters(draw_lookup,
		draw_individual,
		draw_generations,
		verbose_level,
		map_size_x,
		map_size_y,
		map_filename,
		static_map_file,
		static_map_valid_pos_file);
//#endif

	// set glasm parameters
	initialize_invGray_table();
//    std::cout << "[glasm.init                ] setting genetic parameters:" <<
//            nbitx << ":" <<
//            nbity << ":" <<
//            nbitrot << ":" <<
//            popsize << ":" <<
//            maxruns << ":" <<
//            maxgen << ":" <<
//            pcross << ":" <<
//            pmutation << std::endl;
	setGeneticParameters(
			nbitx,
			nbity,
			nbitrot,
			popsize,
			maxruns,
			maxgen,
			pcross,
			pmutation);
	// todo: where to center the lookup table: previously always in 0,0. Now? Dynamically?
//    std::cout << "[glasm.init                ] setting lookup parameters:" << lookup_rows << ":" <<
//            lookup_columns << ":" <<
//            lookup_size_x << ":" <<
//            lookup_size_y << ":" <<
//            corr_distance << std::endl;
	setLookupTableParameters(lookup_rows,lookup_columns,lookup_size_x,lookup_size_y,0,0,corr_distance); // only params, the lookup table will be recalculated each step
//    std::cout << "[glasm.init                ] setting search area:" <<
//            " x from " << MinX << " to " << MaxX <<
//            "  y from " << MinY << " to " << MaxY <<
//            "  rot from " << MinRot << " to " << MaxRot << std::endl;
	setSearchArea(MinX,MaxX,MinY,MaxY,MinRot,MaxRot);


	if (gradient)
	{   //std::cout << "[glasm.init                ] lookup type is gradient" << std::endl;
		initialize_radial_gradient(RGmean,RGvariance,RGextension,RGweight,RGscale,RGbitmapsize); // RGbitmapsize is equivalent to GLASM_CORR_DISTANCE_TRESHOLD in the old GLASM
	} else
	{   //std::cout << "[glasm.init                ] lookup type is binary" << std::endl;
	}
	if (static_map)
	{   //std::cout << "[glasm.init                ] lookup created from static_map" << std::endl;
		initialize_lookup_table_from_bitmap(static_map_file.c_str(), lookup_size_x, lookup_size_y, border_size);
		initialize_valid_table_from_bitmap(static_map_valid_pos_file.c_str());
	} else
	{   //std::cout << "[glasm.init                ] lookup will be created from sref" << std::endl;
	}
}

void pre_match(const scan &sref, const scan &snew)
{
	// temp
	// todo: check why it works if both scans are set in origin (global) ref frame
	// but it should also work if both scans are in same ref frame (not necesserally  global ref frame)

// p1=sref.pos;
// p2=snew.pos;

//    temp
//    scan s1,s2;
//    s1.copy(sref);
//    scanChangeRefFrame(s1,{0.0,0.0,0.0});
//    s1.readings.sort(readingsAngleCompare);
//
//    s2.copy(snew);
//    scanChangeRefFrame(s2,{0.0,0.0,0.0});
//    s2.readings.sort(readingsAngleCompare);
//    setRefScan(s1);
//    setNewScan(s2);
//    end temp


//std::cout << "[pre_match] map_filename: " << map_filename << std::endl;


	setRefScan(sref);
	setNewScan(snew);
	// end temp

	if (!static_map) {
		if (gradient)   initialize_gradient_lookup_table();
		else            initialize_binary_lookup_table();
	}
	if (global_search)  setCenterOfSearchArea(sref.pos);
	else                setCenterOfSearchArea(snew.pos);

}


matching_result match()
{	matching_result mres;
	mres.status=0;
	mres.numiter=0;
	mres.p=glasm(&outfitness);
	return mres;
}

void post_match()
{   if (!static_map) {
		delete_lookup_table();
	}
}

void deinit()
{   if (static_map) delete_lookup_table();
}






} // namespace GLASM
