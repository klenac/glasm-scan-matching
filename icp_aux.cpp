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


// icp_aux.cpp functions help interfacing to sm_cli simulator


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

//#include "utils.h"
#include "icp.h"

namespace ICP
{


scan ws1, ws2;

void pre_match(const scan &sref, const scan &snew)
{
	ws1.copy(sref);
	ws2.copy(snew);
}


matching_result match()
{	 matching_result mres;
	mres.status=0;
	mres.status=0;
	mres.p=icp(ws1, ws2, mres.numiter);
	return mres;
}

void post_match()
{
}

void deinit()
{
}


void init(const char* filename)
{

	bool draw_iterations;
	unsigned verbose_level;
	std::string map_filename;
	double map_size_x;
	double map_size_y;

	double corr_angle_treshold;
	double sensor_range;
	double stepangle;
	double maxfloat;
	double corr_distance_treshold;// distanza entro la quale cercare i corrispondenti
	unsigned max_iteration;
	double min_ddisp;
	double min_rdisp;
	unsigned corr_search_method;

	try {

	std::ifstream config_stream(filename);
	namespace po = boost::program_options;
	po::options_description my_options("Options");

	my_options.add_options()
	("draw_iterations", po::value<bool>(&draw_iterations)->default_value("false"), "set true to draw iterations")
	("verbose_level", po::value(&verbose_level)->default_value(0), "verbose_level")
	("map_filename", po::value<std::string>(&map_filename)->default_value("none"), "map filename")
	("map_size_x", po::value(&map_size_x)->default_value(15.0), "map_size_x")
	("map_size_y", po::value(&map_size_y)->default_value(15.0), "map_size_y")

	("corr_angle_treshold", po::value(&corr_angle_treshold)->default_value(0.148352986), "corr_angle_treshold")
	("sensor_range", po::value(&sensor_range)->default_value(8.0), "sensor_range")
	("stepangle", po::value(&stepangle)->default_value(0.015707963), "stepangle")
	("maxfloat", po::value(&maxfloat)->default_value(0.20), "maxfloat")
	("corr_distance_treshold", po::value(&corr_distance_treshold)->default_value(0.30), "corr_distance_treshold")
	("max_iteration", po::value(&max_iteration)->default_value(50), "max_iteration")
	("min_ddisp", po::value(&min_ddisp)->default_value(0.002), "min_ddisp")
	("min_rdisp", po::value(&min_rdisp)->default_value(0.001), "min_rdisp")
	("corr_search_method", po::value(&corr_search_method)->default_value(1), "corr_search_method")
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

setParameters(
		corr_angle_treshold,
		sensor_range,
		stepangle,
		maxfloat,
		corr_distance_treshold,
		max_iteration,
		min_ddisp,
		min_rdisp);
setCorrSearchMethod(corr_search_method);
setDebugParameters(draw_iterations,
		verbose_level,
		map_size_x,
		map_size_y,
		map_filename);


} // init

} // namespace
