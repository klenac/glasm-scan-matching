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

/*
 * t.cpp
 *
 *  Created on: Jul 15, 2010
 *      Author: Kristijan Lenac
 */


#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <dirent.h>

//#include <fstream>

#define ERRvsITER

#include "simulation.h"
#include "mydraw.h"
#ifdef DRAW_PNG
	#include "mydraw.h"
#endif

#include "utils.h"

#include "hglasm_aux.h"
#include "glasm_aux.h"
#include "polarga_aux.h"
#include "pgambicp_aux.h"
#include "icp_aux.h"
#include "mbicp_aux.h"
#include "gcp_aux.h"

// empty algorithm functions
namespace EMPTY
{
	scan s1;
	void init(const char* filename) {}					// one time initialization. Reads and sets parameters.
	void pre_match(const scan &sref, const scan &snew)	// before matching
	{	s1.copy(snew);
	}
	matching_result match()								// matching
	{	matching_result mres;
		mres.status=0;
		mres.numiter=0;
		mres.p=s1.pos;
		return mres;
	}
	void post_match() {}                                // after matching
	void deinit() {}                                    // one time deinitialization
}




namespace po = boost::program_options;

// an array with alg_info structures is filled from simulation parameters file
// for example look in examples/sm.ini file (each experiment should be defined in one
// such file)

struct alg_info {
	std::string name;
	std::string base_algorithm;
	std::string parameters_file;
	double mtime_success;
	double vtime_success;
	double mtime_failure;
	double vtime_failure;
	double mdist;
	double mrot;
	double vdist;
	double vrot;
	double miter_success;
	double miter_failure;
	double viter_success;
	double viter_failure;
	unsigned successes;
	unsigned failures;
	unsigned false_positives;
	unsigned true_negatives;
	unsigned false_negatives;
	unsigned true_positives;

	// c style callback functions. TODO: rewrite this in C++ style
	// pointers to base_algorithm dependent functions:
	void (*init)(const char *);	// read algorithm parameters from file and initialize algorithm
	void (*pre_match)(const scan &snew, const scan &sref); // before matching
	matching_result (*match)();	// matching
	void (*post_match)(); // after matching
	void (*deinit)(); // one time deinitialization
};

unsigned  NumAlgorithms; // number of algorithms to compare

std::vector<alg_info> algorithm;

/********************************
 * VARIABLES
 ********************************/

// general
scan			ts1,ts2;	// true scans, read from file at beginning of a test they never change afterwards
scan			ws1,ws2;	// working scans
scan			sws1;		// used to store ws1 which must be restored before each test (algorithms may alter ws1)
scan			sws2;		// used to store ws1 which must be restored before each test (algorithms may alter ws1)
scan            rs;			// holds the resulting scan for drawing
matching_result mres;		// resulting position after the matching

double ex;
double ey;
double er;
const double estep=0.05;

/********************************
 * PARAMETERS
 * (read from parameters file or
 * cmd line options)
 ********************************/

// general
const char *general_parameters_file = "experiments/current/sm.ini";

std::string     map_filename;
double          map_size_x;
double          map_size_y;
std::string     scan_path;
std::string     scan_basename;
std::vector<std::string> rawseeds_scans;
std::string     scan_format;
std::string     ref_scan_filename;
std::string     new_scan_filename;
int             new_scan_offset;	// when chosen with offset method, new scan is loaded from file
									// with index at this offset from ref scan
unsigned        loop_start;			// to test one scan pair set loop_start and loop_end to 1
unsigned        loop_end;			// to test multiple scan pairs, ref scan is loaded starting from
									// index loop_start, the loop repeats until loop_end.
unsigned        runs;				// number of repetitions of the same test -> this is only to
									// average computation time if we have few tests only.
unsigned        skip_first_runs;	// we can skip the initial runs which may be perturbed until OS
									// activity settles ...(typ 5)

bool            draw_scans;
bool            verbose_level;

unsigned        ref_scan_color_R;
unsigned        ref_scan_color_G;
unsigned        ref_scan_color_B;

unsigned        new_scan_color_R;
unsigned        new_scan_color_G;
unsigned        new_scan_color_B;

unsigned        res_scan_color_R;
unsigned        res_scan_color_G;
unsigned        res_scan_color_B;

unsigned        min_scan_readings;	// minimum number of readings that should be in the scans before
									// matching to proceed with matching test

bool            resample_scans;
double          resample_distance;
bool            median_filter_scans;
unsigned        median_filter_order;
double          median_filter_distance;

bool            ref_scan_odd_readings_only;
bool            ref_scan_even_readings_only;
bool            new_scan_odd_readings_only;
bool            new_scan_even_readings_only;

std::string     fakeError;
double          maxFakeErrorRot;
double          maxFakeErrorX;
double          maxFakeErrorY;

double          success_treshold_rotation;
double          success_treshold_traslation;


void read_simulation_parameters()
{
	const unsigned MAXALGS = 10;                 // max number of algorithms to compare at a time
	alg_info a[MAXALGS];

	if (verbose_level>0)
	{	std::cout << "[read_simulation_parameters] reading simulation parameters based on file: " << general_parameters_file << std::endl;
	}
	try {
		std::ifstream config_stream(general_parameters_file);
		namespace po = boost::program_options;
		po::options_description my_options("Options");

		my_options.add_options()
			("general.map_filename", po::value<std::string>(&map_filename)->default_value("none"), "map filename")
			("general.map_size_x", po::value(&map_size_x)->default_value(8.0), "map size in meters (x axis)")
			("general.map_size_y", po::value(&map_size_y)->default_value(8.0), "map size in meters (y axis)")
			("general.ref_scan", po::value<std::string>(&ref_scan_filename)->default_value("scan3"), "reference scan filename")
			("general.new_scan", po::value<std::string>(&new_scan_filename)->default_value("scan4"), "new scan filename")
			("general.scan_path", po::value<std::string>(&scan_path)->default_value("./scans/"), "path to folder containing scans")
			("general.scan_basename", po::value<std::string>(&scan_basename)->default_value("scan"), "base name for scan file/s. For sm_cli format basename is the base name"
				" for files like <basename>xxx where xxx is unsigned int, For rawseeds format basename is the name of the file containing scans")
			("general.scan_format", po::value<std::string>(&scan_format)->default_value("sm_cli"), "format used to save the scans in file.")
			("general.new_scan_offset", po::value(&new_scan_offset)->default_value(1), "new scan is read from file new_scan_offset positions ahead of ref scan")
			("general.loop_start", po::value(&loop_start)->default_value(1), "loop start")
			("general.loop_end", po::value(&loop_end)->default_value(1), "loop end")
			("general.runs", po::value(&runs)->default_value(1), "number of repetitions of the same test")
			("general.skip_first_runs", po::value(&skip_first_runs)->default_value(1), "number of initial repetitions of the same test to skip")
			("general.draw_scans", po::value<bool>(&draw_scans)->default_value("true"), "set true to enable drawing of scans")
			("general.verbose_level", po::value(&verbose_level)->default_value(0), "verbose_level")
			("general.ref_scan_color_R", po::value(&ref_scan_color_R)->default_value(0), "color to draw ref_scan")
			("general.ref_scan_color_G", po::value(&ref_scan_color_G)->default_value(0), "color to draw ref_scan")
			("general.ref_scan_color_B", po::value(&ref_scan_color_B)->default_value(30000), "color to draw ref_scan")
			("general.new_scan_color_R", po::value(&new_scan_color_R)->default_value(0), "color to draw ref_scan")
			("general.new_scan_color_G", po::value(&new_scan_color_G)->default_value(0), "color to draw ref_scan")
			("general.new_scan_color_B", po::value(&new_scan_color_B)->default_value(60000), "color to draw ref_scan")
			("general.res_scan_color_R", po::value(&res_scan_color_R)->default_value(0), "color to draw ref_scan")
			("general.res_scan_color_G", po::value(&res_scan_color_G)->default_value(60000), "color to draw ref_scan")
			("general.res_scan_color_B", po::value(&res_scan_color_B)->default_value(0), "color to draw ref_scan")
			("general.min_scan_readings", po::value(&min_scan_readings)->default_value(8), "minimum number of readings that should be in the scans before matching to proceed with matching test")
			("general.resample_scans", po::value<bool>(&resample_scans)->default_value("true"), "set true to enable resampling of scans ")
			("general.resample_distance", po::value(&resample_distance)->default_value(0.03), "resample distance")
			("general.median_filter_scans", po::value<bool>(&median_filter_scans)->default_value("true"), "set true to enable filtering of scans with median filter")
			("general.median_filter_distance", po::value(&median_filter_distance)->default_value(0.06), "median filter distance")
			("general.median_filter_order", po::value(&median_filter_order)->default_value(3), "median filter order")
			("general.ref_scan_odd_readings_only", po::value<bool>(&ref_scan_odd_readings_only)->default_value("false"), "set true to take only every second reading starting from reading 1, thus reducing ref_scan number of readings by half")
			("general.ref_scan_even_readings_only", po::value<bool>(&ref_scan_even_readings_only)->default_value("false"), "set true to take only every second reading starting from reading 1, thus reducing ref_scan number of readings by half")
			("general.new_scan_odd_readings_only", po::value<bool>(&new_scan_odd_readings_only)->default_value("false"), "set true to take only every second reading starting from reading 1, thus reducing new_scan number of readings by half")
			("general.new_scan_even_readings_only", po::value<bool>(&new_scan_even_readings_only)->default_value("false"), "set true to take only every second reading starting from reading 1, thus reducing new_scan number of readings by half")
			("general.fakeError", po::value<std::string>(&fakeError)->default_value("random"), "type of fake Error to introduce in new scan position")
			("general.maxFakeErrorX", po::value(&maxFakeErrorX)->default_value(0.2), "maximum added error (initial perturbation) in X direction (meters)")
			("general.maxFakeErrorY", po::value(&maxFakeErrorY)->default_value(0.2), "maximum added error (initial perturbation) in X direction (meters)")
			("general.maxFakeErrorRot", po::value(&maxFakeErrorRot)->default_value(0.2), "maximum added error (initial perturbation) in Rot (radians)")
			("general.success_treshold_traslation", po::value(&success_treshold_traslation)->default_value(0.1), "A matching is considered a failure when the solution is larger than this value in translation (meters)")
			("general.success_treshold_rotation", po::value(&success_treshold_rotation)->default_value(0.01), "A matching is considered a failure when the solution is larger than this value in rotation (rad)")
			;

		unsigned i=0;                                // this part of code is ugly but it's the first thing that came to my mind
		while (i<MAXALGS) {                          //  to add options for a group of algorithms in ini file
			const unsigned NAMESIZE=255;
			char name[NAMESIZE];
			char base_algorithm[NAMESIZE];
			char parameters_file[NAMESIZE];
			char base[NAMESIZE];
			memset(name,0,NAMESIZE);
			memset(base_algorithm,0,NAMESIZE);
			memset(parameters_file,0,NAMESIZE);
			memset(base,0,NAMESIZE);
			strncpy(base,"algorithm_",10);
			sprintf(&base[10],"%d",i);
			strcpy(name,base);
			strcpy(&name[strlen(name)],".name"); // name="algorithm_i.name"
			strcpy(base_algorithm,base);
			strcpy(&base_algorithm[strlen(base_algorithm)],".base_algorithm"); // base_algorithm="algorithm_i.base_algorithm"
			strcpy(parameters_file,base);
			strcpy(&parameters_file[strlen(parameters_file)],".parameters_file"); // parameters_file="algorithm_i.parameters_file"

			my_options.add_options()
				(name, po::value<std::string>(&(a[i].name))->default_value("-1"), "algorithm name")
				(base_algorithm, po::value<std::string>(&(a[i].base_algorithm))->default_value("-1"), "algorithm from which this one is derived")
				(parameters_file, po::value<std::string>(&(a[i].parameters_file))->default_value("-1"), "parameters_file from which to load algorithm parameters")
				;

			i++;
		} // end while

		po::variables_map vm;
		po::store(po::parse_config_file(config_stream, my_options), vm);
		po::notify(vm);
		// option values are now in their variables, also accessible by vm["SectionName.option"].as<int>()
	} // end try
	catch(std::exception& e)
	{   std::cerr << "error: " << e.what() << "\n";
		return;
	}
	catch(...)
	{   std::cerr << "Exception of unknown type!\n";
	}

	if (verbose_level>0)
	{	std::cout << "[read_simulation_parameters] simulation parameters read" << std::endl;
		//std::cout << "[read_simulation_parameters] map_filename: " << map_filename << std::endl;
	}


	NumAlgorithms=0;                             // (set after parsing program options, currenty max 10 algorithms at a time)
	algorithm.resize(MAXALGS);
	unsigned i=0;
	while (i<MAXALGS)
	{	/*************************************
		*
		* SUPPORTED BASE ALGORITHMS (as of 20110515)
		* - glasm
		* - icp
		* - gcp     [todo enable it here ..]
		* - polarga
		* - mbicp
		* - wrsm    [todo enable it here ..]
		* - heglasm [todo enable it here ..]
		* - empty
		*
		* - to add a new type of algorithm define callback functions in
		*  newalgname_aux[.cpp .h] files and associate them here
		*
		*************************************/


		//std::cout << " a["<<i<<"].name:"<<a[i].name << std::endl;
		//std::cout << " a["<<i<<"].base_algorithm:"<<a[i].base_algorithm << std::endl;
		//std::cout << " a["<<i<<"].parameters_file:"<<a[i].parameters_file << std::endl;

		// associate the right functions to call (depending on the base algorithm)
		if (a[i].base_algorithm.compare("empty")==0 || a[i].base_algorithm.compare("\"empty\"")==0)
		{   algorithm[NumAlgorithms].init            = EMPTY::init;
			algorithm[NumAlgorithms].pre_match       = EMPTY::pre_match;
			algorithm[NumAlgorithms].match           = EMPTY::match;
			algorithm[NumAlgorithms].post_match      = EMPTY::post_match;
			algorithm[NumAlgorithms].deinit          = EMPTY::deinit;
		} else
		if (a[i].base_algorithm.compare("glasm")==0 || a[i].base_algorithm.compare("\"glasm\"")==0)
		{   algorithm[NumAlgorithms].init            = GLASM::init;
			algorithm[NumAlgorithms].pre_match       = GLASM::pre_match;
			algorithm[NumAlgorithms].match           = GLASM::match;
			algorithm[NumAlgorithms].post_match      = GLASM::post_match;
			algorithm[NumAlgorithms].deinit          = GLASM::deinit;
		} else
		if (a[i].base_algorithm.compare("icp")==0 || a[i].base_algorithm.compare("\"icp\"")==0)
		{
			algorithm[NumAlgorithms].init            = ICP::init;
			algorithm[NumAlgorithms].pre_match       = ICP::pre_match;
			algorithm[NumAlgorithms].match           = ICP::match;
			algorithm[NumAlgorithms].post_match      = ICP::post_match;
			algorithm[NumAlgorithms].deinit          = ICP::deinit;
		} else
		if (a[i].base_algorithm.compare("mbicp")==0 || a[i].base_algorithm.compare("\"mbicp\"")==0)
		{
			algorithm[NumAlgorithms].init            = MBICP::init;
			algorithm[NumAlgorithms].pre_match       = MBICP::pre_match;
			algorithm[NumAlgorithms].match           = MBICP::match;
			algorithm[NumAlgorithms].post_match      = MBICP::post_match;
			algorithm[NumAlgorithms].deinit          = MBICP::deinit;
		} else
		if (a[i].base_algorithm.compare("hglasm")==0 || a[i].base_algorithm.compare("\"hglasm\"")==0)
		{   algorithm[NumAlgorithms].init            = HGLASM::init;
			algorithm[NumAlgorithms].pre_match       = HGLASM::pre_match;
			algorithm[NumAlgorithms].match           = HGLASM::match;
			algorithm[NumAlgorithms].post_match      = HGLASM::post_match;
			algorithm[NumAlgorithms].deinit          = HGLASM::deinit;
		} else
		if (a[i].base_algorithm.compare("pgambicp")==0 || a[i].base_algorithm.compare("\"pgambicp\"")==0)
		{   algorithm[NumAlgorithms].init			= PGAMBICP::init;
			algorithm[NumAlgorithms].pre_match		= PGAMBICP::pre_match;
			algorithm[NumAlgorithms].match			= PGAMBICP::match;
		algorithm[NumAlgorithms].post_match			= PGAMBICP::post_match;
		algorithm[NumAlgorithms].deinit				= PGAMBICP::deinit;
		} else
		if (a[i].base_algorithm.compare("polarga")==0 || a[i].base_algorithm.compare("\"polarga\"")==0)
		{
			algorithm[NumAlgorithms].init			= POLARGA::init;
			algorithm[NumAlgorithms].pre_match		= POLARGA::pre_match;
			algorithm[NumAlgorithms].match			= POLARGA::match;
			algorithm[NumAlgorithms].post_match		= POLARGA::post_match;
			algorithm[NumAlgorithms].deinit			= POLARGA::deinit;
		}
		if (!a[i].base_algorithm.compare("-1")==0) {	// if there was an algorithm_i in ini file

			algorithm[NumAlgorithms].name 			= a[i].name;
			algorithm[NumAlgorithms].base_algorithm 	= a[i].base_algorithm;
			algorithm[NumAlgorithms].parameters_file = a[i].parameters_file;
			// todo: trim a[i].parameters_file from surrounding "". For now do not surround filename with "" in sm.ini
			if (verbose_level>0)
			{	std::cout << "[read_simulation_parameters] added algorithm: " << algorithm[NumAlgorithms].name
				<< " (based on " << algorithm[NumAlgorithms].base_algorithm << ")"
				<< " with parameters in: " << algorithm[NumAlgorithms].parameters_file
				<< std::endl;
			}


			NumAlgorithms++;
		}

		i++;
	}
}


int filecount(const char * path)
{	int file_count = 0;
	DIR * dirp;
	struct dirent * entry;

	dirp = opendir(path); /* There should be error handling after this */
	while ((entry = readdir(dirp)) != NULL)
	{	if (entry->d_type == DT_REG)  /* If the entry is a regular file */
		{	file_count++;
		}
	}
	closedir(dirp);
	return file_count;
}

void rawseeds_write_random_sample(const char *filename, unsigned N)
// take N random scans from the vector rawseeds_scans and save them to disc in rawseeds format
{	unsigned TOTALSCANS=rawseeds_scans.size();
	std::cout << TOTALSCANS << std::endl;
	std::ofstream file(filename);
	for(unsigned int i=0; i<N; i++)
	{	unsigned r=rand() % TOTALSCANS;
		std::cout << " i:" << i << " r:" << r << std::endl;
		if (i>=TOTALSCANS)
		{	std::cout << "rawseeds_scans size is less than required!"
			"The number of saved scans is less than requested.";
			break;
		}
		if (verbose_level>0)
		{	std::cout << rawseeds_scans[r] << std::endl;
		}
		file << rawseeds_scans[r] << std::endl;
	}
	file.close();
}

void load_scan_pair(unsigned pair)
{
/*
Determine which scans to load from disk based on content of ref_scan_filename and
new_scan_filename parameters.
2 file formats are supported: the original one scan one file format and the Rawseeds project format
Depending of the parameters in sm.ini many variants for choosing the scan pair are implemented.
Look at the source for all supported variants. Some possibilities for choosing ref and new scans:
ref_scan_filename format: { "sequence" | "sequence_repeat_8" | "random" | "random_repeat_8" | "<N>" }
new_scan_filename format: { "random" | "offset" | "<N>" }
For fixed scan "<N>" put loop_start==loop_end, otherwise you just repeat the same test, wasting time
To match a scan with itself put offset = 0

Example: if
	scan_path= "./scans/"
	scan_basename= "scan"
	ref_scan_filename = "15"
	new_scan_filename = "offset"
	new_scan_offset = 1
then the files would be: "./scans/scan15" for reference scan and
							"./scans/scan16" for new scan
*/
	if (verbose_level>0) std::cout << "---" <<std::endl;

	std::string sr = ref_scan_filename;
	std::string sn = new_scan_filename;
	unsigned N;
	if (scan_format.compare("sm_cli")==0)
	{	N=filecount(scan_path.c_str());
		if (verbose_level>0)
		{	std::cout << "[load_scan_pair            ] Number of scans found in " << scan_path << " folder: " << N << std::endl;
			std::cout << "[load_scan_pair            ] ";
		}
		// ref_scan_filename
		if (sr.compare("random") == 0)
		{   unsigned r=rand() % N;
			if (r+new_scan_offset>=N) r-=new_scan_offset;  // to be sure that new scan will be in too
			std::stringstream st;
			st << r; // convert to string
			sr=scan_path+scan_basename+st.str();
			if (verbose_level>0) std::cout << "using random ref_scan: " << sr<< " : ";
		}
		else if (sr.compare("random_repeat_8") == 0)
		{	static unsigned r=rand() % N;
			if (verbose_level>0)  std::cout << "in randomrepeat8 pair%8==0 " << pair%8 << " pair "<< pair << std::endl;
			if (pair%8==0) {
				r=rand() % N;
				if (r+new_scan_offset>=N) r-=new_scan_offset;  // to be sure that new scan will be in too
			}
			std::stringstream st;
			st << r; // convert to string
			sr=scan_path+scan_basename+st.str();
			if (verbose_level>0) std::cout << "using random_repeat_8 ref_scan: " << sr<< " : ";
		}
		else if (sr.compare("sequence_repeat_8") == 0)
		{	std::stringstream st;
			st << int(pair/8.0); // convert to string
			if (verbose_level>0)  std::cout << "in sequencerepeat8 pair:" << pair << " st"<< st << std::endl;
			sr=scan_path+scan_basename+st.str();
			if (verbose_level>0)  std::cout << "using next in sequence repeat ref_scan: " << sr<< " : ";
		}
		else if (sr.compare("sequence") == 0)
		{   std::stringstream st;
			st << pair; // convert to string
			sr=scan_path+scan_basename+st.str();
			if (verbose_level>0)  std::cout << "using next in sequence ref_scan: " << sr<< " : ";
		}
		else
		{   // ref_scan is fixed
			sr=scan_path+scan_basename+sr;
			//std::cout << "scan_path: " << scan_path << " : "<< "scan_basename: " << scan_basename << " : ";
			if (verbose_level>0)  std::cout << "using ref_scan: " << sr<< " : ";
		}

		// new_scan_filename
		if (sn.compare("offset") == 0)
		{   // new_scan is at some offset from ref_scan
			// find ref_scan index
			std::string s=sr;
			s.erase(0,scan_path.length()+scan_basename.length());
			unsigned i = atoi(s.c_str());
			std::stringstream st;    // convert to string
			st << (i+new_scan_offset);
			sn=scan_path+scan_basename+st.str();
			if (verbose_level>0)  std::cout << "using offset new_scan: " << sn<<std::endl;;
		}
		else
		{   // new_scan is fixed
			sn=scan_path+scan_basename+sn;
			if (verbose_level>0)  std::cout << "using new_scan: " << sn<<std::endl;;
		}

		// read scans
		//
		// typically ground truth is known - it is saved in scan.pos in scans on disk
		// (be sure to save the ground truth position when preparing folder with scans)
		// for simulated scans this is easy. For real scans ground truth position
		// is not known. In this case match against same scan (odd even readings) or establish
		// ground truth with other method.
		ts1.read(sr.c_str());
		ts2.read(sn.c_str());

	} else if (scan_format.compare("rawseeds") == 0)
	{
		static bool file_loaded = false;
		if (!file_loaded)
		{	if (verbose_level>0) std::cout << "Loading from rawseeds format file: " << scan_path+scan_basename << std::endl;
			rawseeds_scans.clear();
			std::ifstream file((scan_path+scan_basename).c_str());
			std::string mystr;
			long nn=0;
			while (getline(file,mystr) && nn<200000)
			{	rawseeds_scans.push_back(mystr);
				nn++;
			}
			if (verbose_level>0)  std::cout << "Loaded " << nn << " lines." << std::endl;
			file_loaded=true;
		}

		N=rawseeds_scans.size();

		//std::cout << "[load_scan_pair            ] Number of scans found in " << scan_path << scan_basename << ": " << N << std::endl;
		if (verbose_level>0)  std::cout << "[load_scan_pair            ] ";

		static unsigned r=0;
		static unsigned n=0;

		// ref_scan
		if (sr.compare("random") == 0)
		{	r=rand() % N;
			if (r+new_scan_offset>=N) r-=new_scan_offset;  // to be sure that new scan will be in too
			if (verbose_level>0)  std::cout << "using random ref_scan: " << r<< " : ";
		}
		else if (sr.compare("random_repeat_8") == 0)
		{	if (verbose_level>0) std::cout << "in randomrepeat8 pair%8==0 " << pair%8 << " pair "<< pair << std::endl;
			if (pair%8==0)
			{	r=rand() % N;
				if (r+new_scan_offset>=N) r-=new_scan_offset;  // to be sure that new scan will be in too
			}
			if (verbose_level>0) std::cout << "using random_repeat_8 ref_scan: " << r<< " : ";
		}
		else if (sr.compare("sequence_repeat_8") == 0)
		{	r= int(pair/8.0);
			if (verbose_level>0) std::cout << "using next in sequence repeat ref_scan: " << r << " : ";
		}
		else if (sr.compare("sequence") == 0)
		{	r=pair;
			if (verbose_level>0) std::cout << "using next in sequence ref_scan: " << r<< " : ";
		}
		else
		{	r=atoi(sr.c_str());
			if (verbose_level>0) std::cout << "using ref_scan: " << r<< " : ";
		}

		// new_scan
		if (sn.compare("offset") == 0)
		{   // new_scan is at some offset from ref_scan
			n=r+new_scan_offset;
			if (verbose_level>0) std::cout << "using offset new_scan: " << n<<std::endl;
		}
		else
		{   // new_scan is fixed
			n=atoi(sn.c_str());
			if (verbose_level>0) std::cout << "using new_scan: " << n<<std::endl;
		}

		// read scans
		if (r<rawseeds_scans.size())
		{
			ts1.read(rawseeds_scans[r]);
			ts2.read(rawseeds_scans[n]);
			//std::cout << "read scan ts1 size: " << ts1.readings.size()<<std::endl;
			//std::cout << "read scan ts2 size: " << ts2.readings.size()<<std::endl;
		}
		else
		{	if (verbose_level>0) std::cout << "r > rawseeds_scans size: " << rawseeds_scans.size()<<std::endl;
			exit(1);
		}

	} // if rawseeds scan format

	if (ref_scan_odd_readings_only)  scanDecimate(ts1,2,1); else
	if (ref_scan_even_readings_only) scanDecimate(ts1,2,0);

	if (new_scan_odd_readings_only)  scanDecimate(ts2,2,1); else
	if (new_scan_even_readings_only) scanDecimate(ts2,2,0);
	//std::cout << "read scan ts1 size: " << ts1.readings.size()<<std::endl;
	//std::cout << "read scan ts2 size: " << ts2.readings.size()<<std::endl;
} // load_scan_pair


void add_fake_error_pos(scan &s)
{
/* add fake error - initial error that we add to new scan to test robustness
 * The error, depending on chosen method, may span between [-maxFakeError,maxFakeError]
 * There are two ways to do add the error:
 * 1) change scan.pos
 * 2) change every reading keeping scan.pos intact
 * This procedure uses the first simple way.
 */
	if (verbose_level>0) std::cout << "[add_fake_error            ] ";
	if (fakeError.compare("random") == 0)
	{	s.pos.rot=angleUnwrap(s.pos.rot+randDouble(-maxFakeErrorRot,maxFakeErrorRot));
		s.pos.x=s.pos.x+randDouble(-maxFakeErrorX,maxFakeErrorX);
		s.pos.y=s.pos.y+randDouble(-maxFakeErrorY,maxFakeErrorY);
		if (verbose_level>0) std::cout << "added random (-maxFakeError,maxFakeError) error to new_scan position";
	}
	else if (fakeError.compare("range") == 0)
	{	s.pos.rot=angleUnwrap(s.pos.rot+er);
		s.pos.x=s.pos.x+ex;
		s.pos.y=s.pos.y+ey;
		if (verbose_level>0) std::cout << "added fixed ("<<ex<<","<<ey<<","<<er<<") error to new_scan position";
		ex+=estep;		// each call changes error allowing to sweep entire range (for same scan typically)
		if (ex>maxFakeErrorX)
		{	ey+=estep; ex=-maxFakeErrorX;
			if (ey>maxFakeErrorY)
			{	er+=estep; ey=-maxFakeErrorY;
			}
		}
	}
	else if (fakeError.compare("fixed_all_combinations") == 0)
	{	// used to test one scan with all 8 combinations of fixed error (+-ex,+-ey,+-er)
		s.pos.rot=angleUnwrap(s.pos.rot+er);
		s.pos.x=s.pos.x+ex;
		s.pos.y=s.pos.y+ey;
		if (verbose_level>0) std::cout << "added fixed_all_combinations ("<<ex<<","<<ey<<","<<er<<") error to new_scan position";
		ex=-ex;	// change ex each time
		if (ex<0)
		{	ey=-ey; // change ey every second time
			if (ey<0)
			{	er=-er; // change er every fourth time
			}
		}
	}
	else
	{	// new_scan is fixed
		s.pos.rot=angleUnwrap(s.pos.rot+maxFakeErrorRot);
		s.pos.x=s.pos.x+maxFakeErrorX;
		s.pos.y=s.pos.y+maxFakeErrorY;
		if (verbose_level>0) std::cout << "added fixed ("<<maxFakeErrorX<<","<<maxFakeErrorY<<","<<maxFakeErrorRot<<") error to new_scan position";
	}
	std::cout << std::endl;

	// both changing scan ref frame (ie. scan position) or manipulating readings can change
	// the order of scans which are expected ordered by angle, so we must restore it
	s.readings.sort(readingsAngleCompare);
} // add_fake_error_pos


void add_fake_error_readings(scan &s)
{
	/* add fake error - initial error that we add to new scan to test robustness
	 * The error, depending on chosen method, may span between [-maxFakeError,maxFakeError]
	 * There are two ways to do add the error:
	 * 1) change scan.pos
	 * 2) change every reading keeping scan.pos intact
	 * This procedure uses the second way.
	 */
	if (verbose_level>0) std::cout << "[add_fake_error_readings   ] ";
	if (fakeError.compare("random") == 0)
	{   er=randDouble(-maxFakeErrorRot,maxFakeErrorRot);
		ex=randDouble(-maxFakeErrorX,maxFakeErrorX);
		ey=randDouble(-maxFakeErrorY,maxFakeErrorY);
		readingsChangeRefFrame(s.readings, {0,0,0}, {ex,ey,er});
		if (verbose_level>0) std::cout << "added random error to readings of scan";
	}
	else if (fakeError.compare("range") == 0)
	{	readingsChangeRefFrame(s.readings, {0,0,0}, {ex,ey,er});
		if (verbose_level>0) std::cout << "added fixed ("<<ex<<","<<ey<<","<<er<<") error to readings of scan";
		ex+=estep;	// each call changes error allowing to sweep entire range (for same scan typically)
		if (ex>maxFakeErrorX)
		{	ey+=estep; ex=-maxFakeErrorX;
			if (ey>maxFakeErrorY)
			{	er+=estep; ey=-maxFakeErrorY;
			}
		}
	}
	else if (fakeError.compare("fixed_all_combinations") == 0)
	{	// used to test one scan with all 8 combinations of fixed error (+-ex,+-ey,+-er)
		ex=-ex;	// change ex each time
		if (ex<0)
		{	ey=-ey; // change ey every second time
			if (ey<0)
			{	er=-er; // change er every fourth time
			}
		}
		readingsChangeRefFrame(s.readings, {0,0,0}, {ex,ey,er});
		if (verbose_level>0) std::cout << "added fixed_all_combinations ("<<ex<<","<<ey<<","<<er<<") error to readings of scan";
	}
	else
	{   // new_scan is fixed
		readingsChangeRefFrame(s.readings, {0,0,0}, {maxFakeErrorX,maxFakeErrorY,maxFakeErrorRot});
		if (verbose_level>0) std::cout << "added fixed ("<<maxFakeErrorX<<","<<maxFakeErrorY<<","<<maxFakeErrorRot<<") error to readings of scan";
	}
	if (verbose_level>0) std::cout << std::endl;

	// both changing scan ref frame (ie. scan position) or manipulating readings can change
	// the order of scans which are expected ordered by angle, so we must restore it
	s.readings.sort(readingsAngleCompare);
} // add_fake_error_readings



void pre_process_scans()
{
	// optional (but recommended): resample and filter scans
	if (median_filter_scans) scanMedianFilter(ws1, median_filter_order, median_filter_distance);
	if (resample_scans)      scanResample(ws1, resample_distance);
	if (median_filter_scans) scanMedianFilter(ws2, median_filter_order, median_filter_distance);
	if (resample_scans)      scanResample(ws2, resample_distance);

	if (median_filter_scans)
	{ 	if (verbose_level>0)
		{	std::cout << "[pre_process_scans         ] median filter of order "
			<<  median_filter_order << " and distance " << median_filter_distance
			<< " applied to sref and snew" << std::endl;
		}
	}

	if (resample_scans)
	{	if (verbose_level>0)
		{	std::cout << "[pre_process_scans         ] resampling with resample distance "
			<< resample_distance << " applied to sref and snew" << std::endl;
		}
	}


	// optional (but recommended): filter out (for 2nd time) readings out of range
	SIT LI=ws1.readings.begin();
	while (LI!=ws1.readings.end())
	{   if (LI->distance>maxRange()) LI=ws1.readings.erase(LI); else LI++;
	}
	LI=ws2.readings.begin();
	while (LI!=ws2.readings.end())
	{   if (LI->distance>maxRange()) LI=ws2.readings.erase(LI); else LI++;
	}

	if (resample_scans)
		if (verbose_level>0) std::cout << "[pre_process_scans         ] filtering out readings out of range" << std::endl;

} // pre_process_scans


bool is_failure(double dist, double rot, unsigned iter)
{	if (rot>M_PI) rot-=2*M_PI;
	rot=fabs(rot);
	//std::cout << "isfailure dist=" << dist << " rot = " << rot << std::endl;
	if (isnan(dist) || isnan(rot)
	|| dist>success_treshold_traslation
	|| rot>success_treshold_rotation
	|| dist<0.0 )
	{
		return true;
	}
	return false;
}


int main(int argc, char **argv)
{	srand(time(NULL));
	srand(10);

	read_simulation_parameters();

	if (NumAlgorithms < 1)
	{   std::cout<<"Could not read any algorithm in simulation parameters file" << std::endl;
		return -1;
	}

	// set initial fake error
	ex=maxFakeErrorX;	// it gets inverted when first add_fake_error is called in **_all_combinations variant
	ey=maxFakeErrorY;
	er=maxFakeErrorRot;

//// temp
//// use this to save a sample of scans that rawseeds_scans vector is holding
//load_scan_pair(0);
//rawseeds_write_random_sample("/home/klenac/scans/Bovisa_2008-09-01-SICK_FRONT_random_sample_1.csv",1000);
//exit(1);
//// end temp


	// initializing variables for calculation of average and variance
	const unsigned	NumPairs=loop_end-loop_start+1;

	double	avg_time[NumAlgorithms][NumPairs];
	double	errx[NumAlgorithms][NumPairs];		// store added errors for later processing
	double	erry[NumAlgorithms][NumPairs];
	double	errdist[NumAlgorithms][NumPairs];
	double	errrot[NumAlgorithms][NumPairs];
	double	x[NumAlgorithms][NumPairs];
	double	y[NumAlgorithms][NumPairs];
	double	dist[NumAlgorithms][NumPairs];
	double	rot[NumAlgorithms][NumPairs];
	double	iter[NumAlgorithms][NumPairs];
	int		status[NumAlgorithms][NumPairs];

	unsigned numiter;

	for (unsigned i=0; i<NumAlgorithms; i++)
	{   for (unsigned j=0; j<NumPairs; j++)
		{   avg_time[i][j]=0.0;
			errx[i][j]=0.0;
			erry[i][j]=0.0;
			errdist[i][j]=0.0;
			errrot[i][j]=0.0;
			dist[i][j]=0.0;
			rot[i][j]=0.0;
			x[i][j]=0.0;
			y[i][j]=0.0;
			status[i][j]=0;
		}
	}

	unsigned NumValidPairs=0;
	//for (unsigned pair=loop_start; pair<=loop_end; pair++)
	// carefull: we are using NumValidPairs instead of pair to exit from loop!
	for (unsigned pair=loop_start; NumValidPairs<=loop_end; pair++)
	{   /* main testing loop:
		 * 1) we first read one scan pair from disk/file and store it ts1 (ref scan) and ts2 (new scan)
		 * The scan.pos typically contains the true position (ground truth)
		 * ts1 and ts2 are never modified by the program
		 * 2) ts1 and ts2 are then copied in ws1 and ws2 (work scans) which can be modified
		 * 3) a fake error is then added to ws2
		 * 4) ws1 and ws2 are optionally preprocessed
		 * 5) ws1 and ws2 are copied in sws1 and sws2 so they can be restored for each test
		 * 6) algorithm is called passing ws1 and ws2 as reference
		 */
		//std::cout << pair << std::endl;
		load_scan_pair(pair);                       // load one scan pair according to preferences in sm.ini

		// don't care about stored scan position. We are interested only in readings (well, in these tests at least...)
		ts1.pos.rot=ts1.pos.x=ts1.pos.y=0.0;
		ts2.pos.rot=ts2.pos.x=ts2.pos.y=0.0;

		// working scans (leave ts1 and ts2 intact)
		ws1.copy(ts1);
		ws2.copy(ts2);

		// add fake error according to preferences in sm.ini
		add_fake_error_readings(ws2);

		// pre_process scans according to preferences in sm.ini
		pre_process_scans();

		// skip scans with less than min_scan_readings readings
		if (ws1.readings.size()<min_scan_readings || ws2.readings.size()<min_scan_readings)
		{   if (verbose_level>0) std::cout << "[main loop                 ] this scan pair has too few readings, skipping ..." << std::endl;
			continue;
		}

		// store the working scans
		sws1.copy(ws1);
		sws2.copy(ws2);

		// loop through algorithms defined in sm.ini
		for(unsigned alg=0; alg<NumAlgorithms; alg++)
		{	// restore working scans - all algorithms use the same scans
			ws1.copy(sws1);
			ws2.copy(sws2);

			// store the fake error added to scan (in order to plot number of iterations against error)
			errx[alg][NumValidPairs]=ex;
			erry[alg][NumValidPairs]=ey;
			errdist[alg][NumValidPairs]=sqrt(ex*ex+ey*ey);
			errrot[alg][NumValidPairs]=er;

			struct timeval	tm_init,
					tm_pre_match,
					tm_match,
					tm_post_match,
					tm_deinit,
					tm_end;

			gettimeofday(&tm_init, NULL);

			// read algorithm parameters and initialize algorithm structures - one time initialization
			// todo: take init/deinit out of the pairs loop (now it gets executed for every pair)
			algorithm[alg].init(algorithm[alg].parameters_file.c_str());

			double avgruntime=0;
			for (unsigned j=0; j<runs;j++) // these runs are just to average chronometer
			{   long oneruntime, seconds, useconds;
				numiter=0;
				gettimeofday(&tm_pre_match, NULL);
				algorithm[alg].pre_match(ws1,ws2);	// prepare to match reference with new scan
				gettimeofday(&tm_match, NULL);
				mres=algorithm[alg].match();		// matching!
				gettimeofday(&tm_post_match, NULL);
				algorithm[alg].post_match();		// maybe there is something to do after matching

				seconds  = tm_post_match.tv_sec  - tm_match.tv_sec;
				useconds = tm_post_match.tv_usec - tm_match.tv_usec;
				oneruntime = (seconds) * 1000000 + useconds;
				if (j>=skip_first_runs) avgruntime+=double(oneruntime);
			}
			gettimeofday(&tm_deinit, NULL);
			algorithm[alg].deinit();				// deinitialize algorithm and structures
			gettimeofday(&tm_end, NULL);
			avgruntime/=double(runs-skip_first_runs);

			if (verbose_level>0)
			{
				// print resulting position and computation time
				std::cout << "Algorithm "<< algorithm[alg].name << std::endl;
				std::cout << "\tEstimated position: \t\t\t(" << mres.p.x << "," << mres.p.y << "," << mres.p.rot << ")" << std::endl;
				std::cout << "\tComputation time (microseconds): \t" << avgruntime << std::endl;
			}
			// store the distance and rotation error and average time for later calculation of averages and variances
			avg_time[alg][NumValidPairs] = avgruntime;
			status[alg][NumValidPairs]=mres.status;                        // for now only mbicp uses this (analysis of true and false positives)

			if (!isnan(mres.p.x) && !isnan(mres.p.y) && !isnan(mres.p.rot))
			{

#ifdef DRAW_PNG
				if (draw_scans)
				{	const unsigned SIZE = 255;
					char filename[SIZE];
					snprintf(filename,SIZE,"./images/alg%d_%s.png",pair,algorithm[alg].name.c_str());
					drawing img(filename, map_size_x, map_size_y,-map_size_x/2.0,-map_size_y/2.0,map_filename.c_str());
					img.grid();
					img.setpointsize(6);
					img.setcolor(0,0,30000);
					img.readings(ts1,ts1.pos);
					img.setpointsize(4);
					img.setcolor(0,0,50000);
					img.readings(ts2,ts2.pos);
					rs.copy(ws2);
					rs.pos=mres.p;                // rs.pos contains corrected position (ie motion estimation)
					img.setpointsize(4);
					img.setcolor(0,30000,0);
					img.readings(ws2,ws2.pos);
					img.setcolor(0,50000,0);
					img.setpointsize(2);
					img.readings(rs,rs.pos);
					img.write();
				}
#endif

				x[alg][NumValidPairs]=mres.p.x;
				y[alg][NumValidPairs]=mres.p.y;
				dist[alg][NumValidPairs]=sqrt((mres.p.x-ex)*(mres.p.x-ex)+(mres.p.y-ey)*(mres.p.y-ey));
				rot[alg][NumValidPairs]=angleUnwrap(mres.p.rot-er);
				if (verbose_level>0) std::cout << "\tDistance from true pos (dx,dy,drot): \t("<<mres.p.x-ex<<","<<mres.p.y-ey<<","<<rot[alg][NumValidPairs]<<")  dist:"<<dist[alg][NumValidPairs]<<std::endl;
			}
			else
			{   dist[alg][NumValidPairs]=-1.0;    // we use negative distance to indicate invalid matching. Used later when calculating averages and variances
				rot[alg][NumValidPairs]=0.0;
			}
			iter[alg][NumValidPairs]=numiter;


		} // end algorithms loop
		NumValidPairs++;
	} // end pairs loop

	if (NumValidPairs==0) return 0;

	double SuccessRatio[NumAlgorithms];

	if (verbose_level>0) std::cout << "---" <<std::endl;

#ifdef ERRvsITER
	// write number of iterations against IPE in file for later drawing/analysis
	const char *filename = "errvsiter.csv";
	FILE *myff;
	if ((myff = fopen(filename, "w")) == NULL)
	{   std::cout<<"write iter: Error opening file"<< std::endl;
	}
#endif

	for (unsigned alg=0; alg<NumAlgorithms; alg++)
	{   algorithm[alg].mtime_success=0.0;
		algorithm[alg].vtime_success=0.0;
		algorithm[alg].mtime_failure=0.0;
		algorithm[alg].vtime_failure=0.0;
		algorithm[alg].mdist=0.0;
		algorithm[alg].mrot=0.0;
		algorithm[alg].vdist=0.0;
		algorithm[alg].vrot=0.0;
		algorithm[alg].miter_success=0.0;
		algorithm[alg].viter_success=0.0;
		algorithm[alg].miter_failure=0.0;
		algorithm[alg].viter_failure=0.0;
		SuccessRatio[alg]=0.0;
		algorithm[alg].successes=0;
		algorithm[alg].failures=0;
		algorithm[alg].false_positives=0;
		algorithm[alg].true_negatives=0;
		algorithm[alg].false_negatives=0;
		algorithm[alg].true_positives=0;

		std::cout << std::endl << std::endl << "Algorithm: " << algorithm[alg].name << std::endl;
		for (unsigned pair=0; pair<NumValidPairs; pair++)
		{   //std::cout<<"dist["<<alg<<"]["<<pair<<"]: "<<dist[alg][pair]<<"  rot["<<alg<<"]["<<pair<<"]: "<<rot[alg][pair];

			if (is_failure(dist[alg][pair],rot[alg][pair],iter[alg][pair]))
			{	algorithm[alg].failures++;
				if (status[alg][pair]==1) algorithm[alg].false_positives++; else algorithm[alg].true_negatives++;
				algorithm[alg].mtime_failure+=avg_time[alg][pair];
				algorithm[alg].miter_failure+=iter[alg][pair];
				//std::cout<<" -> failure dist"<<dist[alg][pair]<<" rot"<<rot[alg][pair]<<"errdist"<<errdist[alg][pair]<<" errrot"<<errrot[alg][pair]<<std::endl;
				//std::cout<<" -> failure"<<std::endl;

			} else
			{   algorithm[alg].successes++;
				if (status[alg][pair]==2) algorithm[alg].false_negatives++; else algorithm[alg].true_positives++;
				algorithm[alg].mdist+=dist[alg][pair];
				algorithm[alg].mrot+=fabs(rot[alg][pair]);
				algorithm[alg].mtime_success+=avg_time[alg][pair];
				algorithm[alg].miter_success+=iter[alg][pair];
				//std::cout<<" -> success"<<std::endl;
			}
		}
		if (verbose_level>0)
		{
			std::cout<<"successes("<<algorithm[alg].successes<<")  failures("<<algorithm[alg].failures<<")"<<std::endl;
			std::cout<<"false_positives("<<algorithm[alg].false_positives<<")  true_positives("<<algorithm[alg].true_positives<<")"<<std::endl;
			std::cout<<"false_negatives("<<algorithm[alg].false_negatives<<")  true_negatives("<<algorithm[alg].true_negatives<<")"<<std::endl;
		}
#ifdef ERRvsITER
		if (fakeError.compare("range") == 0) {
			std::cout << std::endl << std::endl << "Sensitivity to IPE error analysis: " << std::endl;
			std::cout << std::endl << std::endl << "edist,erot,1=success/0=failure" << std::endl;
			std::cout.precision(4);
			for (unsigned pair=0; pair<NumValidPairs; pair++)
			{
				//std::cout<<pair<<", "<<errdist[alg][pair]<<", "<<errrot[alg][pair]<<", "<<errx[alg][pair]<<", "<<erry[alg][pair]<<", ";
				std::cout<<errdist[alg][pair]<<", "<<errrot[alg][pair]<<", ";
				if (is_failure(dist[alg][pair],rot[alg][pair],iter[alg][pair]))
				{	std::cout<<"0"<<std::endl;

				} else
				{   std::cout<<"1"<<std::endl;
				}
			}
		}
#endif

		if (algorithm[alg].failures>0)
		{	algorithm[alg].mtime_failure/=algorithm[alg].failures;
			algorithm[alg].miter_failure/=algorithm[alg].failures;
		}
		if (algorithm[alg].successes>0)
		{	algorithm[alg].mtime_success/=algorithm[alg].successes;
			algorithm[alg].miter_success/=algorithm[alg].successes;
			algorithm[alg].mdist/=algorithm[alg].successes;
			algorithm[alg].mrot/=algorithm[alg].successes;
		}
		SuccessRatio[alg]=algorithm[alg].successes/double(NumValidPairs); // NumValidPairs is guaranteed to be > 0, otherwise we wouldn't be in here

		// variances
		for(unsigned pair=0;pair<NumValidPairs;pair++)
		{	if (is_failure(dist[alg][pair],rot[alg][pair],iter[alg][pair]))
			{   // skip distance and rotation variance calculation (consider variances only for successful matchings)
				algorithm[alg].vtime_failure+=(avg_time[alg][pair]-algorithm[alg].mtime_failure)*(avg_time[alg][pair]-algorithm[alg].mtime_failure);
				algorithm[alg].viter_failure+=(iter[alg][pair]-algorithm[alg].miter_failure)*(iter[alg][pair]-algorithm[alg].miter_failure);
			}
			else
			{   algorithm[alg].vdist+=(dist[alg][pair]-algorithm[alg].mdist)*(dist[alg][pair]-algorithm[alg].mdist);
				algorithm[alg].vrot+=(rot[alg][pair]-algorithm[alg].mrot)*(rot[alg][pair]-algorithm[alg].mrot);
				algorithm[alg].vtime_success+=(avg_time[alg][pair]-algorithm[alg].mtime_success)*(avg_time[alg][pair]-algorithm[alg].mtime_success);
				algorithm[alg].viter_success+=(iter[alg][pair]-algorithm[alg].miter_success)*(iter[alg][pair]-algorithm[alg].miter_success);
#ifdef ERRvsITER
				if (fprintf(myff,"%f,%f,%f\n",errdist[alg][pair],errrot[alg][pair],iter[alg][pair])<3)
				{	std::cout<<"write iter: Error writing iter"<< std::endl;
				}
#endif
			}
		}

		if (algorithm[alg].failures>0)
		{	algorithm[alg].vtime_failure/=algorithm[alg].failures;
			algorithm[alg].viter_failure/=algorithm[alg].failures;
		}
		if (algorithm[alg].successes>0)
		{	algorithm[alg].vtime_success/=algorithm[alg].successes;
			algorithm[alg].viter_success/=algorithm[alg].successes;
			algorithm[alg].vdist/=algorithm[alg].successes;
			algorithm[alg].vrot/=algorithm[alg].successes;
		}

		if (verbose_level>0)
		{
			// print resulting position and computation time
			std::cout << std::endl;
			std::cout << "Valid pairs: " << NumValidPairs <<
				" SuccessRatio(" << SuccessRatio[alg] << ")" << std::endl;
			std::cout << "computation time (microseconds) for succesfull matchings" <<
				" average: " << algorithm[alg].mtime_success <<
				" variance: " << algorithm[alg].vtime_success <<
				" std dev: " << sqrt(algorithm[alg].vtime_success) << std::endl;
			std::cout << "number of iterations for succesfull matchings" <<
				" average: " << algorithm[alg].miter_success <<
				" variance: " << algorithm[alg].viter_success <<
				" std dev: " << sqrt(algorithm[alg].viter_success) << std::endl;
			std::cout << "computation time (microseconds) for failed matchings" <<
				" average: " << algorithm[alg].mtime_failure <<
				" variance: " << algorithm[alg].vtime_failure <<
				" std dev: " << sqrt(algorithm[alg].vtime_failure) << std::endl;
			std::cout << "number of iterations for failed matchings" <<
				" average: " << algorithm[alg].miter_failure <<
				" variance: " << algorithm[alg].viter_failure <<
				" std dev: " << sqrt(algorithm[alg].viter_failure) << std::endl;
			std::cout << "accuracy (meters)" <<
				" average(" << algorithm[alg].mdist << "," << algorithm[alg].mrot << ")" <<
				" variance(" << algorithm[alg].vdist << "," << algorithm[alg].vrot << ")" <<
				" std dev(" << sqrt(algorithm[alg].vdist) << "," << sqrt(algorithm[alg].vrot) << ")" << std::endl;
		}

		std::cout << std::endl;
		std::cout << std::fixed;
		std::cout.precision(1);
		std::cout << SuccessRatio[alg]*100 << " %" << std::endl;
		std::cout.precision(2);
		std::cout << algorithm[alg].mtime_success/1000.0 << " ("
				<< sqrt(algorithm[alg].vtime_success)/1000.0 << ") ms" << std::endl;
		std::cout << algorithm[alg].mdist*100 << " (" << sqrt(algorithm[alg].vdist)*100 << ") mm" << std::endl;
		std::cout << algorithm[alg].mrot << " (" << sqrt(algorithm[alg].vrot) << ") mrad" << std::endl;
	} // end algorithms loop

#ifdef ERRvsITER
	if (fclose(myff) !=  0)
	{	std::cout<<"write iter: Error closing file"<< std::endl;
	}
#endif

	return 1;
} // main()
// -----------------------------------------------------------------------------
