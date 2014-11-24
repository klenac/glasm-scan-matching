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


// mbicp_aux.cpp functions help interfacing to sm_cli simulator


#include <iostream>
#include <fstream>
#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <math.h>

#include "mbicp_aux.h"

#ifdef DRAW_PNG
	#include <pngwriter.h>
	#include "mydraw.h"
#endif

#include "mbicp/MbICP.h"

namespace MBICP
{

// mbicp
Tpfp            laserK[MAXLASERPOINTS];        // used by mbicp
Tpfp            laserK1[MAXLASERPOINTS];
Tsc             sensorMotion;
Tsc             solution;


bool			draw_iterations;
std::string 	map_filename;
double			map_size_x;
double 			map_size_y;

float           max_laser_range;
float           Bw;
float           Br;
float           L;
int             laserStep;
int             ProjectionFilter;
float           MaxDistInter;
float           filter;
float           AsocError;
int             MaxIter;
float           errorRatio;
float           errx_out;
float           erry_out;
float           errt_out;
int             IterSmoothConv;



void scan2mbicpscan(const scan &s, Tpfp *lK)
{	for (int i=0; i<MAXLASERPOINTS; i++) {
				lK[i].r=max_laser_range;
				lK[i].t=2.0*M_PI/MAXLASERPOINTS;
	}

	CSIT si=s.readings.begin();
	int i=0;
	while (si!=s.readings.end())
	{	i=si->angle/(2.0*M_PI/MAXLASERPOINTS);
		lK[i].r=si->distance;
		lK[i].t=si->angle;
		si++; i++;
	}
}


void init(const char* filename)
{



	try {

	std::ifstream config_stream(filename);
	namespace po = boost::program_options;
	po::options_description my_options("Options");

	my_options.add_options()
#ifdef DRAW_PNG
	("draw_iterations", po::value<bool>(&draw_iterations)->default_value("false"), "set true to draw the mbicp iterations")
	("map_filename", po::value<std::string>(&map_filename)->default_value("none"), "map filename")
	("map_size_x", po::value(&map_size_x)->default_value(15.0), "map_size_x")
	("map_size_y", po::value(&map_size_y)->default_value(15.0), "map_size_y")
#endif
		("L", po::value(&L)->default_value(3.0), "L")
		("max_laser_range", po::value(&max_laser_range)->default_value(7.9), "max_laser_range")
		("Bw", po::value(&Bw)->default_value(0.523333333), "Bw")
		("Br", po::value(&Br)->default_value(0.3), "Br")
		("laserStep", po::value(&laserStep)->default_value(1), "laserStep")
		("MaxDistInter", po::value(&MaxDistInter)->default_value(0.5), "MaxDistInter")
		("filter", po::value(&filter)->default_value(0.95), "filter")
		("ProjectionFilter", po::value(&ProjectionFilter)->default_value(0), "ProjectionFilter")
		("AsocError", po::value(&AsocError)->default_value(0.1), "AsocError")
		("MaxIter", po::value(&MaxIter)->default_value(50), "MaxIter")
		("errorRatio", po::value(&errorRatio)->default_value(0.0001), "errorRatio")
		("errx_out", po::value(&errx_out)->default_value(0.0001), "errx_out")
		("erry_out", po::value(&erry_out)->default_value(0.0001), "erry_out")
		("errt_out", po::value(&errt_out)->default_value(0.0001), "errt_out")
		("IterSmoothConv", po::value(&IterSmoothConv)->default_value(2), "IterSmoothConv")
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

// set parameters
Init_MbICP_ScanMatching(
			max_laser_range,
			Bw,
			Br,
			L,
			laserStep,
			MaxDistInter,
			filter,
			ProjectionFilter,
			AsocError,
			MaxIter,
			errorRatio,
			errx_out,
			erry_out,
			errt_out,
			IterSmoothConv,
			draw_iterations);
} // init




#ifdef DRAW_PNG
	// drawing the matching process
	void draw_mbicp_associations(int numIteration)
	{   const int FILENAME_SIZE = 19;
		char filename[FILENAME_SIZE]="images/mite_50.tmp";
		snprintf(filename,FILENAME_SIZE,"images/mite_%02d.png",numIteration);
		drawing img(filename,map_size_x,map_size_y,-map_size_x/2.0,-map_size_y/2.0,map_filename.c_str());
		img.grid();
		FILE *in;
		float x1,y1,x2,y2;

		// draw associations
		snprintf(filename,FILENAME_SIZE,"images/mite_%02d.tmp",numIteration);
		if ((in = fopen(filename, "r")) == NULL)
		{   std::cout<<"read_mbicp_associations: Error opening file "<<filename<< std::endl;
			return;
		}
		while (fscanf(in,"%f,%f,%f,%f",&x1,&y1,&x2,&y2)==4)
		{   img.line(x1,y1,x2,y2);
		}
		if (fclose(in) !=  0)
		{   std::cout<<"read_mbicp_associations: Error closing file"<< std::endl;
		}
		int status;
		if( (status = remove(filename)) == 0 )
		{	// printf("%s file deleted successfully.\n",filename);
		}
		else
		{	printf("Unable to delete the file\n");
			perror("Error");
		}

		/* todo: to draw scans of mbicp during iterations we should first export ptosNewRef in EStep() in MbICP.c in every iteration
		 *

		// draw reference scan
		snprintf(filename,FILENAME_SIZE,"images/msc1_%d.tmp",numIteration);
		if ((in = fopen(filename, "r")) == NULL)
		{   std::cout<<"read_mbicp_scan1: Error openinlinuxg file"<< std::endl;
		}
		img.setpointsize(4);
		img.setcolor(30000,0,0);
		while (fscanf(in,"%f,%f",&x1,&y1)==2)
		{   img.point(x1,y1);
			std::cout<<"drawingpointscan1: "<< x1 << " " << y1 << std::endl;
		}
		if (fclose(in) !=  0)
		{   std::cout<<"read_mbicp_scan1: Error closing file"<< std::endl;
		}
		if( (status = remove(filename)) == 0 )
		{
			// printf("%s file deleted successfully.\n",filename);
		}
		else
		{
			printf("Unable to delete the file\n");
			perror("Error");
		}

		// draw new scan
		snprintf(filename,FILENAME_SIZE,"images/msc2_%d.tmp",numIteration);
		if ((in = fopen(filename, "r")) == NULL)
		{   std::cout<<"read_mbicp_scan2: Error opening file"<< std::endl;
		}
		img.setcolor(0,30000,0);
		img.setpointsize(2);
		while (fscanf(in,"%f,%f",&x1,&y1)==2)
		{   img.point(x1,y1);
			std::cout<<"drawingpointscan2: "<< x1 << " " << y1 << std::endl;
		}
		if (fclose(in) !=  0)
		{   std::cout<<"read_mbicp_scan2: Error closing file"<< std::endl;
		}
		if( (status = remove(filename)) == 0 )
		{	// printf("%s file deleted successfully.\n",filename);
		}
		else
		{	printf("Unable to delete the file\n");
			perror("Error");
		}
		*/

		// close image
		img.write();
	}
#endif


void pre_match(const scan &sref, const scan &snew)
{	/*  It seems that MBICP does not work if both scans are
	 *  not in the same reference frame. Even so, for some reference frames (origins) mbicp
	 *  also failes to converge giving NaN,NaN,NaN as result, while just changing the ref frame
	 *  for input scans (i.e. same points in global ref frame) makes it work again. */

	// convert scans in arrays used by mbicp algorithm
	scan2mbicpscan(sref,MBICP::laserK);
	scan2mbicpscan(snew,MBICP::laserK1);

	sensorMotion.tita=angleUnwrap(snew.pos.rot-sref.pos.rot);
	sensorMotion.x=snew.pos.x-sref.pos.x;
	sensorMotion.y=snew.pos.y-sref.pos.y;

/*
	// temp just for debug - draw image with scans input to mbicp algorithm
	drawing img1("images/blabla.png", map_size_x, map_size_y,-map_size_x/2.0,-map_size_y/2.0,map_filename.c_str());
	img1.grid();
	img1.setpointsize(6);
	img1.setcolor(0,0,30000);
	img1.readings(sref,sref.pos);
	img1.setpointsize(4);
	img1.setcolor(0,0,50000);
	img1.readings(snew,snew.pos);
	img1.setpointsize(2);
	img1.setcolor(0,30000,0);
	for(unsigned i=0; i<MAXLASERPOINTS;i++)
	{	if (laserK[i].r<maxRange())
		{	img1.point(laserK[i].r*cos(laserK[i].t),laserK[i].r*sin(laserK[i].t));
		}
	}
	img1.setpointsize(1);
	img1.setcolor(0,50000,0);
	for(unsigned i=0; i<MAXLASERPOINTS;i++)
	{	if (laserK1[i].r<maxRange())
		{	img1.point(laserK1[i].r*cos(laserK1[i].t),laserK1[i].r*sin(laserK1[i].t));
		}
	}
	img1.write();
	// end temp
*/
}


matching_result match()
{   matching_result mres;
	mres.status=MbICPmatcher(laserK, laserK1, &sensorMotion, &solution, &mres.numiter);

#ifdef DRAW_PNG
	if (draw_iterations)
	{	for(int i=0; i<mres.numiter+1; i++) draw_mbicp_associations(i);
	}

#endif

	if (mres.status==1)
	{   // converged
		mres.p.x=solution.x;
		mres.p.y=solution.y;
		mres.p.rot=solution.tita;
	} else if (mres.status==2)
	{   // did not converge (reached max number of iterations)
		mres.p.x=solution.x;
		mres.p.y=solution.y;
		mres.p.rot=solution.tita;
	} else if (mres.status==-1)
	{   // did not converge (too few correspondences at one point)
		mres.p.x=NAN;
		mres.p.y=NAN;
		mres.p.rot=NAN;
	} else if (mres.status==-2)
	{   // did not converge (could not perform calculation)
		mres.p.x=NAN;
		mres.p.y=NAN;
		mres.p.rot=NAN;
	}
	return mres;
}

void post_match()
{
}

void deinit()
{
}

} // namespace
