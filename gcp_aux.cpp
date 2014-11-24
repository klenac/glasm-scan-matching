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


// gcp_aux.cpp functions help interfacing to sm_cli simulator


#include <iostream>
#include <fstream>
#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "gcp_aux.h"

#ifdef DRAW_PNG
        #include "mydraw.h"
#endif

namespace GCP
{

double          outfitness;

/*
void read_parameters(const char* filename)
{
    int     nbitx;
    int     nbity;
    int     nbitrot;
    int     pcross;
    double  pmutation;
    int     popsize;
    int     maxruns;
    int     maxgen;
    int     RisulArea;
    double  SearchAreaDist;
    double  SearchAreaRot;

    double  MinX;
    double  MaxX;
    double  MinY;
    double  MaxY;
    double  MinFI;
    double  MaxFI;

    unsigned        lookup_rows;        // number of rows= 3000;
    unsigned        lookup_columns;     // number of columns= 3000;
    double          lookup_size_x;      // dimensione lato x area in metri =10.0
    double          lookup_size_y;      // dimensione lato y area in metri =10.0
    double          corr_distance;      // distanza entro la quale cercare i corrispondenti
    double          lookup_ox;          // coord x dell'angolo inf sin dell'area a ScanRef.pos.x-lookup_size_x/2.0
    double          lookup_oy;          // ScanRef.pos.y-lookup_size_y/2.0



    try {

        std::ifstream config_stream(filename);
        namespace po = boost::program_options;
        po::options_description my_options("Options");

        // todo: option global_search                   = true

        my_options.add_options()
            // search area (specify either exact range or by distance and rotation. the former takes precedence if both specified)
           ("MinX", po::value(&MinX)->default_value(0), "SearchArea MinX")
           ("MaxX", po::value(&MaxX)->default_value(1), "SearchArea MaxX")
           ("MinY", po::value(&MinY)->default_value(0), "SearchArea MinY")
           ("MaxY", po::value(&MaxY)->default_value(1), "SearchArea MaxY")
           ("MinFI", po::value(&MinFI)->default_value(0), "SearchArea MinFI")
           ("MaxFI", po::value(&MaxFI)->default_value(1), "SearchArea MaxFI")

           ("SearchAreaDist", po::value(&SearchAreaDist)->default_value(0), "SearchArea Distance")
           ("SearchAreaRot", po::value(&SearchAreaRot)->default_value(1), "SearchArea Rotation")

           // genetic parameters
           ("nbitx", po::value(&nbitx)->default_value(8), "nbitx")
           ("nbity", po::value(&nbity)->default_value(8), "nbity")
           ("nbitrot", po::value(&nbitrot)->default_value(8), "nbitrot")
           ("pcross", po::value(&pcross)->default_value(8), "pcross")
           ("pmutation", po::value(&pmutation)->default_value(0.012), "pmutation")
           ("popsize", po::value(&popsize)->default_value(100), "popsize")
           ("maxruns", po::value(&maxruns)->default_value(1), "maxruns")
           ("maxgen", po::value(&maxgen)->default_value(10), "maxgen")

           // lookup table parameters
           ("lookup_rows", po::value(&lookup_rows)->default_value(3000), "lookup_rows")
           ("lookup_columns", po::value(&lookup_columns)->default_value(3000), "lookup_columns")
           ("lookup_size_x", po::value(&lookup_size_x)->default_value(40.0), "lookup_size_x")
           ("lookup_size_y", po::value(&lookup_size_y)->default_value(40.0), "lookup_size_y")
           ("corr_distance", po::value(&corr_distance)->default_value(0.20), "corr_distance")
           ("lookup_ox", po::value(&lookup_ox)->default_value(-20.0), "lookup_ox")
           ("lookup_oy", po::value(&lookup_oy)->default_value(-20.0), "lookup_oy")
           ;

        po::variables_map vm;
        po::store(po::parse_config_file(config_stream, my_options), vm);
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
    initializeInvGrayTable();
    setGeneticParameters(
            nbitx,
            nbity,
            nbitrot,
            popsize,
            maxruns,
            maxgen,
            pcross,
                pmutation);
    setLookupTableParameters(RisulArea,RisulArea,lookup_size_x,lookup_size_y,0,0,0.20); // only params, the lookup table will be recalculated each step
    initializeRadialGradient(RGmean,RGvariance,RGextension,RGweight,RGscale,RGbitmapsize); // RGbitmapsize is equivalent to GLASM_CORR_DISTANCE_TRESHOLD in the old GLASM
    initializeLookupTableStaticBitmap("./cfg/bitmaps/cave_highres.png",16.0,16.0,2.0);
    // todo: set numerical parameters in ini file


} // read_parameters


void init()
{
}

void pre_match(const scan &snew, const scan &sref)
{
    setRefScan(sref);
    setNewScan(snew);
    initializeLookupTable();
    setCenterOfSearchArea(snew.pos);     // todo: if global search center in sref.pos
}


position match()
{
    return gcp(&outfitness);
}

void post_match()
{
}

void deinit()
{
    deleteLookupTable();
}




*/

} // namespace GCP
