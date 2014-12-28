glasm-scan-matching
===================

Implementation of some laser scan matching algorithms and a program to compare and test them.


Contents
--------

1. Introduction
2. Contents of the repository
3. Setting up a development environment with Eclipse
4. Building the sm_cli program
5. Preparing for the experiments
6. Running the tests
7. Analysing the results
8. Developing



## Introduction ##

This repository contains source code that I wrote while researching scan matching algorithms, mainly in the period from 2004-2009.
The code of the *mbicp* algorithm was available (Here you will find a patched version). All the other code was written by myself using the description available in published articles as reference.

The main program allows to batch test multiple scan matching algorithms from command line. The results are printed to screen and optionaly PNG images are produced.

Example scans, both real and simulated are provided.

For each algorithm example default configurations are provided in their respective .ini files.

Implemented algorithms include GLASM and its variants, GCP, PolarGA, ICP, MBICP and hybrid variants which combine some genetic and ICP variants for further refinement of the solution.
I implemented other scan matching algorithms in different projects so there is a possibility of adding them here to if I return to the topic.

I am publishing this code in the hope that someone may find it useful. I never had the time to rewrite it properly and to clean it up. Maybe this was the reason why I didn't publish it years ago.
All code used to compile and work a few years ago, but that may not be the case anymore.
If this code helps you in your research then please cite the following article:
"K. Lenac, E. Mumolo, M. Nolich,  Fast Genetic Scan Matching using Corresponding Point Measurements in Mobile Robotics, Proc. EvoIASP 2007"


## Contents of the repository ##

experiments/ - folder where you will setup the experiments writing your configurations. Some example configurations are inside.
images/ - the images produced during simulation are placed here
lib/ - contains the libpngwriter library installation packages copied here for convenience since it is not included anymore in ubuntu repository
src/ - source code. The main program is sm_cli.cpp. Files with "_aux" are wrappers for algorithms providing uniform interface to sm_cli. utils.cpp contains generic utility functions. All other files are algorithm implementations.


## Setting up a development environment with Eclipse ##

### Quick description ###

Install Eclipse-cdt, ecc...
Install required libraries (png, boost_program_options, pngwriter, freetype,...).
Create new project in eclipse (you will need to add required libraries and paths).
Build and run (you may want to configure the Run configurations in eclipse for that)


### Detailed description ###

As example we will install everything starting from a clean fresh Ubuntu 14.04.1 virtual machine. After the installation open the terminal window (Ctrl+Alt+T).

Preliminary step - update system

$ sudo apt-get update
$ sudo apt-get upgrade

Install eclipse with C/C++ development tools

$ sudo apt-get install eclipse-cdt
$ sudo apt-get install build-essential

Install git

$ sudo apt-get install git

Clone the glasm-scan-matching repository from github in your home folder (you can clone it wherever you want. It will be easier to create the eclipse project if you don't clone it in the eclipse workspace folder.)

$ cd
$ git clone https://github.com/klenac/glasm-scan-matching.git

Install required libraries

$ sudo apt-get install boost-program-options-dev libfreetype6-dev

From version Ubuntu 10.10 libpngwriter libraries are not available in the ubuntu repositories. They can be installed manually. Download the last available package for your architecture (http://packages.ubuntu.com/lucid/libpngwriter0-dev) (the packages for amd64 architecture are already included in the lib directory for convenience, so you don't have to download them).

$ cd lib
$ sudo dpkg -i libpngwriter0* 


Create a new project in eclipse

File -> New -> Project -> C++ Project -> Executable -> Empty Project
Untick "Use default location"
Select "Linux GCC" as toolchain
Select the glasm-scan-matching folder cloned from github
choose a name for your project (glasm-scan-matching is fine)

Click Next
Click Advanced settings...

In C/C+ +Build:
Cofiguration: [All configurations]
Builder type: Internal builder


In C/C++ General -> Paths and Symbols
Cofiguration: [All configurations]
Add the following Includes [to all languages]:
/usr/include 
/usr/include/freetype2

Add the following libraries:
png 
pngwriter
boost_program_options
m
freetype

Add the following library paths:
/usr/lib/x86_64-linux-gnu/
/usr/lib



## Building the sm_cli program ##

Set up Run configurations in Eclipse. By default the program will be built automatically when you run it. 
Once built you can also launch the program from command line outside elipse.


## Preparing for the experiments ##


### Type of experiments ###

TODO: describe in detail different types of experiments

ref frame <=> global pose from which scan readings were made (more precisely position of sensors' local ref frame origin)

Inputs to algorithms (different SM algorithms expect the scans to be in different formats):

possible inputs:  
    - pair of scans in same ref frame
    - pair of scans each in it's own ref frame
        - both ref frames given
        - relative displacement of one ref frame to the other given
        - no ref frame given at all, only readings

possible outputs:
    - pose (of new scan, after matching)
    - displacement (motion estimation)
    - multiple hypothesis
        - likelyhood of each (covariance matrix?)



### Preparation (short version) ###

0. make sure the scans are where you expect them to be
1. set up all the options in ./experiments/current/sm.ini
2. for each chosen algorithm set further options in his *.ini file (filename is set in sm.ini)
3. Run

Observe the output on stdout (Console if run from eclipse) and look for images in ./images folder.

### Preparation (detailed version) ###

Typical steps for a new set of tests are as follows.

#### Preliminaries ####


If the ./experiments/current folder is not empty, make it empty. Just copy, backup, archive or delete contents already present there.

If the ./images folder is not empty, make it empty. It may contain images from previous experiments.

#### Prepare the scans ####

Scans may be in one of the supported formats (look the examples or the source code for supported formats). In the basic format scan.pos (x,y,z) holds the true position from which the scan was taken (ground truth) while scan readings hold the readings sorted with increasing angle. The ground truth position is necessary for some types of experiments, but not for all.

Put them in their folder (name it as you like, as long as the you set the path in sm.ini correctly), typically in ./current/scans

#### Configure the simulation (sm.ini) ####

When launched the program reads all the configurations from the 

./experiments/current/

folder. The sm.ini file is read first and depending on the options found there the other configuration files may also be read.

Some example configurations are provided. For start, just copying all the files from one of the example folders in the current/ folder should give you a working configuration (the scans must be in the folder specified in sm.ini).

For an explanation of various options of each configuration file refer to ./experiments/doc/ folder. For example you can start with

./experiments/doc/sm.ini.explained


#### Configure the algorithms ####

For each chosen algorithm set further options in his *.ini file (filename is set in sm.ini).

Refer to documentation in ./experiments/doc/ folder or look in the source code on how to set it up.

#### Drawing ####

To enable the creation of png images in various parts of the project (ie simulation wide) enable DRAW_PNG in "simulation.h" this way we avoid compiling and linking drawing code when DRAW_PNG is not enabled.

When adding new code in project, include "simulation.h" to get DRAW_PNG preference and encircle drawing code with #ifdef DRAW_PNG to include it only if necessary



## Running the tests ##

Run it from Eclipse or launch the program from command line outside elipse.
The results are printed on standard output (console). The images will be created in ./images folder only if you configured their creation in sm.ini and drawing is enabled in simulation.h.


## Analysing the results ##

Once finished running your tests, comment and save the results, explain the experiment. Annotate the changeset / revision of the program used to perform them. Also save the configurations (*.ini) used, the scans, the results and (optional) images, and anything else necessary to reproduce the experiments.

Now you can analyse the results.

## Developing ##

Git is used for revision control. When making changes to source code remember that some files and directories are ignored when committing. Look at .gitignore file. For example the experiments folder is ignored.



