glasm-scan-matching
===================

Implementation of some laser scan matching algorithms and a program to compare and test them.

## Introduction ##

This repository contains source code that I wrote while researching scan matching algorithms, mainly in the period from 2004-2009.
The code of the *mbicp* algorithm was available (Here you will find a patched version). All the other code was written by myself using the description available in published articles as reference.

The command line program for batch testing the scan matching algorithms is called sm_cli. Read the QUICKSTART file for usage instructions.

Example scans, both real and simulated are also provided.

For each algorithm example default configurations are provided in their respective .ini files.



Implemented algorithms include GLASM and its variants, GCP, PolarGA, ICP, MBICP and hybrid variants which combine some genetic and ICP variants for further refinement of the solution.
I implemented several other scan matching algorithms so their code may be added quickly as well if I return to the topic (e-mail me if interested).

I am publishing this code in the hope that someone may find it useful. I never had the time to rewrite it properly and to clean it up. Maybe this was the reason why I didn't publish it years ago.
All code used to compile and work a few years ago, but that may not be the case anymore.
If this code helps you in your research then please cite the following article:
"K. Lenac, E. Mumolo, M. Nolich,  Fast Genetic Scan Matching using Corresponding Point Measurements in Mobile Robotics, Proc. EvoIASP 2007"

