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

#include "gcp.h"
#include "sga.h"
#include <iostream>
#include "simulation.h"
#ifdef DRAW_PNG
    #include <pngwriter.h>
#endif

namespace GCP
{


//---------------------------------------------------------------------------
//---------------------------- GCP ------------------------------------------
//---------------------------------------------------------------------------

// local functions declarations
void GCPobjfun(struct individual *critter);
void GCPchromosome2pos(unsigned c);
void FindCorrPolarFitness(scan& S2, const double corr_distance, double& sumE, unsigned& counter);

// genetic params
int GCPnbitx;
int GCPnbity;
int GCPnbitrot;
int GCPpopsize;
int GCPmaxruns;
int GCPmaxgen;
double GCPpcross;
double GCPpmutation;

unsigned GCPxmask;
unsigned GCPymask;

// scans
scan GCPScanRef;
scan GCPScanNew;

// search area
double GCPMaxX;
double GCPMinX;
double GCPMaxY;
double GCPMinY;
double GCPMaxFI;
double GCPMinFI;;
double GCPStepX;
double GCPStepY;
double GCPStepFI;

// Correspondences Grid (equivalent to lookup table in GLASM)
unsigned CorrGrid_ROWS;			// number of rows (~1000 seems ok for usual env);
unsigned CorrGrid_COLUMNS;		// number of columns;
double CorrGrid_SIZE_X;			// in meters, side x
double CorrGrid_SIZE_Y;			// in meters, side y
double GCP_CORR_DISTANCE_TRESHOLD;	// distance inside which search for correspondence points
double CorrGrid_STEP_X;			// CorrGrid_SIZE_X/CorrGrid_COLUMNS;
double CorrGrid_STEP_Y;			// CorrGrid_SIZE_Y/CorrGrid_ROWS;
double CorrGrid_OX; 			// coord x lower left corner goes in GCPScanRef.pos.x-CorrGrid_SIZE_X/2.0
double CorrGrid_OY; 			// coord y lower left corner goes in GCPScanRef.pos.y-CorrGrid_SIZE_Y/2.0

const float BIGFLOAT=10000.0;
float **CorrGrid;


void setGCPGeneticParameters(int unbitx, int unbity, int unbitrot, int upopsize, int umaxruns, int umaxgen, double upcross, double upmutation)
{	GCPnbitx=unbitx;
	GCPnbity=unbity;
	GCPnbitrot=unbitrot;
	//lchrom=nbitx+nbity+nbitrot;
	GCPpopsize=upopsize;
	GCPmaxruns=umaxruns;
	GCPmaxgen=umaxgen;
	GCPpcross=upcross;
	GCPpmutation=upmutation;
}
//---------------------------------------------------------------------------


void initializeGCPCorrGrid()
{   //bool showMessageBox=true;
	// ALGORITMO SIMILE A YAMANY, IMPLEMENTATO VELOCEMENTE PER CONFRONTO. USO LA FITNESS DI POLARGA DI MARTINEZ (MOLTO SIMILE A YAMANY)
	// discretizziamo una porzione dello spazio globale nell'intorno di GCPScanRef (centrato su GCPScanRef.pos e con grandezza che dipende da SENSOR_RANGE e MAXSTEPDISTANCE (EXPLORATION_STEP_SIZE))
	// la risoluzione della griglia (dimensione della cella) e' un compromesso tra velocita' di preparazione della lookup table (questa procedura - meno celle piu veloce sara')
	// e l'errore di discretizzazione nell'associazione del corrispondente piu' vicino, come anche dell'occupazione di memoria della griglia.

	// allocazione memoria per array (tabella di lookup)

	try {
		CorrGrid = new float*[CorrGrid_ROWS];            // righe dell'array
		for (unsigned j = 0; j < CorrGrid_ROWS; j++)
			CorrGrid[j] = new float[CorrGrid_COLUMNS];   // colonne dell'array
	}
	catch (std::bad_alloc& ba)
	{
	   // cerr << "bad_alloc caught: " << ba.what() << endl;
	}

	// inizializzazione array di lookup
	for (unsigned i = 0; i < CorrGrid_ROWS; i++)
		for (unsigned j = 0; j < CorrGrid_COLUMNS; j++)
		   CorrGrid[i][j] = BIGFLOAT;

	CorrGrid_STEP_X=CorrGrid_SIZE_X/CorrGrid_COLUMNS;
	CorrGrid_STEP_Y=CorrGrid_SIZE_Y/CorrGrid_ROWS;

		// per ogni reading di scan S1 cresci l'area intorno marcando le celle
	// calcola limiti discreti area da crescere (num celle per avere un area che circoscrive corr_distance_trashold) e variabili ausiliarie
	unsigned NX=GCP_CORR_DISTANCE_TRESHOLD/CorrGrid_STEP_X;
	unsigned NY=GCP_CORR_DISTANCE_TRESHOLD/CorrGrid_STEP_Y;
	unsigned NX2=NX/2;
	unsigned NY2=NY/2;
	unsigned MX=CorrGrid_SIZE_X/CorrGrid_STEP_X-NX2;// maximum X, cioe last X cell dove inizia il bordo
	unsigned MY=CorrGrid_SIZE_Y/CorrGrid_STEP_Y-NY2;// maximum Y, cioe last Y cell

	// grow
	SIT i1=GCPScanRef.readings.begin();
	while (i1!=GCPScanRef.readings.end())
	{
		// global coordinates of the reading of S1
		double X=GCPScanRef.pos.x+i1->distance*cos(GCPScanRef.pos.rot+i1->angle);
		double Y=GCPScanRef.pos.y+i1->distance*sin(GCPScanRef.pos.rot+i1->angle);

		// discretization
		unsigned icx=floor((X-CorrGrid_OX)/CorrGrid_STEP_X+0.5);	// X+CorrGrid_OX is dicretized (X espressed in local frame coord)
		unsigned icy=floor((Y-CorrGrid_OY)/CorrGrid_STEP_Y+0.5);	// -

		// some auxiliary variables
		unsigned ix=(icx>NX2?icx-NX2:0);		// initial cell
		unsigned fx=(icx<MX?icx+NX2:MX+NX2);	// final cell
		unsigned iy=(icy>NY2?icy-NY2:0);
		unsigned fy=(icy<MY?icy+NY2:MY+NY2);

		if (fx>=CorrGrid_ROWS) fx=CorrGrid_ROWS-1;
		if (ix<0) ix=0;
		if (fy>=CorrGrid_COLUMNS)	fy=CorrGrid_COLUMNS-1;
		if (iy<0) iy=0;

		for (unsigned i=ix; i<fx; i++)
				for (unsigned j=iy; j<fy; j++)
				{
					float tX=CorrGrid_OX+CorrGrid_STEP_X*i+CorrGrid_STEP_X/2.0;                      // global cell coordinates
					float tY=CorrGrid_OY+CorrGrid_STEP_Y*j+CorrGrid_STEP_Y/2.0;
					float dist=sqrt((X-tX)*(X-tX)+(Y-tY)*(Y-tY));

					if (CorrGrid[i][j]>dist)
					{   CorrGrid[i][j]=dist;
					}
				}
		i1++;
	}
}
//---------------------------------------------------------------------------

void deleteGCPCorrGrid()
{  for (unsigned i = 0; i < CorrGrid_ROWS;  i++)
	   delete[] CorrGrid[i];                 // delete the columns
   delete[] CorrGrid;                        // delete the rows
}
//---------------------------------------------------------------------------

void setGCPSearchArea(double uMaxX, double uMinX, double uMaxY, double uMinY,double uMaxRot, double uMinRot)
{	GCPMaxX=uMaxX;
	GCPMinX=uMinX;
	GCPMaxY=uMaxY;
	GCPMinY=uMinY;
	GCPMaxFI=uMaxRot;
	GCPMinFI=uMinRot;
}
//---------------------------------------------------------------------------

void GCPobjfun(struct individual *critter)
{	GCPchromosome2pos(critter->chrom);
	unsigned counter=0;
	double sumE=0.0;

	FindCorrPolarFitness(GCPScanNew,0.8,sumE,counter);

//	critter->fitness = sumE*(GCPScanNew.readings.size()+1)/counter/counter;
//temp
	if (sumE>0.0) critter->fitness = counter*counter/sumE/(GCPScanNew.readings.size()+1); else critter->fitness=0.0;
// end temp



}
//---------------------------------------------------------------------------

void GCPchromosome2pos(const unsigned c)											// interpreta posizione dal cromosoma   (sapendo che il chromosoma e' < di 32 bit (ci sta in un unsigned di 32 bit) ho semplificato un po')
{   unsigned x,y,fi;
	x=(c&GCPxmask);                                            // primi nbitx bit
	y=c;
	y=y>>GCPnbitx;
	y=(y&GCPymask);
	fi=c;
	fi=fi>>(GCPnbitx+GCPnbity);

	GCPScanNew.pos.x=GCPMinX+x*GCPStepX;
	GCPScanNew.pos.y=GCPMinY+y*GCPStepY;
	GCPScanNew.pos.rot=GCPMinFI+fi*GCPStepFI;
}

//---------------------------------------------------------------------------

void FindCorrPolarFitness(scan& S2, const double corr_distance, double& sumE, unsigned& counter)
{   sumE=0.0;
	counter=0;
	SIT i2=S2.readings.begin();
	while (i2!=S2.readings.end())
	{   // X,Y: S2 reading global coordinates
		double X=S2.pos.x+i2->distance*cos(S2.pos.rot+i2->angle);
		double Y=S2.pos.y+i2->distance*sin(S2.pos.rot+i2->angle);

		// discretisation
		unsigned int ix;
		unsigned int iy;

		// temp
//		double tix=(X-CorrGrid_OX)/CorrGrid_STEP_X+0.5;
//		double tiy=(Y-CorrGrid_OY)/CorrGrid_STEP_Y+0.5;

		// temp debug
//		if (floor(tix)<0.0) 	Application->MessageBox("ix<0",FloatToStrF(tix,ffFixed,7,6).c_str(),0);
//		if (floor(tix)>=CorrGrid_COLUMNS) 	Application->MessageBox("ix>",FloatToStrF(tix,ffFixed,7,6).c_str(),0);
//		if (floor(tiy)<0.0) 	Application->MessageBox("iy<0",FloatToStrF(tiy,ffFixed,7,6).c_str(),0);
//		if (floor(tiy)>=CorrGrid_ROWS)	Application->MessageBox("iy>",FloatToStrF(tiy,ffFixed,7,6).c_str(),0);
		// end temp debug

//		if (tix<0.0) ix=0; else if (tix>=CorrGrid_COLUMNS) ix=CorrGrid_COLUMNS-1; else ix=floor(tix);
//		if (tiy<0.0) iy=0; else if (tiy>=CorrGrid_ROWS) iy=CorrGrid_ROWS-1; else iy=floor(tiy);
		// end temp

		// discretisation
		ix=floor((X-CorrGrid_OX)/CorrGrid_STEP_X+0.5);
		iy=floor((Y-CorrGrid_OY)/CorrGrid_STEP_Y+0.5);

		// lookup
		if (CorrGrid[ix][iy]<corr_distance) {
			counter++;
			sumE+=CorrGrid[ix][iy];
		}
		i2++;
	}
}
//---------------------------------------------------------------------------

void setGCPCorrGridParameters(int uCorrGrid_ROWS, int uCorrGrid_COLUMNS, double uCorrGrid_SIZE_X, double uCorrGrid_SIZE_Y, double uCorrGrid_CENTERX, double uCorrGrid_CENTERY ,double uCORR_DISTANCE_TRESHOLD)
{	CorrGrid_COLUMNS=uCorrGrid_COLUMNS;
	CorrGrid_ROWS=uCorrGrid_ROWS;
	CorrGrid_SIZE_X=uCorrGrid_SIZE_X;
	CorrGrid_SIZE_Y=uCorrGrid_SIZE_Y;
	CorrGrid_OX=uCorrGrid_CENTERX-CorrGrid_SIZE_X/2.0;  // imposto la coord x dell'angolo inf sin dell'area
	CorrGrid_OY=uCorrGrid_CENTERY-CorrGrid_SIZE_Y/2.0;  // imposto la coord y dell'angolo inf sin dell'area
	GCP_CORR_DISTANCE_TRESHOLD=uCORR_DISTANCE_TRESHOLD;
}
//---------------------------------------------------------------------------

void setGCPRefScan(const scan& usref)
{	GCPScanRef=usref;
}
//---------------------------------------------------------------------------

void setGCPNewScan(const scan& usnew)
{	GCPScanNew=usnew;
}
//---------------------------------------------------------------------------

position GCP()
// input parameters (must be set prior to calling GCP):
//	reference and new scan (the provided position of the new scan is not used and can be whatever, it will be estimated by the algorithm.
//		this is different in local scan matching algorithms where this position is used as the initial position estimate and is only corrected by the algorithm)
//	search area
// 	GCP genetic parameters (nbitx,nbity,nbitrot,popsize,maxruns,maxgen,pcross,pmutation)
//	GCP corr grid parameters
{
	sga_parameters(GCPnbitx,GCPnbity,GCPnbitrot,GCPpopsize,GCPmaxruns,GCPmaxgen,GCPpcross,GCPpmutation);	// set genetic parameters

	GCPxmask=(1<<GCPnbitx)-1;
	GCPymask=(1<<GCPnbity)-1;

	unsigned ResX=1<<GCPnbitx;
	GCPStepX=(GCPMaxX-GCPMinX)/ResX;
	unsigned ResY=1<<GCPnbity;
	GCPStepY=(GCPMaxY-GCPMinY)/ResY;
	unsigned ResFI=1<<GCPnbitrot;
	GCPStepFI=(GCPMaxFI-GCPMinFI)/ResFI;

	set_sga_objfun(&GCPobjfun);

	bestever bestfit=sga();
	GCPchromosome2pos(bestfit.chrom);

	return GCPScanNew.pos;
}
//---------------------------------------------------------------------------

void initializeGCPCorrGrid_ORIGINAL()
{
// This procedure is never used in testing, apart for testing the speed of the original Yamany implementation

	try {
		CorrGrid = new float*[CorrGrid_ROWS];            // righe dell'array
		for (unsigned j = 0; j < CorrGrid_ROWS; j++)
			CorrGrid[j] = new float[CorrGrid_COLUMNS];   // colonne dell'array
	}
	catch (std::bad_alloc& ba) {
	}

	// inizializzazione array di lookup
	for (unsigned i = 0; i < CorrGrid_ROWS; i++)
		for (unsigned j = 0; j < CorrGrid_COLUMNS; j++)
		   CorrGrid[i][j] = BIGFLOAT;

	CorrGrid_STEP_X=CorrGrid_SIZE_X/CorrGrid_COLUMNS;
	CorrGrid_STEP_Y=CorrGrid_SIZE_Y/CorrGrid_ROWS;

	// grow
	SIT i1=GCPScanRef.readings.begin();
	while (i1!=GCPScanRef.readings.end())
	{
		// global coordinates of the reading of S1
		double X=GCPScanRef.pos.x+i1->distance*cos(GCPScanRef.pos.rot+i1->angle);
		double Y=GCPScanRef.pos.y+i1->distance*sin(GCPScanRef.pos.rot+i1->angle);

		// grow around the reading of S1
		for (unsigned i=0; i<CorrGrid_ROWS; i++)
				for (unsigned j=0; j<CorrGrid_COLUMNS; j++)
				{
					float tX=CorrGrid_OX+CorrGrid_STEP_X*i+CorrGrid_STEP_X/2.0;                      // global cell coordinates
					float tY=CorrGrid_OY+CorrGrid_STEP_Y*j+CorrGrid_STEP_Y/2.0;
					float dist=sqrt((X-tX)*(X-tX)+(Y-tY)*(Y-tY));
					if (CorrGrid[i][j]>dist)
					{   CorrGrid[i][j]=dist;
					}
				}
		i1++;
	}
}
//---------------------------------------------------------------------------


void initializeGCPCorrGrid_SUGGESTED()
{
// This procedure is never used in testing, apart for testing the speed of the suggested Yamany implementation
// Yamany suggested to speed up the GCP transform preparation (initialization of the CorrGrid) by applying different cell sizes when near and far from the scan


// actually not of much interest for our SM since the area around the scan to which apply the gcp transform is never small in comparison to the size of the corrgrid

	try {
		CorrGrid = new float*[CorrGrid_ROWS];            // righe dell'array
		for (unsigned j = 0; j < CorrGrid_ROWS; j++)
			CorrGrid[j] = new float[CorrGrid_COLUMNS];   // colonne dell'array
	}
	catch (std::bad_alloc& ba) {
	}

	// inizializzazione array di lookup
	for (unsigned i = 0; i < CorrGrid_ROWS; i++)
		for (unsigned j = 0; j < CorrGrid_COLUMNS; j++)
		   CorrGrid[i][j] = BIGFLOAT;

	CorrGrid_STEP_X=CorrGrid_SIZE_X/CorrGrid_COLUMNS;
	CorrGrid_STEP_Y=CorrGrid_SIZE_Y/CorrGrid_ROWS;



// TODO:

}
//---------------------------------------------------------------------------
} // namespace
