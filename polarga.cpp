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

#include "polarga.h"
#include "sga.h"
#include <iostream>

#ifdef DRAW_PNG
	#ifndef PNGWRITER_H
		#include <pngwriter.h>
	#endif

#include "mydraw.h" 	// for drawing
const int FILENAME_SIZE = 100;
#endif


namespace POLARGA
{


#ifdef DRAW_PNG
	bool draw_individual;
	bool draw_generations;
	unsigned verbose_level;

	drawing img;

	void setDebugParameters(bool udraw_individual, bool udraw_generations, unsigned uverbose_level)
	{   draw_individual=udraw_individual;
		draw_generations=udraw_generations;
		verbose_level=uverbose_level;
	}

#endif



//---------------------------------------------------------------------------
//------------------------ PolarGA --------------------------------------------
//---------------------------------------------------------------------------

// local proc declarations
//void PolarGAobjfunc(struct individual *critter);
void PolarGAchromosome2pos(unsigned c);
//string printChrom(unsigned chrom);
//---------------------------------------------------------------------------

// genetic params and proc to set them from outside
int PolarGAnbitx;
int PolarGAnbity;
int PolarGAnbitrot;
int PolarGApopsize;
int PolarGAmaxruns;
int PolarGAmaxgen;
double PolarGApcross;
double PolarGApmutation;

double dangle;
//---------------------------------------------------------------------------

void setPolarGAGeneticParameters(int unbitx, int unbity, int unbitrot, int upopsize, int umaxruns, int umaxgen, double upcross, double upmutation)
{	PolarGAnbitx=unbitx;
	PolarGAnbity=unbity;
	PolarGAnbitrot=unbitrot;
	PolarGApopsize=upopsize;
	PolarGAmaxruns=umaxruns;
	PolarGAmaxgen=umaxgen;
	PolarGApcross=upcross;
	PolarGApmutation=upmutation;
}
//---------------------------------------------------------------------------
double PolarGA_CORR_DISTANCE_TRESHOLD;  // distanza entro la quale cercare i corrispondenti
unsigned PolarGA_NREADINGS; // number of readings that laser is capable of in one scan
double PolarGA_FOV;
bool PolarGA_RELAXED_CORR_SEARCH;	// set to true to get more correspondences

void setPolarGAParameters(unsigned numreadingsinscan, double field_of_view, double corr_distance_treshold, bool relaxed_corr_search)
{   PolarGA_NREADINGS=numreadingsinscan;
	PolarGA_CORR_DISTANCE_TRESHOLD=corr_distance_treshold;
	PolarGA_FOV=field_of_view;
	PolarGA_RELAXED_CORR_SEARCH=relaxed_corr_search;
}



unsigned PolarGAxmask;
unsigned PolarGAymask;
unsigned PolarGArotmask;
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
// scans and proc to set them from outside
scan PolarGAScanRef;
scan PolarGAScanNew;
double sref[1000];
double snew[1000];

void setPolarGARefScan(const scan& usref)
{
	PolarGAScanRef.copy(usref);
	for(unsigned i=0; i<PolarGA_NREADINGS;i++) sref[i]=0.0;

	double dangle=PolarGA_FOV/PolarGA_NREADINGS;
	CSIT si;
	si=PolarGAScanRef.readings.begin();
	while (si!=PolarGAScanRef.readings.end())
	{
		unsigned index = int(angleUnwrap(si->angle)/dangle+0.5);
		sref[index]=si->distance;
		si++;
	}
}


void setPolarGANewScan(const scan& usnew)
{
	PolarGAScanNew.copy(usnew);
	for(unsigned i=0; i<PolarGA_NREADINGS;i++) snew[i]=0.0;

	double dangle=PolarGA_FOV/PolarGA_NREADINGS;
	CSIT si;
	si=PolarGAScanNew.readings.begin();
	while (si!=PolarGAScanNew.readings.end())
	{
		unsigned index = int(angleUnwrap(si->angle)/dangle+0.5);
		snew[index]=si->distance;
		si++;
	}
}

// ATTN: readings are not copied during assignment, only pointer to readings, so
// after the assignment both scans use the same list of readings.
// The readings are not modified by PolarGA so it is ok, but it is not nice to use
// assignment like this. We should define the assignment operator in scan class
// that copies the list of readings.

//---------------------------------------------------------------------------
// search area and proc to set it from outside
double PolarGAMaxX;
double PolarGAMinX;
double PolarGAMaxY;
double PolarGAMinY;
double PolarGAMaxFI;
double PolarGAMinFI;;
double PolarGAStepX;
double PolarGAStepY;
double PolarGAStepFI;


void setPolarGASearchArea(double uMinX, double uMaxX, double uMinY, double uMaxY,double uMinRot, double uMaxRot)
{	PolarGAMaxX=uMaxX;
PolarGAMinX=uMinX;
PolarGAMaxY=uMaxY;
PolarGAMinY=uMinY;
if (uMaxRot<uMinRot) {
	PolarGAMinFI=uMaxRot;
	PolarGAMaxFI=uMinRot;
} else {
	PolarGAMaxFI=uMaxRot;
	PolarGAMinFI=uMinRot;
}
}

// temp
unsigned ofcounter;
// end temp

#ifdef DRAW_PNG


void PolarGAdrawfun(int run, int gen, int popsize, struct individual *newpop, struct bestever *bestfit)
{
	if (draw_generations)
	{   char filename[FILENAME_SIZE];
			snprintf(filename,FILENAME_SIZE,"./images/polarga_run_%02d_gen_%03d.png",run, gen);

			// create image
			drawing img(filename,25.0,25.0,-12.5,-12.5,"./cfg/bitmaps/empty_1000x1000.png");
			img.grid();
			for (int ii=0; ii<popsize; ii++) {
				PolarGAchromosome2pos(newpop[ii].chrom);
				img.setcolor(10000,10000,1000);
				img.readings(PolarGAScanNew, PolarGAScanNew.pos);
			}
			PolarGAchromosome2pos(bestfit->chrom);
			img.setcolor(60000,0,20000);
			img.readings(PolarGAScanNew, PolarGAScanNew.pos);
			img.write();
			std::cerr << "PolarGAMobj bestfit = " << bestfit->fitness << std::endl;
	}
}


#endif






void PolarGAobjfun(struct individual *critter)
{

// todo: questa funzione si puo implementare piu efficacemente: invece di copiare
// polargascannew (per lavorare con la copia) e poi usare l'inefficiente
// readingsChangeRefFrame(), calcolare le corrispondenze direttamente con i
// calcoli nell'array locale senza modificare polargascannew
// verificare con oprofile

#ifdef DRAW_PNG
bool test=ofcounter<200; // limit to 200 images

	if (draw_individual && test)
	{   char filename[FILENAME_SIZE];
		snprintf(filename,FILENAME_SIZE,"./images/polarGAobjfun%03d.png",ofcounter);
		img.init(filename,25.0,25.0,-12.5,-12.5,"./cfg/bitmaps/empty_1000x1000.png");
		img.grid();
		img.setcolor(10000,10000,30000);
		img.readings(PolarGAScanRef, PolarGAScanRef.pos);
		//img.setcolor(50000,50000,50000);
		//img.readings(PolarGAScanNew, PolarGAScanNew.pos);
	}
#endif
	//std::cout << "starting polarGAobjfun pnew.pos(" << PolarGAScanNew.pos.x << "," << PolarGAScanNew.pos.y << "," << PolarGAScanNew.pos.rot << std::endl;
position oldPos=PolarGAScanNew.pos;
PolarGAchromosome2pos(critter->chrom);

// per mantenere i readings di PolarGAScanNew inalterati in quanto servono per tutta la popolazione uso un'altra scan
scan tScan;
tScan.copy(PolarGAScanNew);

// riposiziona le letture perche scannew.pos e' cambiata valutando il cromosoma
	//readingsChangeRefFrame(tScan.readings, oldPos, tScan.pos);
#ifdef DRAW_PNG
	if (draw_individual && test)
	{   img.setcolor(0,0,10000);
		img.readings(tScan, tScan.pos);
		//std::cout << "tScan.pos(" << tScan.pos.x << "," << tScan.pos.y << "," << tScan.pos.rot << std::endl;
		img.setcolor(30000,30000,30000);
		img.readings(PolarGAScanNew, PolarGAScanNew.pos);
	}
#endif

	// riposiziona le letture rispetto allo stesso sist di rif di ScanRef
	readingsChangeRefFrame(tScan.readings, tScan.pos, PolarGAScanRef.pos);
#ifdef DRAW_PNG
	if (draw_individual && test)
	{   img.setcolor(50000,10000,0000);
		img.readings(tScan, PolarGAScanRef.pos);
	}
#endif

unsigned counter=0;
double sumE=0.0;

for(unsigned i=0; i<PolarGA_NREADINGS*2;i++) snew[i]=0.0;
double dangle=PolarGA_FOV/PolarGA_NREADINGS;

	// calcola per ogni lettura di scannew distance e angle nel sist di rif di scanref e salva in snew[]
	CSIT i2=tScan.readings.begin();
	while ( i2 != tScan.readings.end() )
	{   unsigned index = int(angleUnwrap(i2->angle)/dangle+0.5);
		snew[index]=i2->distance;
		i2++;
	}

	// relax
	if (PolarGA_RELAXED_CORR_SEARCH)
	{
		for(unsigned i=1; i<PolarGA_NREADINGS*2-1;i++)
		{   if (snew[i]==0.0 && snew[i+1]>0.0 && snew[i-1]>0.0) snew[i]=(snew[i-1]+snew[i+1])/2.0;
			if (sref[i]==0 && sref[i+1]>0.0 && sref[i-1]>0.0) sref[i]=(sref[i-1]+sref[i+1])/2.0;
		}
	}

	for(unsigned i=0; i<PolarGA_NREADINGS*2;i++)
	{
		if (snew[i]==0 || sref[i]==0) continue;
		double dist=snew[i]-sref[i];
		if (dist<0.0) dist= -dist; // absolute value
		//std::cout << "  sref[" << i << "]=" << sref[i] << "  snew[" << i << "]=" << snew[i];
		//if (dist>0.0 && dist<PolarGA_CORR_DISTANCE_TRESHOLD) std::cout << "---->  bingo" << std::endl; else std::cout << std::endl;
		if (dist>0.0 && dist<PolarGA_CORR_DISTANCE_TRESHOLD)
		{
			counter++;
			sumE+=dist;

#ifdef DRAW_PNG
		if (draw_individual && test)
		{
			double X1=PolarGAScanRef.pos.x+snew[i]*cos(PolarGAScanRef.pos.rot+i*dangle);
			double Y1=PolarGAScanRef.pos.y+snew[i]*sin(PolarGAScanRef.pos.rot+i*dangle);
			double X2=PolarGAScanRef.pos.x+sref[i]*cos(PolarGAScanRef.pos.rot+i*dangle);
			double Y2=PolarGAScanRef.pos.y+sref[i]*sin(PolarGAScanRef.pos.rot+i*dangle);

//             double X1=snew[i]*cos(i*dangle);
//             double Y1=snew[i]*sin(i*dangle);
//             double X2=sref[i]*cos(i*dangle);
//             double Y2=sref[i]*sin(i*dangle);



			img.setcolor(10000,50000,10000);
			img.line(PolarGAScanRef.pos.x,PolarGAScanRef.pos.y,X2,Y2);
			img.setcolor(50000,50000,30000);
			img.line(PolarGAScanRef.pos.x,PolarGAScanRef.pos.y,X1,Y1);
			img.setcolor(50000,10000,10000);
			img.line(X1,Y1,X2,Y2);
			img.setpointsize(1);
	img.setcolor(10000,10000,10000);
	img.point(X1,Y1);
	img.setcolor(40000,40000,40000);
	img.point(X2,Y2);
		}
		//std::cerr << "angle diff = " << i1->angle-i2->angle << "distance = " << dist << std::endl;
		//flag=true;
#endif
		//std::cout << "PolarGAMobj fitness = " << counter << " sum " << sumE << std::endl;
		//std::cout << "counter = " << counter << " sum " << sumE << std::endl;
		}
	}

//if (counter>20) critter->fitness = sumE*PolarGA_NREADINGS/counter/counter;    // original
//if (counter>20) critter->fitness = sumE*PolarGA_NREADINGS/counter;
//if (counter>20) critter->fitness = 1-sumE/counter;
	if (counter>2) critter->fitness = 1/(sumE*PolarGA_NREADINGS/counter/counter);    // original
	else {
	critter->fitness = 0.001;
	counter=1;
	}

//  std::cout << ofcounter << ": Fitness = " << critter->fitness
//      << " counter " << counter << " sum " << sumE << " pos: "
//      << PolarGAScanNew.pos.x << "," << PolarGAScanNew.pos.y << "," << PolarGAScanNew.pos.rot << std::endl;
// std::cout << tScan.pos.x << "," << tScan.pos.y << "," << tScan.pos.rot << "," << std::endl;


#ifdef DRAW_PNG
	if (draw_individual && test)
	{
		const unsigned NAMESIZE=255;
		char tmp[NAMESIZE];
		memset(tmp,0,NAMESIZE);
		strncpy(tmp,"PolarGAMobj fitness = ",22);
		sprintf(&tmp[22],"%d - fitness = %f, counter = %d, sum = %f, pos = %f,%f,%f",ofcounter, critter->fitness,counter,sumE,PolarGAScanNew.pos.x,PolarGAScanNew.pos.y,PolarGAScanNew.pos.rot);
		img.setcolor(0,0,0);
		img.text(300,805,tmp);
		img.write();
		//if (flag) exit(1);
		//std::cout << "PolarGAMobj fitness = " << critter->fitness <<" counter " << counter << " sum " << sumE << std::endl;
		//char ch; std::cout << "Press key";
		//std::cin >> ch;
	}
#endif
	ofcounter++;

}

//---------------------------------------------------------------------------




void PolarGAobjfun2(struct individual *critter)
{
PolarGAchromosome2pos(critter->chrom);

unsigned counter=0;
double sumE=0.0;


// calcola per ogni lettura di scannew distance e angle nel sist di rif di scanref e salva in snew[]
// le letture di snew sono locali a snew e vanno prima spostate in

CSIT i2=PolarGAScanNew.readings.begin();
while ( i2 != PolarGAScanNew.readings.end() )
{

	// coordinate X,Y gobali del reading di S2
	double X=PolarGAScanNew.pos.x+i2->distance*cos(PolarGAScanNew.pos.rot+i2->angle);
	double Y=PolarGAScanNew.pos.y+i2->distance*sin(PolarGAScanNew.pos.rot+i2->angle);

	//  std::cout << "X,Y" << X << " " << Y << std::endl;

	double dist=sqrt(X*X+Y*Y);
	double angle=atan2(Y,X);

	//std::cout << "dist,angle,dangle" << dist << " " << angle << " " << dangle << std::endl;


	unsigned index = int(angleUnwrap(angle)/dangle+0.5);

	double rdist=dist-sref[index];
	if (rdist<0.0) rdist= -rdist; // absolute value
	if (rdist>0.0 && rdist<PolarGA_CORR_DISTANCE_TRESHOLD)
	{
	counter++;
	sumE+=rdist;
	}
	i2++;
}

	if (counter>2) critter->fitness = 1/(sumE*PolarGA_NREADINGS/counter/counter);    // original
	else {
	critter->fitness = 0.001;
	counter=1;
	}
	ofcounter++;
}

//---------------------------------------------------------------------------






void PolarGAchromosome2pos(const unsigned c)// interpreta posizione dal cromosoma (sapendo che il chromosoma e' < di 32 bit (ci sta in un unsigned di 32 bit) ho semplificato un po')
{
	unsigned x,y,fi;
	//std::cerr << "PolarGA2pos c: " << c << std::endl;
	x=(c&PolarGAxmask);// primi nbitx bit
	y=c;
	y=y>>PolarGAnbitx;
	y=(y&PolarGAymask);
	//std::cerr << "PolarGA2pos-0001--->" << std::endl;
	fi=c;
	fi=fi>>(PolarGAnbitx+PolarGAnbity);
	fi=(fi&PolarGArotmask);

	PolarGAScanNew.pos.x=PolarGAMinX+x*PolarGAStepX;
	PolarGAScanNew.pos.y=PolarGAMinY+y*PolarGAStepY;
	PolarGAScanNew.pos.rot=angleUnwrap(PolarGAMinFI+fi*PolarGAStepFI);

	//std::cerr << "PolarGAMinX " << PolarGAMinX << "PolarGAMinY " << PolarGAMinY << "PolarGAMinFI " << PolarGAMinFI << "PolarGAStepX " << PolarGAStepX << std::endl;

	//std::cerr << "c: " << c << " x,y,rot = (" << x << "," << y << "," << fi << ")" << " Position x,y,rot = (" << PolarGAScanNew.pos.x << "," << PolarGAScanNew.pos.y << "," << PolarGAScanNew.pos.rot << ")" << std::endl;

}
//---------------------------------------------------------------------------

position PolarGA(double* outfitness)
// input parameters (must be set prior to calling PolarGA):
//	reference and new scan (the provided position of the new scan is not used and can be whatever, it will be estimated by the algorithm.
//		this is different in local scan matching algorithms where this position is used as the initial position estimate and is only corrected by the algorithm)
//	search area
// 	PolarGA genetic parameters (nbitx,nbity,nbitrot,popsize,maxruns,maxgen,pcross,pmutation)
{


sga_parameters(PolarGAnbitx,PolarGAnbity,PolarGAnbitrot,PolarGApopsize,PolarGAmaxruns,PolarGAmaxgen,PolarGApcross,PolarGApmutation);	// set genetic parameters

PolarGAxmask=(1<<PolarGAnbitx)-1;
PolarGAymask=(1<<PolarGAnbity)-1;
PolarGArotmask=(1<<PolarGAnbitrot)-1;

unsigned ResX=1<<PolarGAnbitx;
PolarGAStepX=(PolarGAMaxX-PolarGAMinX)/ResX;

unsigned ResY=1<<PolarGAnbity;
PolarGAStepY=(PolarGAMaxY-PolarGAMinY)/ResY;
unsigned ResFI=1<<PolarGAnbitrot;
double anglediff=PolarGAMaxFI-PolarGAMinFI;
	if (anglediff<0) anglediff+=2*M_PI;
	PolarGAStepFI=anglediff/ResFI;
set_sga_objfun(&PolarGAobjfun2);

dangle=PolarGA_FOV/PolarGA_NREADINGS;

#ifdef DRAW_PNG
set_sga_drawfun(&PolarGAdrawfun);
// temp
ofcounter=0;
// endtemp
#endif
//std::cout << "starting sga" << std::endl;
bestever bestfit=sga();
  *outfitness = bestfit.fitness;
//std::cout << "Best fitness = " << bestfit.fitness << " Position x,y,rot = (" << PolarGAScanNew.pos.x << "," << PolarGAScanNew.pos.y << "," << PolarGAScanNew.pos.rot << ")" << std::endl;
PolarGAchromosome2pos(bestfit.chrom);
//std::cout << "Best fitness = " << bestfit.fitness << " Position x,y,rot = (" << PolarGAScanNew.pos.x << "," << PolarGAScanNew.pos.y << "," << PolarGAScanNew.pos.rot << ")" << std::endl;
return PolarGAScanNew.pos;
}


} // namespace
