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


//#include <iostream>
//using namespace std;


#include <stdlib.h>

#include "sga.h"
#include "utils.h"


// OTTIMIZZAZIONI POSSIBILI:
//- ipotizzando numero di bit per x,y,z <=8 si puo usare sempre una word di 32 bit per memorizzare un cromosoma. Togliere chromsize ecc ...
//- invece della memoria dinamica usare var globali
//- invece del codice diviso in sub e funzioni, scriverlo inline in un unico blocco
//- usare funzioni random ottimizzate
//- fare lookup per trasformazione coordinate polari in cartesiane in obj
// - allocates auxiliary memory for stochastic remainder selection. E' necessario?




//---------------------------------------------------------------------------
//------------------------ sga ----------------------------------------------
//---------------------------------------------------------------------------



void generation(void);
void initialize(void);
void initmalloc();// memory allocation of space for global data structures
void freeall(); // A routine to free all the space dynamically allocated in initspace()
void mutation(unsigned &child); // Mutate an allele w/ pmutation, count # of mutations
int crossover (const unsigned &parent1, const unsigned &parent2, unsigned &child1, unsigned &child2);// Cross 2 parent strings, place in 2 child strings
void select_memory(); // allocates auxiliary memory for stochastic remainder selection
void select_free(); // frees auxiliary memory for stochastic remainder selection
void preselect(); // preselection for stochastic remainder method
int select(); // selection using remainder method
void statistics(struct individual *pop);// Calculate population statistics

// definizioni globali
int nbitx;
int nbity;
int nbitrot;
int popsize;
int maxruns;
int maxgen;
double pcross;
double pmutation;

const int BITS_PER_BYTE = 8;// number of bits per byte on this machine
//const int UINTSIZE = BITS_PER_BYTE*sizeof(unsigned int);// # of bits in unsigned

int lchrom;

struct individual *oldpop;/* last generation of individuals */
struct individual *newpop;/* next generation of individuals */
struct bestever bestfit; /* fittest individual so far */
double sumfitness;/* summed fitness for entire population */
double mymax;/* maximum fitness of population */
double avg;/* average fitness of population */
double mymin;/* minumum fitness of population */

int gen;/* current generation number */
int run;/* current run number */
int nmutation;/* number of mutations */
int ncross;/* number of crossovers */


void (*sga_objfun)(struct individual *critter);		// pointer to application dependent objective function

#ifdef DRAW_PNG

	void (*sga_drawfun)(int run, int gen, int popsize, struct individual *critter, struct bestever *bestfit);	// pointer to application dependent drawing function

	void set_sga_drawfun(void (*drawfun)(int run, int gen, int popsize, struct individual *critter, struct bestever *bestfit))
	{	sga_drawfun=drawfun;
	}

#endif

//---------------------------------------------------------------------------
//------------------------ sga ----------------------------------------------
//---------------------------------------------------------------------------

void set_sga_objfun(void (*objfun)(struct individual *critter))
{	sga_objfun=objfun;
}

void sga_parameters(int unbitx, int unbity, int unbitrot, int upopsize, int umaxruns, int umaxgen, double upcross, double upmutation)
{	nbitx=unbitx;
	nbity=unbity;
	nbitrot=unbitrot;
	lchrom=nbitx+nbity+nbitrot;

	popsize=upopsize;
	if((popsize%2) != 0) popsize++;
	maxruns=umaxruns;
	maxgen=umaxgen;
	pcross=upcross;
	pmutation=upmutation;
}
//---------------------------------------------------------------------------

bestever sga(void) // la funzione principale da chiamare dopo aver settato tutti i parametri del
{	initmalloc(); // malloc space for global data structures
	nmutation = 0;// initialize global counters/values
	ncross = 0;
	bestfit.fitness = 0.0;
	bestfit.generation = 0;
	for(run=1; run<=maxruns; run++)
	{	initialize();
#ifdef DRAW_PNG
		//sga_drawfun(run, 0, popsize, oldpop, &bestfit);
#endif
		for(gen=0; gen<maxgen; gen++)
		{	generation(); // create a new generation
#ifdef DRAW_PNG
			sga_drawfun(run, gen, popsize, oldpop, &bestfit);
#endif
			statistics(newpop); // compute fitness statistics on new populations
			struct individual *temp = oldpop; // advance the generation
			oldpop = newpop;
			newpop = temp;


		}
	}
	freeall();
	return bestfit;
}
//---------------------------------------------------------------------------

void generation()
{	int mate1, mate2, jcross, j = 0;
	preselect();// perform any preselection actions necessary before generation
	do// select, crossover, and mutation
	{ 	mate1 = select(); // pick a pair of mates
		mate2 = select();
		jcross = crossover(oldpop[mate1].chrom, oldpop[mate2].chrom, newpop[j].chrom, newpop[j+1].chrom);// Crossover and mutation
		mutation(newpop[j].chrom);
		mutation(newpop[j+1].chrom);
		sga_objfun(&(newpop[j]));// Decode string, evaluate fitness, & record parentage date on both children
		newpop[j].parent[0] = mate1+1;
		newpop[j].xsite = jcross;
		newpop[j].parent[1] = mate2+1;
		sga_objfun(&(newpop[j+1]));
		newpop[j+1].parent[0] = mate1+1;
		newpop[j+1].xsite = jcross;
		newpop[j+1].parent[1] = mate2+1;
		j = j + 2;// Increment population index */
	} while(j < (popsize-1));
}
//---------------------------------------------------------------------------

void initialize() // Initialization Coordinator
{	int j, k;// initialize population
	unsigned mask = 1;
	for(j = 0; j < popsize; j++)
	{ 	oldpop[j].chrom = 0;
		for(k = 0; k < lchrom; k++)
		{	oldpop[j].chrom = (oldpop[j].chrom<<1);
			if(_random(1000)<499) oldpop[j].chrom = oldpop[j].chrom|mask;
		}
		oldpop[j].parent[0] = 0;// Initialize parent info
		oldpop[j].parent[1] = 0;
		oldpop[j].xsite = 0;
		sga_objfun(&(oldpop[j]));// Evaluate initial fitness
	}
	statistics(oldpop);
}
//---------------------------------------------------------------------------

void initmalloc() // memory allocation of space for global data structures
{ 	unsigned nbytes;

	nbytes = popsize*sizeof(struct individual); // memory for old and new populations of individuals
	if((oldpop = (struct individual *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per oldpop
	if((newpop = (struct individual *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per newpop

	select_memory();// allocate any auxiliary memory for selection
}
//---------------------------------------------------------------------------

void freeall() // A routine to free all the space dynamically allocated in initspace()
{ 	free(oldpop);
	free(newpop);
	select_free();// free any auxiliary memory needed for selection
}
//---------------------------------------------------------------------------

void mutation(unsigned &child)// Mutate an allele w/ pmutation, count # of mutations
{ 	int k;
	unsigned mask, temp = 1;

	for(k = 0; k < lchrom; k++)
	{ 	mask = 0;
		if(_random(1000)<pmutation*1000)
		{ 	mask = mask|(temp<<k);
			nmutation++;
		}
		child = child^mask;
	}
}
//---------------------------------------------------------------------------

int crossover (const unsigned &parent1, const unsigned &parent2, unsigned &child1, unsigned &child2)// Cross 2 parent strings, place in 2 child strings
{	int jcross;
	unsigned mask;

	if(_random(1000)<pcross*1000)// Do crossover with probability pcross
	{ 	jcross = _random(lchrom - 1)+1; // Cross between 1 and l-1
		ncross++;
		mask=1;
		mask<<jcross;
		mask-=1;
		child1 = (parent1&mask)|(parent2&(~mask));
		child2 = (parent1&(~mask))|(parent2&mask);
	} else
	{ 	child1 = parent1;
		child2 = parent2;
		jcross = 0;
	}
	return jcross;
}
//---------------------------------------------------------------------------

static int *choices, nremain;
static float *fraction;

void select_memory()// allocates auxiliary memory for stochastic remainder selection
{ 	unsigned nbytes;
	//char *malloc();
	nbytes = popsize*sizeof(int);
	if((choices = (int *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per choices
	nbytes = popsize*sizeof(float);
	if((fraction = (float *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per fraction
}
//---------------------------------------------------------------------------

void select_free()// frees auxiliary memory for stochastic remainder selection
{ 	//free(choices); // temp per evitare con polarga segfault
    free(choices);
	free(fraction);
}
//---------------------------------------------------------------------------

void preselect()// preselection for stochastic remainder method
{
	int j, jassign, k;
	float expected;
	////////////////////////////////


	if (avg == 0)
	{
		for(j = 0; j < popsize; j++) choices[j] = j;
	}
	else
	{	j = 0;
		k = 0;
		do// Assign whole numbers
		{
			expected = (float)(oldpop[j].fitness)/avg;
			jassign = (int)expected;
			fraction[j] = expected - jassign; // note that expected is automatically truncated
			while(jassign > 0)
			{	jassign--;
				choices[k] = j;
				k++;
			}
			j++;
		} while(j < popsize);

		j = 0;
		while(k < popsize)// Assign fractional parts
		{ 	if(j >= popsize) j = 0;
			if(fraction[j] > 0.0)
			{ 	if(_random(1000)<fraction[j]*1000) // A winner if true
				{ 	choices[k] = j;
					fraction[j] = fraction[j] - 1.0;
					k++;
				}
			}
			j++;
		}
	}
	nremain = popsize - 1;
}
//---------------------------------------------------------------------------

int select()// selection using remainder method
{	int jpick, slect;
	jpick = _random(nremain);
	slect = choices[jpick];
	choices[jpick] = choices[nremain];
	nremain--;
	return(slect);
}
//---------------------------------------------------------------------------

void statistics(struct individual *pop)// Calculate population statistics
{ 	int j;
	sumfitness = 0.0;
	mymin = pop[0].fitness;
	mymax = pop[0].fitness;
	for(j = 0; j < popsize; j++)// Loop for mymax, mymin, sumfitness
	{	sumfitness = sumfitness + pop[j].fitness;// Accumulate
		if(pop[j].fitness > mymax) mymax = pop[j].fitness;// New maximum
		if(pop[j].fitness < mymin) mymin = pop[j].fitness;// New minimum
		if(pop[j].fitness > bestfit.fitness)// new global best-fit individual
		{ 	bestfit.chrom = pop[j].chrom;
			bestfit.fitness= pop[j].fitness;
			bestfit.generation = gen;
		}
	}
	avg = sumfitness/popsize;// Calculate average
}
//---------------------------------------------------------------------------

