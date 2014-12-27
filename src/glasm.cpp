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


#include "glasm.h"
#include "sga.h"
#include <iostream>

#ifdef DRAW_PNG
#include <pngwriter.h>
#include "mydraw.h"     // for drawing
const int FILENAME_SIZE = 100;
#endif

/* TODO:
	ideas for speed optimization:
	a) sin/cos	lookup (No, because probably there will be no speed gain with FPU already very fast
		and slow memory)
	b) faster float	to	int (look at the function float2int in the source. Worth a try)
	c) integer or fixed point arithmetic for all fitness evaluations (discretise all possible values
	that a scan reading may assume (OK, worth a try).  (reasonable for a rough estimation when used
	in hybrid algorithm.)
	d) after discretization: precalculation with lookup table for cos(nrp.angle) and
	sin (nrp.angle)
	e) parallelization (easy for all GA distributing the population to multiple cores)
	f) procedure for (semi)automatic genetic parameters adjustment/optimization
	g) depending on application: Search space size can be reduced if a mask of valid positions, ie
	chromosomes, is supplied. It is a binary yes/no associative array/direct address/lookup. For
	invalid chromosomes one simple solution is to just skip them assigning 0/negative fitness
	h) chromosome cache - if already evaluated use stored fitness (good for dense areas. Since it
	typically happens around high fitness areas, only high fitness chromosomes could be cached.
	We should first study the number hits/misses for typical cases and see if it pays off.
	i) SSE single instruction multiple data (needs exploring ...)
	j) reordering of execution to study effects on cache lines, pipeline ... (needs exploring..)
 */



namespace GLASM
{
//---------------------------------------------------------------------------
//------------------------ GLASM --------------------------------------------
//---------------------------------------------------------------------------

// local proc declarations
//void objfunc(struct individual *critter);
bool chromosome2pos(unsigned c);
//void initializeInvGrayTable(void);
//void initializeLookupTable(void);
//void deleteLookupTable(void);
//void initializeGrayTable();
//string printChrom(unsigned chrom);
//---------------------------------------------------------------------------

// todo: put this as parameter
const unsigned outlier_treshold=200;

//#ifdef DRAW_PNG
bool draw_lookup;
bool draw_individual;
bool draw_generations;
unsigned verbose_level;
std::string lookup_bitmap;
std::string valid_bitmap;
double map_size_x;
double map_size_y;
std::string map_filename;

void setDebugParameters(bool udraw_lookup,
	bool udraw_individual,
	bool udraw_generations,
	unsigned uverbose_level,
	double umap_size_x,
	double umap_size_y,
	std::string umap_filename,
	std::string ulookup_bitmap,
	std::string uvalid_bitmap)
{   draw_lookup=udraw_lookup;
	draw_individual=udraw_individual;
	draw_generations=udraw_generations;
	verbose_level=uverbose_level;
	lookup_bitmap=ulookup_bitmap;
	valid_bitmap=uvalid_bitmap;
	map_filename=umap_filename;
	map_size_x=umap_size_x;
	map_size_y=umap_size_y;
}

//#endif


// genetic params and proc to set them from outside
int nbitx;
int nbity;
int nbitrot;
int popsize;
int maxruns;
int maxgen;
double pcross;
double pmutation;
//---------------------------------------------------------------------------

void setGeneticParameters(int unbitx, int unbity, int unbitrot, int upopsize, int umaxruns, int umaxgen, double upcross, double upmutation)
{   nbitx=unbitx;
nbity=unbity;
nbitrot=unbitrot;
popsize=upopsize;
maxruns=umaxruns;
maxgen=umaxgen;
pcross=upcross;
pmutation=upmutation;
}
//---------------------------------------------------------------------------

unsigned xmask;
unsigned ymask;
unsigned rotmask;
//---------------------------------------------------------------------------

// lookup table, its parameters and proc to set them from outside
unsigned char **lookup;
unsigned char **valid;
double          lookup_step_x;      // lookup_size_x/lookup_columns;
double          lookup_step_y;      // lookup_size_y/lookup_rows;
unsigned char **RGbmp;
unsigned        lookup_rows;        // number of rows= 3000;
unsigned        lookup_columns;     // number of columns= 3000;
double          lookup_size_x;      // dimensione lato x area in metri =10.0
double          lookup_size_y;      // dimensione lato y area in metri =10.0
double          corr_distance;      // distanza entro la quale cercare i corrispondenti
double          lookup_ox;          // coord x dell'angolo inf sin dell'area a ScanRef.pos.x-lookup_size_x/2.0
double          lookup_oy;          // ScanRef.pos.y-lookup_size_y/2.0
unsigned        NX,NY,NX2,NY2,MX,MY; // variabili aux per indicare i limiti dell'area da crescere intorno ai punti di sref
//---------------------------------------------------------------------------

void setLookupTableParameters(int ulookup_rows, int ulookup_columns, double ulookup_size_x, double ulookup_size_y,  double ulookup_centerx, double ulookup_centery, double ucorr_distance)
{   lookup_columns=ulookup_columns;
lookup_rows=ulookup_rows;
lookup_size_x=ulookup_size_x;
lookup_size_y=ulookup_size_y;
lookup_step_x=lookup_size_x/lookup_columns;
lookup_step_y=lookup_size_y/lookup_rows;

// todo: check if ulookup_centerx and y can always be 0. If yes remove them as parameters
lookup_ox=ulookup_centerx-lookup_size_x/2.0;  // coord x lower left corner
lookup_oy=ulookup_centery-lookup_size_y/2.0;  // coord y lower left corner
corr_distance=ucorr_distance;
}

// gray code lookup tables and their initialization (only invGray is used)
const unsigned GraySIZE = 1024;    // todo: attention! max 10 bit per x,y or fi
unsigned Gray[GraySIZE];
unsigned invGray[GraySIZE];

void initialize_Gray_table()
{    for (unsigned i=0; i < GraySIZE; i++)
	Gray[i]=i ^ (i >> 1);
}
void initialize_invGray_table()
{    for (unsigned i=0; i < GraySIZE; i++)
{    int ish;
	unsigned long ans,idiv;
	ish=1; // This is the more complicated direction: In hierarchical stages, starting with a one-bit right shift, cause each bit to be XORed with all more significant bits.
	ans=i;
	for (;;)
	{    ans ^= (idiv=ans >> ish);
	if (idiv <= 1 || ish == 16) break;
	ish <<= 1; // Double the amount of shift on the next cycle.
	}
	invGray[i]=ans;
}
}

// alternative way ...
inline unsigned long unGray(unsigned long r) {
r ^= r>>1;
r ^= r>>2;
r ^= r>>4;
r ^= r>>8;
r ^= r>>16;
return r;
}


void initialize_invGray_table2()
{   for (unsigned i=0; i < GraySIZE; i++)
{   invGray[i]=unGray(i);
}
}

//---------------------------------------------------------------------------
// scans and proc to set them from outside
scan ScanRef;
scan ScanNew;
void setRefScan(const scan& usref)
{
ScanRef.pos=usref.pos;
ScanRef.copy(usref);
}

void setNewScan(const scan& usnew)
{
ScanNew.pos=usnew.pos;
ScanNew.copy(usnew);
}


void add_dither_to_scan(scan &S)
{
const double LSIGMA_L = 0.03;
const double LSIGMA_THETA = 0.003;
const int NumDitherPoints = 5;
// [todo] put these constants as parameters
// [todo] uniformare con radial gradient -> radial gradient specificarlo in coord polari come qui e dipendente dalla distanza

SIT i2=S.readings.begin();

while ( i2 != S.readings.end() )
{
	for(int i=0;i<NumDitherPoints;i++) {
	double eps_l=gauss(LSIGMA_L,0.0);    //[todo] si potrebbe agg. dipendenza dalla distanze dell'errore
	double eps_theta=gauss(LSIGMA_THETA,0.0);
	reading p = {i2->angle+eps_theta,i2->distance+eps_l};
	S.readings.insert(i2,p);
	}
	i2++;
}

}


//---------------------------------------------------------------------------
// search area and proc to set it from outside
// todo: setcenterofsearcharea (search area distance, search area rotation)
double MaxX;
double MinX;
double MaxY;
double MinY;
double MaxFI;
double MinFI;;
double StepX;   // one step through the search space will be derived from the resolution (ie nbitx) and search space size
double StepY;
double StepFI;

void setSearchArea(double uMinX, double uMaxX, double uMinY, double uMaxY,double uMinRot, double uMaxRot)
{   MaxX=uMaxX;
MinX=uMinX;
MaxY=uMaxY;
MinY=uMinY;
if (uMaxRot<uMinRot) {
	MinFI=uMaxRot;
	MaxFI=uMinRot;
} else {
	MaxFI=uMaxRot;
	MinFI=uMinRot;
}
// temp
//MaxFI=2*M_PI;
//MinFI=0.0;
// endtemp
}


void setCenterOfSearchArea(position pos)
{   // todo
}

#ifdef DRAW_PNG

void draw_lookup_table(const char * filename, int type)
{   static bool first_occurance = true; // draw only for first scan pair. currently it is
	if (first_occurance)                //  of little use to draw for every scan pair
	{
	pngwriter img(lookup_columns,lookup_rows,1.0,filename);
	for (unsigned i = 0; i < lookup_columns; i++)
		for (unsigned j = 0; j < lookup_rows; j++)
		{
		if (type==0) img.plot(i,j,(int)lookup[i][j]*50000,(int)lookup[i][j]*50000,(int)lookup[i][j]*50000);
		else if (type==1) img.plot(i,j,(int)lookup[i][j]*256,(int)lookup[i][j]*256,(int)lookup[i][j]*256);
		}
	img.close();
	}
}

void draw_valid_table(const char * filename)
{   static bool first_occurance = true; // draw only for first scan pair. currently it is
	if (first_occurance)                //  of little use to draw for every scan pair
	{
	pngwriter img((1<<nbitx),(1<<nbity),1.0,filename);
	for (unsigned i = 0; i < (1<<nbitx); i++)
		for (unsigned j = 0; j < (1<<nbity); j++)
		{   img.plot(i,j,(int)valid[i][j]*50000,(int)valid[i][j]*50000,(int)valid[i][j]*50000);
		}
	img.close();
	}
}
#endif


void initialize_binary_lookup_table()
{
try
{    // allocazione memoria per array (tabella di lookup)

	lookup= new unsigned char*[lookup_columns]; // sizeX (sizeX puntatori a colonne)

	for ( unsigned j = 0 ; j < lookup_columns ; j++ )
	lookup[j]= new unsigned char[lookup_rows]; // sizeY (sizeY elementi di una colonna)
}
catch (std::bad_alloc &e)
{
	throw std::string("Bad allocation, probabilmente non ho sufficente memoria");
}


// inizializzazione array di lookup
for (unsigned i = 0; i < lookup_columns; i++)
	for (unsigned j = 0; j < lookup_rows; j++)
	lookup[i][j] = 0;

lookup_step_x=lookup_size_x/lookup_columns;
lookup_step_y=lookup_size_y/lookup_rows;

// calcola limiti discreti area da crescere (num celle per avere un area che circoscrive corr_distance_trashold) e variabili ausiliarie
NX= (unsigned int)(corr_distance/lookup_step_x);
NY= (unsigned int)(corr_distance/lookup_step_y);
NX2= NX/2;
NY2= NY/2;
MX= (unsigned int)(lookup_size_x/lookup_step_x-NX2); // maximum X, ie last X cell where border starts
MY= (unsigned int)(lookup_size_y/lookup_step_y-NY2); // maximum Y, ie last Y cell


// cresci
SIT i1=ScanRef.readings.begin();

while ( i1!=ScanRef.readings.end() )
{
	// coordinate X,Y gobali del reading di S1
	double X=ScanRef.pos.x+i1->distance*cos(ScanRef.pos.rot+i1->angle);
	double Y=ScanRef.pos.y+i1->distance*sin(ScanRef.pos.rot+i1->angle);


	// discretizzazione
	unsigned icx=(unsigned int)(floor((X-lookup_ox)/lookup_step_x+0.5)); // discretizzo X+lookup_ox cioe X espresso in coord locali nel frame
	unsigned icy=(unsigned int)(floor((Y-lookup_oy)/lookup_step_y+0.5)); // -

	// praparazione variabili ausiliarie
	unsigned ix=(icx>NX2?icx-NX2:0);// ix e' la cella iniziale, fx la cella finale. IF nelle parentesi per stare attenti ai bordi
	unsigned fx=(icx<MX?icx+NX2:MX+NX2);
	unsigned iy=(icy>NY2?icy-NY2:0);
	unsigned fy=(icy<MY?icy+NY2:MY+NY2);

	if (fx>=lookup_rows) fx=lookup_rows-1;
	if (ix<0) ix=0;
	if (fy>=lookup_columns)    fy=lookup_columns-1;
	if (iy<0) iy=0;


	// crescita nell'intorno di reading di S1
	for ( unsigned i= ix ; i < fx ; i++ )

		for ( unsigned j= iy ; j < fy ; j++ )
		{
		lookup[i][j]=1;
		}
	i1++;
}

#ifdef DRAW_PNG
static bool first=true;
if (draw_lookup && first) {
	draw_lookup_table("./images/lookup_binary.png",0);
	std::cout << "./images/lookup_binary.png" << std::endl;
	first=false;
}
#endif

}
//---------------------------------------------------------------------------

// immagine bitmap e' convertita in tabella di lookup proiettando pixel neri in valori TRUE dell'array
// Purpose of LOOKUP TABLE: x,y real values are discretized, ie mapped into one pixel of LOOKUP TABLE to find out quickly if the position is close to reference scan or map.
void initialize_lookup_table_from_bitmap(std::string pBitmapFileName, // nomefile immagine per la fitness lookup
					double BITMAP_SIZE_X,        // dimensione lato x in metri
					double BITMAP_SIZE_Y,        // dimensione lato y in metri                                            double BORDER_SIZE)          // dimensione bordo che viene aggiunto alla bitmap (lookup che sia un po piu grande)
					double BORDER_SIZE)
{

#ifdef DRAW_PNG
pngwriter bmp_lookup;

bmp_lookup.readfromfile(pBitmapFileName.c_str());
unsigned BITMAP_ROWS=bmp_lookup.getheight();
unsigned BITMAP_COLUMNS=bmp_lookup.getwidth();
// deduco tutti i parametri della lookup dalla risoluzione della bitmap e fisso l'origine del sist. di coord in mezzo alla lookup
unsigned BORDER_COLUMNS=BITMAP_COLUMNS/BITMAP_SIZE_X*BORDER_SIZE;
unsigned BORDER_ROWS=BITMAP_ROWS/BITMAP_SIZE_Y*BORDER_SIZE;
lookup_columns=BITMAP_COLUMNS+BORDER_COLUMNS*2;
lookup_rows=BITMAP_ROWS+BORDER_ROWS*2;
lookup_size_x=BITMAP_SIZE_X+BORDER_SIZE*2.0;
lookup_size_y=BITMAP_SIZE_Y+BORDER_SIZE*2.0;
lookup_step_x=lookup_size_x/lookup_columns;
lookup_step_y=lookup_size_y/lookup_rows;
lookup_ox=-lookup_size_x/2.0;  // coord x lower left corner (ATTENZIONE: ORIGINE SIST. COORD E' SEMPRE IN MEZZO ALLA LOOKUP
lookup_oy=-lookup_size_y/2.0;  // coord y lower left corner   -//-
//std::cout << "PNG READ GLASM " << lookup_rows << "R: " << lookup_columns << std::endl;


// TEMP verifica - plot su schermo come verrebbe letto
// coordinata (x,y)=(0,0) e' in angolo inferiore sinistro dell'immagine .png cosi come vista in gimp
printf("\nBITMAP: %s",pBitmapFileName.c_str());
printf("\nrows: %d",BITMAP_ROWS);
printf("\ncolumns: %d",BITMAP_COLUMNS);
printf("\nborder columns: %d",BORDER_COLUMNS);
printf("\nborder rows: %d",BORDER_ROWS);
printf("\nlookup_size_x: %f",lookup_size_x);
printf(" lookup_size_y: %f",lookup_size_y);
printf(" lookup_rows: %d",lookup_rows);
printf(" lookup_columns: %d",lookup_columns);
printf(" lookup_ox: %f",lookup_ox);
printf(" lookup_oy: %f\n",lookup_oy);
printf(" BITMAP_SIZE_X: %f",BITMAP_SIZE_X);
printf(" BITMAP_SIZE_Y: %f\n",BITMAP_SIZE_Y);

// allocazione memoria per array (tabella di lookup)
try
{   lookup=new unsigned char*[lookup_columns];
	for(unsigned j=0;j<lookup_columns;j++)
	lookup[j]=new unsigned char[lookup_rows];
}
catch (std::bad_alloc &e)
{    throw std::string("Bad allocation, probabilmente non ho sufficente memoria");
}

// inizializzazione array di lookup
for (unsigned y = 0; y < lookup_rows; y++)
{  for (unsigned x = 0; x < lookup_columns; x++)
{   if (x<BORDER_COLUMNS || x>=lookup_columns-BORDER_COLUMNS) lookup[x][y] = 0; else    // bordi sono sempre false
	if (y<BORDER_ROWS    || y>=lookup_rows-BORDER_ROWS)       lookup[x][y] = 0; else    // bordi sono sempre false
	if (bmp_lookup.read(x+1,y+1)<100)                         lookup[x][y] = 1; else    // true quando il pixel e' verso il colore nero
	lookup[x][y] = 0;
}
}


static bool first=true;
	if (draw_lookup && first ) {
	draw_lookup_table("./images/lookup_binary_frombitmap.png",0);
	first=false;
	}

#endif

}
//---------------------------------------------------------------------------

// tabella per indicare posizioni valide da quelle non valide (masks valid positions for search area)
// Purpose of VALID_POS_LOOKUP TABLE: discrete x,y chromosome positions (ie 8 bit for x => 256 positions in x direction) are
//  mapped as valid or not based on existing bitmap of valid positions. VALID_POS_LOOKUP table is nothing more then original bitmap of valid positions squeezedfitted/covered in the resolution.
void initialize_valid_table_from_bitmap(std::string pBitmapFileName)
{

#ifdef DRAW_PNG
pngwriter bmp_valid;
bmp_valid.readfromfile(pBitmapFileName.c_str());
unsigned valid_rows=bmp_valid.getheight();
unsigned valid_columns=bmp_valid.getwidth();
unsigned rows=(1<<nbity);
unsigned columns=(1<<nbitx);
double step_x=valid_columns/double(columns);
double step_y=valid_rows/double(rows);


// TEMP verifica - plot su schermo come verrebbe letto
// coordinata (x,y)=(0,0) e' in angolo inferiore sinistro dell'immagine .png cosi come vista in gimp
printf("\nValid pos BITMAP  : %s",pBitmapFileName.c_str());
printf("\nbitmap rows       : %d",valid_rows);
printf("\nbitmap columns    : %d",valid_columns);
printf("\nrows (1<<nbity)   : %d",rows);
printf("\ncolumns (1<<nbitx): %d",columns);
printf("\nstep_x             : %f",step_x);
printf("\nstep_y             : %f",step_y);

// allocazione memoria per array (tabella di valid)
try
{   valid=new unsigned char*[columns];
	for(unsigned j=0;j<columns;j++)
	valid[j]=new unsigned char[rows];
}
catch (std::bad_alloc &e)
{    throw std::string("Bad allocation, probabilmente non ho sufficente memoria");
}

// initialization of table of valid positions
for (unsigned y = 0; y < rows; y++)
{  for (unsigned x = 0; x < columns; x++)
{    if (bmp_valid.read(x*step_x,y*step_y)<170) valid[x][y] = 1; else valid[x][y]=0;
}
}



if (draw_lookup)
{   draw_valid_table("./images/valid_pos_lookup.png");
}
#endif

}
//---------------------------------------------------------------------------



void initialize_radial_gradient(double RGmean, double RGvariance, double RGextension, double RGweight, double RGscale, double RGbitmapsize)
{   // prepare radial gradient bitmap
//  this will be impressed in lookup on top of each sref reading
// double RGmean - normal curve mean
// double RGvariance - normal curve variance
// double RGweight - weight given to normalized normal curve when drawing - ie to increase 'color intensity'
// double RGscale - pixels in one unit distance (meter)
// double RGbitmapsize - size in distance units (meters)


// bitmap size in pixels relative to lookup table size
NX = (unsigned int)(RGbitmapsize*RGscale);
NX2= NX/2;
NY=NX;
NY2=NX2;
//std::cout << "bmp" << RGmean << "," << RGvariance << "," << RGextension << "," << RGweight << "," << RGscale << "," <<  RGbitmapsize << std::endl;
try
{   // allocazione memoria per array (RadialGradientBitmap)

	RGbmp= new unsigned char*[NX]; // righe dell'array

	for ( unsigned j = 0 ; j < NX ; j++ )
	RGbmp[j]= new unsigned char[NX]; // colonne dell'array
}
catch (std::bad_alloc  &e)
{
	throw std::string("Bad allocation, probabilmente non ho sufficente memoria");
}


// inizializzazione array di lookup
for (unsigned i = 0; i < NX; i++)
	for (unsigned j = 0; j < NX; j++)
	{      double d=sqrt((i-NX2)*(i-NX2)+(j-NX2)*(j-NX2))/RGextension;
	double t1 = (d -  RGmean)/RGvariance;
	double t2 = -0.5*t1*t1;
	RGbmp[i][j] = (unsigned) (exp(t2)/(RGvariance*sqrt(2*M_PI))*RGweight);
	//if (RGbmp[i][j]>127) RGbmp[i][j]=127; // clipping
	//if (RGbmp[i][j]<0) RGbmp[i][j]=0; // clipping
	}
// temp: checking ...
//    for (unsigned j = 0; j < NX; j++)
//    {   std::cout << "  R=" << (unsigned) RGbmp[NX2][j] << std::endl;
//    }
// end temp

#ifdef DRAW_PNG
if (draw_lookup)
{   static bool first_occurance = true; // draw only for first scan pair. currently it is
	if (first_occurance)                //  of little use to draw for every scan pair
	{   char filename[FILENAME_SIZE];   //  in experiments with multiple scan pairs
	snprintf(filename,FILENAME_SIZE,"./images/RGbmp.png");
	pngwriter img(NX,NX,1.0,filename);
	for (unsigned i = 0; i < NX; i++)
		for (unsigned j = 0; j < NX; j++)
		{   //unsigned tmp = RGbmp[i][j];
		img.plot(i,j,(int)RGbmp[i][j]*256,(int)RGbmp[i][j]*256,(int)RGbmp[i][j]*256);
//                img.plot(i,j,(int)i*200,(int)i*200,(int)j*200);
		}
	img.close();
	first_occurance=false;
	}
}
#endif
} // initialize_radial_gradient

void initialize_gradient_lookup_table()
{
try
{   // allocazione memoria per array (tabella di lookup)

	lookup= new unsigned char*[lookup_columns]; // righe dell'array

	for ( unsigned j = 0 ; j < lookup_columns ; j++ )
	lookup[j]= new unsigned char[lookup_rows]; // colonne dell'array
}
catch (std::bad_alloc &e)
{
	throw std::string("Bad allocation, probabilmente non ho sufficente memoria");
}


// inizializzazione array di lookup
for (unsigned i = 0; i < lookup_columns; i++)
	for (unsigned j = 0; j < lookup_rows; j++)
	lookup[i][j] = 0;

unsigned MX= (unsigned int)(lookup_size_x/lookup_step_x-NX2); // maximum X, ie last X cell where border starts
unsigned MY= (unsigned int)(lookup_size_y/lookup_step_y-NY2); // maximum Y, ie last Y cell


// cresci
SIT i1=ScanRef.readings.begin();

while ( i1!=ScanRef.readings.end() )
{
	// coordinate X,Y gobali del reading di S1
	double X=ScanRef.pos.x+i1->distance*cos(ScanRef.pos.rot+i1->angle);
	double Y=ScanRef.pos.y+i1->distance*sin(ScanRef.pos.rot+i1->angle);


	// discretizzazione
	unsigned icx=(unsigned int)(floor((X-lookup_ox)/lookup_step_x+0.5)); // discretizzo X+lookup_ox cioe X espresso in coord locali nel frame
	unsigned icy=(unsigned int)(floor((Y-lookup_oy)/lookup_step_y+0.5)); // -


	// praparazione variabili ausiliarie
	unsigned ix=(icx>NX2?icx-NX2:0);// ix e' la cella iniziale, fx la cella finale. IF nelle parentesi per stare attenti ai bordi
	unsigned fx=(icx<MX?icx+NX2:MX+NX2);
	unsigned iy=(icy>NY2?icy-NY2:0);
	unsigned fy=(icy<MY?icy+NY2:MY+NY2);

	if (fx>=lookup_rows) fx=lookup_rows-1;
	if (ix<0) ix=0;
	if (fy>=lookup_columns) fy=lookup_columns-1;
	if (iy<0) iy=0;


	// crescita nell'intorno di reading di S1
	unsigned sum;
	for ( unsigned i= ix ; i < fx ; i++ )

		for ( unsigned j= iy ; j < fy ; j++ )
		{    sum=lookup[i][j]+RGbmp[i-ix][j-iy];
		if (sum>255) lookup[i][j]=255; else lookup[i][j]=sum;
		//lookup[i][j]+=RGbmp[i-ix][j-iy];
		}
	i1++;
}

#ifdef DRAW_PNG
static bool first=true;
	if (draw_lookup && first ) {
	draw_lookup_table("./images/lookup_gradient.png",1);
	std::cout << "./images/lookup_gradient.png" << std::endl;
	first=false;
	}

#endif

}
//---------------------------------------------------------------------------


void delete_lookup_table()
{    for (unsigned i = 0; i < lookup_columns; i++)
	delete[] lookup[i]; // delete the columns
delete[] lookup;// delete the rows
}


//---------------------------------------------------------------------------


#ifdef DRAW_PNG

void drawfun(int run, int gen, int popsize, struct individual *newpop, struct bestever *bestfit)
{   static unsigned prevChrom;
	char filename[FILENAME_SIZE];

	if (draw_generations || draw_individual)
	{
	snprintf(filename,FILENAME_SIZE,"./images/GLASM_run_%02d_gen_%02d.png",run,gen);
	//std::cout << filename << " map:" <<map_filename<< " "<< map_size_x << " " << map_size_y << std::endl;
	drawing img(filename,map_size_x,map_size_y,-map_size_x/2.0,-map_size_y/2.0,map_filename.c_str());
	//std::cout << "bbb" << std::endl;
	img.grid();
	img.setpointsize(1);

	// draw population
	for (int ii=0; ii<popsize; ii++)
	{   chromosome2pos(newpop[ii].chrom);
		img.setcolor(10000,10000,10000);
		img.readings(ScanNew, ScanNew.pos);
		if (draw_individual)
		{   snprintf(filename,FILENAME_SIZE,"./images/GLASM_run_%d_gen_%d_ind_%d.png",run,gen,ii);
		//drawing img1(filename,lookup_size_x,lookup_size_y,lookup_ox,lookup_oy,"./cfg/bitmaps/empty_1000x1000.png");
		drawing img1(filename,map_size_x,map_size_y,-map_size_x/2.0,-map_size_y/2.0,map_filename.c_str());
		img1.grid();
		// draw lookup
		img1.setcolor(0,50000,0);
		img1.setpointsize(10);
		for (unsigned i = 0; i < lookup_columns; i++)
			for (unsigned j = 0; j < lookup_rows; j++)
			{   if (lookup[i][j]) img1.point(lookup_step_x*i,lookup_step_y*j);
			}

		// draw readings
		img1.setcolor(10000,10000,10000);
		img1.readings(ScanNew, ScanNew.pos);

		// draw fitness value
		const unsigned NAMESIZE=255;
		char tmp[NAMESIZE];
		memset(tmp,0,NAMESIZE);
		strncpy(tmp,"(pose) fitness: ",16);
		sprintf(&tmp[16],"(%f,%f,%f) %f",ScanNew.pos.x, ScanNew.pos.y, ScanNew.pos.rot, newpop[ii].fitness);
		img1.setcolor(10,20000,10);
		img1.text(10,10,tmp);
		img1.write();
		} // draw_individual
	} // for

	if (draw_generations)
	{   if (bestfit->chrom!=prevChrom)
		{   chromosome2pos(prevChrom);
		img.setcolor(60000,10000,60000);
		img.readings(ScanNew, ScanNew.pos);
		prevChrom=bestfit->chrom;
		}
		chromosome2pos(bestfit->chrom);
		img.setcolor(60000,0,20000);
		img.readings(ScanNew, ScanNew.pos);
		img.write();
	} // draw_generations
	}
}
#endif


// todo: use this function for quick float2int in objfun. Should be faster ...
inline int float2int( double d )
{
	union Cast
	{
		double d;
		long l;
		};
	volatile Cast c;
	c.d = d + 6755399441055744.0;
	return c.l;
}



void objfun(struct individual *critter)
{
	// if (cachehit(critter->chrom)


if (!chromosome2pos(critter->chrom))
{   critter->fitness = 0;                   // invalid position: for now if position is invalid put fitness at 0 (ie low)
	return;
}
unsigned counter=0;

SIT i2=ScanNew.readings.begin();
//    unsigned outlier=0;

while ( i2 != ScanNew.readings.end() )
{
	// coordinate X,Y gobali del reading di S2
	double X=ScanNew.pos.x+i2->distance*cos(ScanNew.pos.rot+i2->angle);
	double Y=ScanNew.pos.y+i2->distance*sin(ScanNew.pos.rot+i2->angle);

	// discretizzazione
	unsigned ix=(unsigned int)(floor( (X - lookup_ox) / lookup_step_x + 0.5) ); // discretizzo X+lookup_ox cioe X espresso in coord locali nel frame
	unsigned iy=(unsigned int)(floor( (Y - lookup_oy) / lookup_step_y + 0.5) ); // -

	if (ix<lookup_columns && ix>=0 && iy<lookup_rows && iy>=0)
	{
	counter+=lookup[ix][iy];     // lookup

//            // lookup, guardiamo se vicino cera un corrispondente
//            if ( lookup[ix][iy] ) counter++;
//            else if (outlier++ > outlier_treshold)
//            {   critter->fitness = 0;                   // too many readings in invalid positions: for now put fitness at 0 (ie low)
//                return;
//            }
	}

	i2++;
}
critter->fitness = counter;
}
//---------------------------------------------------------------------------

// TODO: avoid function call (make it void as before or define inline)
bool chromosome2pos(const unsigned c)// interpreta posizione dal cromosoma (sapendo che il chromosoma e' < di 32 bit (ci sta in un unsigned di 32 bit) ho semplificato un po')
{
unsigned x,y,fi;

x=(c&xmask);// primi nbitx bit
x=invGray[x];// codice Gray inverso

y=c;
y=y>>nbitx;
y=(y&ymask);
y=invGray[y];


// check if chromosome is valid (resulting position, only x and y, is inside search area)
// we do not decode cromosome to pos. It is much faster to perform a lookup in valid positions mask
//    if (valid[x][y])
//    {
//        // temp for debug
////        std::cout << "  invalid pos in x: " << x << " y: " << y << "  real X: " << MinX+x*StepX << " Y: " << MinY+y*StepY << std::endl;
////        img_temp.plot(x,y,50000,0,0);
//        // end temp
//        return false;
//    }



fi=c;
fi=fi>>(nbitx+nbity);
fi=(fi&rotmask);
fi=invGray[fi];

ScanNew.pos.x=MinX+x*StepX;
ScanNew.pos.y=MinY+y*StepY;
ScanNew.pos.rot=MinFI+fi*StepFI;
if (ScanNew.pos.rot>2*M_PI) ScanNew.pos.rot-=2*M_PI; // this line should not be necessary ???
return true;

}
//---------------------------------------------------------------------------

position glasm(double* outfitness)
// input parameters (must be set prior to calling GLASM):
//    reference and new scan (the provided position of the new scan is not used and can be whatever, it will be estimated by the algorithm.
//        this is different in local scan matching algorithms where this position is used as the initial position estimate and is only corrected by the algorithm)
//    search area
//     GLASM genetic parameters (nbitx,nbity,nbitrot,popsize,maxruns,maxgen,pcross,pmutation)
//    GLASM lookup table parameters
{
sga_parameters(nbitx,nbity,nbitrot,popsize,maxruns,maxgen,pcross,pmutation);    // set genetic parameters

xmask=(1<<nbitx)-1;
ymask=(1<<nbity)-1;
rotmask=(1<<nbitrot)-1;

unsigned ResX=1<<nbitx;
StepX=(MaxX-MinX)/ResX;

unsigned ResY=1<<nbity;
StepY=(MaxY-MinY)/ResY;
unsigned ResFI=1<<nbitrot;

angleUnwrap(MinFI);
angleUnwrap(MaxFI);
if (MinFI>MaxFI)
{   // swap
	double t=MaxFI;
	MaxFI=MinFI;
	MinFI=t;
}
double anglediff=MaxFI-MinFI;
StepFI=anglediff/ResFI;
set_sga_objfun(&objfun);

#ifdef DRAW_PNG
set_sga_drawfun(&drawfun);
#endif

bestever bestfit=sga();
  *outfitness = bestfit.fitness;
chromosome2pos(bestfit.chrom);
return ScanNew.pos;
} // objfun

} // namespace
