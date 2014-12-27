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

//---------------------------------------------------------------------------
// This unit contains all the utility and accessory functions not specific to
// algorithms or other units
//
//---------------------------------------------------------------------------

#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

const double SOGLIA_ANGOLO_OUTLIERS     =1.1/180.0*M_PI;
const double SOGLIA_DISTANZA_OUTLIERS   =2.5;
const double CORR_ANGLE_TRESHOLD        =8.5/180.0*M_PI;

double SENSOR_RANGE_DEF=30.0;
//double SENSOR_RANGE_DEF=8.0;

const double SICK_STEP_ANGLE=M_PI/180; // check this value for sick laser (should be 1.0 degree or 0.9 degrees)





void scan::read(const char *filename)
{   FILE *in;
	float f1,f2,f3;
	if ((in = fopen(filename, "r")) == NULL)
	{   std::cout<<"readscan: Error opening file"<< std::endl;
	}

	if (fscanf(in,"scan.pos:%f,%f,%f\n",&f1,&f2,&f3)<3)
	{   std::cout<<"readscan: Error reading scan position"<< std::endl;
	} else
	{   pos.x=f1;
		pos.y=f2;
		pos.rot=f3;
		readings.clear();
		unsigned i=0;
		while (fscanf(in,"%f,%f\n",&f1, &f2)==2)
		{   reading r;
			r.distance=f1;
			r.angle=f2;
			if (r.distance<7.8) readings.push_back(r);
			i++;
		}
	}
	if (fclose(in) !=  0)
	{   std::cout<<"readscan: Error closing file"<< std::endl;
	}
}

void scan::read(std::string &rawseeds_scans)
{	readings.clear();
// use a stringstream to separate the fields out of the line
istringstream ss( rawseeds_scans );
string field;
unsigned i=0;
double angle=0.0;
while (getline( ss, field, ',' ))
{       // for each field we wish to convert it to a double
	stringstream fs( field );
	double f = 0.0;  // (default value is 0.0)
	fs >> f;

	switch (i)
	{       case 0: break; // ignore timestamp
		case 1: pos.x=f; break;
		case 2: pos.y=f; break;
		case 3: pos.rot=f; break;
		case 4: break; // ignore numreadings
		case 5: break; // ignore this field (what is it?)
		default:
		{       reading r;
			r.distance=f;
			r.angle=angle;
			angle+=SICK_STEP_ANGLE;
			readings.push_back( r ); // add the newly-converted field to the end of the record
		}
	}
	i++;
}
}

void scan::write(const char *filename)
{   FILE *out;
	if ((out = fopen(filename, "w")) == NULL)
	{   std::cout<<"writescan: Error opening file"<< std::endl;
		return;
	}
	if (fprintf(out,"scan.pos:%f,%f,%f\n",pos.x,pos.y,pos.rot)<3)
	{   std::cout<<"writescan: Error writing scan position"<< std::endl;
	} else
	{   CSIT si=readings.begin();
		while (si!=readings.end())
		{   fprintf(out,"%f,%f\n",si->distance,si->angle);
			si++;
		}
	}
	if (fclose(out) !=  0)
	{   std::cout<<"writescan: Error closing file"<< std::endl;
	}
}








double maxRange(void)
{

	return SENSOR_RANGE_DEF;

}

void   setMaxRange(double pMRange)
{

	SENSOR_RANGE_DEF=pMRange;

}

void assign_pos(position &posTo, const position &posFrom)
{
	posTo.x=posFrom.x;
	posTo.y=posFrom.y;
	posTo.rot=posFrom.rot;
}


void calc_xy(double newPosRot, double newPosRan, position& pPos)
{

	double newBearing(0);

	newBearing=newPosRot;

	pPos.rot = pPos.rot+newBearing;

	if (pPos.rot>=M_PI)

		pPos.rot=pPos.rot-2*M_PI;

	else if (pPos.rot<=-M_PI)

		pPos.rot=pPos.rot+2*M_PI;

	pPos.x = pPos.x + newPosRan*cos(pPos.rot);
	pPos.y = pPos.y + newPosRan*sin(pPos.rot);

}

void calc_RanRot(position pPosOn, position pPosTo, double& newPosRot, double& newPosRan)
{

	//double  angle,distance;
	double  dx,dy; //,drot;

	dx=pPosTo.x-pPosOn.x;
	dy=pPosTo.y-pPosOn.y;
	//drot=pPosTo.rot-pPosOn.rot;

	/*std::cout << "\t  XYZ :" << pPosTo.x  << " Y :" << pPosTo.y << " X :" << pPosTo.rot << std::endl;
	std::cout << "\t  XYZ :" << pPosOn.x  << " Y :" << pPosOn.y << " X :" << pPosOn.rot << std::endl;
	std::cout << "\t  XYZ :" << dx  << " Y :" << dy << " X :" << drot << std::endl; */

	newPosRan=sqrt(dx*dx+dy*dy);
	//std::cout << "\t  DIS :" << newPosRan << std::endl;
	/*dx=pPosTo.x-pPosOn.x;
	dy=pPosTo.y-pPosOn.y;*/
	//if (dy!=0)
		newPosRot=atan2(dy,dx);
	//else
	//  angle=0;
//std::cout << "RR :" <<  newPosRot << " DIS :" << newPosRan  << std::endl;


	//newPosRan=distance;
	//newPosRot=angle;


}

//---------------------------------------------------------------------------
bool readingsAngleCompare(reading r1, reading r2)
{	if (r1.angle<r2.angle) return true;
	return false;
}
//---------------------------------------------------------------------------

bool compareByS1Angle(corrPoints c1, corrPoints c2)
{   if (c1.pref->angle<c2.pref->angle) return true;
return false;
}
//---------------------------------------------------------------------------

void locToGlob(position pos, double lx, double ly, double& gx, double& gy)
// pos (x,y,rot) (ingresso, noti)               sono il centro e la rotazione del sistema locale rispetto a quello globale
// lx,ly    (ingresso, noti)                    sono le coordinate del punto nel sistema locale (x,y,alfa)
// gx,gy    (uscita, calcolati dalla funzione)  sono le coordinate del punto nel sistema globale

// per ottimizzare il calcolo nelle liste con tanti punti consiglio di evitare chiamate di funzione e
// implementare la trasformazione direttamente copiando le righe di codice
{
gx=pos.x+lx*cos(pos.rot)-ly*sin(pos.rot);
gy=pos.y+lx*sin(pos.rot)+ly*cos(pos.rot);
}
//---------------------------------------------------------------------------

void globToLoc(position pos, double gx, double gy,double& lx,double& ly)
// pos (x,y,rot) (ingresso, noti)              sono il centro e la rotazione del sistema locale rispetto a quello globale
// gx,gy    (ingresso, noti)                    sono le coordinate del punto nel sistema globale
// lx,ly    (uscita, calcolati dalla funzione)  sono le coordinate del punto nel sistema locale (x,y,rot)

// per ottimizzare il calcolo nelle liste con tanti punti consiglio di evitare chiamate di funzione e
// implementare la trasformazione direttamente copiando le righe di codice
{   double tx=gx-pos.x;
double ty=gy-pos.y;
double ag,al,c;
ag=pointAngle(tx,ty);
al=angleUnwrap(ag-pos.rot);    // non serve
c=sqrt(tx*tx+ty*ty);
lx=cos(al)*c;
ly=sin(al)*c;
}
//---------------------------------------------------------------------------

double pointAngle(double x, double y)
{   double a;
if (x==0)
	(y>=0 ? a=M_PI/2 : a=-M_PI/2);  // nel caso x=0, y=0 restituisce PI mezzi, ma qualunque angolo va bene
else
{   a=atan(y/x);
	if (x<0) a+=M_PI;
	else if (y<0) a+=2*M_PI;
}
return a;
}
//---------------------------------------------------------------------------

bool CCW(point a, point b, point c)
// returns true if the points a,b and c are in counter clockwise direction
{
return (c.y-a.y)*(b.x-a.x)>(b.y-a.y)*(c.x-a.x);
}
//---------------------------------------------------------------------------

bool intersect(point a, point b, point c, point d)
// returns true if the segment ab intersects the segment cd
{
return (CCW(a,c,d)!=CCW(b,c,d))&&(CCW(a,b,c)!=CCW(a,b,d));
}
//---------------------------------------------------------------------------

double angleInternalDifference(double arc1, double arc2)
{   if (arc1>arc2 && arc1-arc2>=M_PI) arc1-=2*M_PI; else
if (arc2>arc1 && arc2-arc1>=M_PI) arc2-=2*M_PI;
return fabs(arc1-arc2);
}
//---------------------------------------------------------------------------

bool angleInArc(double alfa, double arc1, double arc2)
{   alfa=angleUnwrap(alfa);
arc1=angleUnwrap(arc1);
arc2=angleUnwrap(arc2);
bool interno;
if (arc1>arc2) { double temp=arc1; arc1=arc2; arc2=temp; }      // swap
if (alfa>arc1 && alfa<arc2) interno=true; else interno=false;   // fin qua funzionerebbe bene su un segmento
if (arc2-arc1>M_PI) interno=!interno;                           // tiene conto del wrapping della fase
return interno;
}
//---------------------------------------------------------------------------

double angleUnwrap(double alfa)
{  //if (alfa>2*M_PI) return alfa-2*M_PI;
//if (alfa<0) return alfa+2*M_PI;
//return alfa;

return alfa-floor(alfa/2.0/M_PI)*2*M_PI;
}
//---------------------------------------------------------------------------

double gauss(double var, double med)
{   double result=0.0;
for (int i=0;i<26;i++)
	result+=_random(1000)/1000.0;
result/=26.0;
return (result-0.5)*sqrt(var)+med;
}
//---------------------------------------------------------------------------

double randDouble(double low, double high)
{	double temp;
/* swap low & high around if the user makes no sense */
if (low > high)
{
	temp = low;
	low = high;
	high = temp;
}

/* calculate the random number & return it */
temp = (rand() / (static_cast<double>(RAND_MAX) + 1.0))	* (high - low) + low;
return temp;
}
//---------------------------------------------------------------------------

double distPolar(double a1, double d1, double a2, double d2)
{   double t=d1*d1+d2*d2-2*d1*d2*cos(a1-a2);
if (t<0) t=0;
return sqrt(t);
}
//---------------------------------------------------------------------------

int _random(unsigned lim)
{	if (lim>0) return rand()%(lim);
return 0;
}
//---------------------------------------------------------------------------

void scanChangeRefFrame(scan& Sref, position newPos)
{   SIT LI = Sref.readings.begin();
while (LI!=Sref.readings.end())
{   // invece di usare due passaggi *usando x,y invece di angolo e dist direttamente) farlo direttamente!!
	double x=cos(LI->angle)*LI->distance;                           // attuali (old x e y)
	double y=sin(LI->angle)*LI->distance;
	double xg,yg;                                                   // globali
	locToGlob(Sref.pos,x,y,xg,yg);                                  // spostiamo il punto prima nelle coordinate globali (si trovano: xg,yg)
	globToLoc(newPos,xg,yg,x,y);                                    //   e poi da globali nel nuovo sist rif locale (si trovano nuovi: x, y)
	LI->angle=pointAngle(x,y);                                     // aggiorniamo l'angolo
	LI->distance=sqrt(x*x+y*y);
	LI++;
}
Sref.pos=newPos;
}
//---------------------------------------------------------------------------

void readingsChangeRefFrame(std::list<reading>& readings, position oldPos, position newPos)
{   SIT LI = readings.begin();
while (LI!=readings.end())
	{   // invece di usare due passaggi *usando x,y invece di angolo e dist direttamente) farlo direttamente!!
		double x=cos(LI->angle)*LI->distance;                           // attuali (old x e y)
		double y=sin(LI->angle)*LI->distance;
		double xg,yg;                                                   // globali
		locToGlob(oldPos,x,y,xg,yg);                                  // spostiamo il punto prima nelle coordinate globali (si trovano: xg,yg)
		globToLoc(newPos,xg,yg,x,y);                                    //   e poi da globali nel nuovo sist rif locale (si trovano nuovi: x, y)
		LI->angle=pointAngle(x,y);                                     // aggiorniamo l'angolo
		LI->distance=sqrt(x*x+y*y);
		LI++;
	}
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


void scanMedianFilter(scan& S, const unsigned N, const double DIST)
// S    - scansione le cui letture saranno trattate con filtro mediano
// N    - numero di elementi vicini (ordine del filtro)
// DIST - distanza tra letture al di sopra della quale il filtro non viene applicato
{
const unsigned n = N/2;  // meta di N

std::vector<position> p(N);
unsigned i;

int distanti=0;                 // solo per analisi / si puo togliere dopo
int non_distanti=0;

// NOTA: i punti vicini degli estremi sono come in lista circolare, ovvero l'intorno del punto finale
// comprende i due punti precedenti e poi il punto iniziale e il sucessivo (caso N=5)
SIT LI=S.readings.begin();
while (LI!=S.readings.end())
{   // passo 1: estrai i punti vicini e mettili in una struttura piu comoda (array)
	// il punto corrente
	p[n].x=LI->distance*cos(LI->angle);
	p[n].y=LI->distance*sin(LI->angle);
	SIT I=LI;
	for (i=0; i<n; i++)                         // i punti a sinistra
	{   if (I==S.readings.begin()) I=S.readings.end();
	I--;
	p[n-i-1].x=I->distance*cos(I->angle);
	p[n-i-1].y=I->distance*sin(I->angle);
	}
		I=LI;
		for (i=0; i<n; i++)                         // i punti a destra
		{   I++;
			if (I==S.readings.end()) I=S.readings.begin();
			p[n+i+1].x=I->distance*cos(I->angle);
			p[n+i+1].y=I->distance*sin(I->angle);

		}

		double Xmedio,Ymedio;
		Xmedio=Ymedio=0.0;
		for (i=0; i<N; i++)
		{   Xmedio+=p[i].x;
			Ymedio+=p[i].y;
		}
		Xmedio/=N;
		Ymedio/=N;

		bool too_distant=false;
		for (i=0; i<N; i++)
		{   double x=p[i].x-Xmedio;
			double y=p[i].y-Ymedio;
			if (sqrt(x*x+y*y)>DIST)
			{   too_distant=true;
			}
		}

		if (!too_distant)
	{   LI->angle=pointAngle(Xmedio,Ymedio);
			LI->distance=sqrt(Xmedio*Xmedio+Ymedio*Ymedio);
			non_distanti++;                                                   // solo per analisi
	} else                                                                //
			distanti++;                                                       //
		LI++;
	}
}

//---------------------------------------------------------------------------

void scanResample(scan& S, const double DIST)
// S    - scansione le cui letture saranno decimate
// DIST - distanza tra letture minima al di sotto della quale saranno decimate
{   SIT LI=S.readings.begin();

	while (LI!=S.readings.end())
	{   SIT I=LI; I++; if (I==S.readings.end()) I=S.readings.begin();
		double x1=LI->distance*cos(LI->angle);
		double y1=LI->distance*sin(LI->angle);
		double x2=I->distance*cos(I->angle);
		double y2=I->distance*sin(I->angle);
	double x=x2-x1;
	double y=y2-y1;
	if (sqrt(x*x+y*y)<DIST)
	I=S.readings.erase(I);
	else
	LI++;
}
}

//---------------------------------------------------------------------------

void scanDecimate(scan& S, const unsigned ORDER, const unsigned OFFSET)
// S     - scansione le cui letture saranno decimate
// ORDER - numero di letture da rimuovere per ogni lettura che resta
// OFFSET - posizione della prima lettura che resta (tipicamente 0<OFFSET<ORDER-1)
{   if (ORDER<2) return;			// nothing to do if ORDER < 2 (ORDER=2 means every second reading will be erased)
SIT LI=S.readings.begin();
unsigned i=0;
while (i<OFFSET && LI!=S.readings.end()) { LI=S.readings.erase(LI); i++; }
i=0;
	while (LI!=S.readings.end())
	{   if (i%ORDER) LI=S.readings.erase(LI); else LI++;
	i++;
	}
}

//---------------------------------------------------------------------------

void readingsInAngleRange(scan& S1, scan& S2, double infAngle, double supAngle, SIT& infReading, SIT& supReading)

// ipotesi: all-inizio infReading e supReading non devono puntare a S2.readings.end()
// infReading e supReading sono const_iterator della scansione S2 di cui si vogliono trovare gli elementi
{   SIT END=S2.readings.end();
SIT BEG=S2.readings.begin();

if (infAngle>supAngle)          // siamo a cavallo (wrapping della fase)
{   infReading=END; infReading--;
	while (infReading->angle>infAngle && infReading!=BEG) infReading--;
	infReading++;
	if (infReading!=END)
	if (infReading->angle<infAngle)
		infReading=END;       // non trovato

	supReading=BEG;
	if (supReading->angle>supAngle) supReading=END; // non ci sono elementi in arco
	else
	{   while (supReading->angle<supAngle && supReading!=END) supReading++;
	supReading--;
	}
}
else
{   if (angleInArc(supReading->angle, infAngle, supAngle))
	{   do supReading++;
	while (supReading!=END && supReading->angle<=supAngle);
	supReading--;
	} else
	{   supReading=BEG;             // caso poco frequente che non sia nell'arco percio posso spendere qualche ciclo facendo ricerca dall'inizio
	if (supReading->angle>supAngle) supReading=END; // non ci sono elementi in arco
	else
	{   while (supReading!=END && supReading->angle<=supAngle) supReading++;
		supReading--; if (supReading->angle<infAngle) supReading=END;
	}
	}
	infReading=supReading;             // caso frequente e allora cominciamo dal supREading che e spesso piu vicino
	while (infReading!=BEG && infReading->angle>infAngle) infReading--;
	if (infReading->angle<infAngle) infReading++;
}
}




//---------------------------------------------------------------------------

//-------------------------- OUTLIERS ---------------------------------------
//---------------------------------------------------------------------------



void outliers(scan& S1, scan& S2)
// questa procedura prende in ingresso la SCAN1-1 e SCAN2-2
{   position r1=S1.pos;                                                 // memorizza posizioni
position r2=S2.pos;

scan T=S1;
scanChangeRefFrame(S1, r2);
POTooFar(S1);
S1.readings.sort(readingsAngleCompare);
POHiddenByOther(S1,S2);

scanChangeRefFrame(S2, r1);
POTooFar(S2);
S2.readings.sort(readingsAngleCompare);
POHiddenByOther(S2,T);

scanChangeRefFrame(S1, r1);
T=S1;
POChangeAngle(T);   // elimina punti che cambiano verso nella nuova posizione
POHiddenBySelf(S1);   // ordina la scansione e poi elimina punti nascosti da altri della stessa scansione (tutto nella nuova posizione)
T.readings.sort(readingsAngleCompare);
fuseScans(S1,T);    // fonde le due scansioni mettendo il risultato nella prima
			//  fusione: prende i punti presenti in tutte e due le scansioni (AND)

scanChangeRefFrame(S2, r2);
T=S2;
POChangeAngle(T);   // elimina punti che cambiano verso nella nuova posizione
POHiddenBySelf(S2);   // ordina la scansione e poi elimina punti nascosti da altri della stessa scansione (tutto nella nuova posizione)
T.readings.sort(readingsAngleCompare);
fuseScans(S2,T);    // fonde le due scansioni mettendo il risultato nella prima
			//  fusione: prende i punti presenti in tutte e due le scansioni (AND)
}
//---------------------------------------------------------------------------

void POChangeAngle(scan& S) // Purge Outliers that Change the Angle Direction
{   SIT i1 = S.readings.begin();
SIT i2 = S.readings.begin(); i2++;
bool erasing=false;
while (i1!=S.readings.end())
{   bool changeAngle=i1->angle-i2->angle>0.0 && i1->angle-i2->angle<M_PI;
	if (erasing && changeAngle)
	{   S.readings.erase(i1);
	}
	erasing=changeAngle;
	i1=i2;
	i2++;
}
}
//---------------------------------------------------------------------------

void POHiddenBySelf(scan& S)  // Purge Outliers that are Hidden by other points of the same scan
{   S.readings.sort(readingsAngleCompare);
SIT i1 = S.readings.begin();
SIT i2 = S.readings.begin(); i2++;
while (i2!=S.readings.end())
{   if (i2->angle-i1->angle<SOGLIA_ANGOLO_OUTLIERS)
	{   if (i1->distance-i2->distance>=SOGLIA_DISTANZA_OUTLIERS)     // elimina il punto piu lontano
		S.readings.erase(i1);
	else if (i2->distance-i1->distance>=SOGLIA_DISTANZA_OUTLIERS)
	{   S.readings.erase(i2);
		i2=i1;
	}
	}
	i1=i2;
	i2++;
}
}
//---------------------------------------------------------------------------

void fuseScans(scan& T1, const scan& T2)
{   SIT i1 = T1.readings.begin();
CSIT i2 = T2.readings.begin();
while (i1!=T1.readings.end())
{   while (i2->angle<i1->angle && i2!=T2.readings.end()) i2++;
	if (i2->angle!=i1->angle) i1=T1.readings.erase(i1);
	else i1++;
}
}
//---------------------------------------------------------------------------

void POHiddenByOther(scan& T1, scan& T2)
{   // 3 - elimina punti nascosti dall'altra scansione
SIT i1 = T1.readings.begin();
SIT i2 = T2.readings.begin();
bool eliminato;
while (i1!=T1.readings.end() && i2!=T2.readings.end())
{   eliminato=false;
	if (fabs(i2->angle-i1->angle)<SOGLIA_ANGOLO_OUTLIERS)
	{   if (i1->distance-i2->distance>SOGLIA_DISTANZA_OUTLIERS)
	{   i1=T1.readings.erase(i1);
		eliminato=true;
	} else if (i2->distance-i1->distance>SOGLIA_DISTANZA_OUTLIERS)
	{   i2=T2.readings.erase(i2);
		eliminato=true;
	}
	}
	if (!eliminato) {
	if (i2->angle<i1->angle) i2++; else i1++;
	}
}
}
//---------------------------------------------------------------------------

void POTooFar(scan& S) // Purge Outliers that are Too Far
{   SIT i = S.readings.begin();
while (i!=S.readings.end())
{   if (i->distance>SENSOR_RANGE_DEF) i=S.readings.erase(i);
	else i++;
}
}


typedef scan record_t;
typedef vector <record_t> data_t;

//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
// Remember, a record is a list of doubles separated by commas ','.
istream& operator >> ( istream& ins, record_t& record )
{
// make sure that the returned record contains only the stuff we read now
record.clear();

// read the entire line into a string (a CSV record is terminated by a newline)
string line;
getline( ins, line );

// now we'll use a stringstream to separate the fields out of the line
stringstream ss( line );
string field;
unsigned i=0;
double angle=0.0;
while (getline( ss, field, ',' ))
{ // for each field we wish to convert it to a double
	// (since we require that the CSV contains nothing but floating-point values)
	stringstream fs( field );
	double f = 0.0;  // (default value is 0.0)
	fs >> f;

	// ignore the first three values (we dont need them)
	// (in case of rawseeds datasets with sick laser: timestamp, numreadings, unknown)
	if (i>=3)
	{	// add the newly-converted field to the end of the record
	reading r;
	r.distance=f;
	r.angle=angle;
	angle+=SICK_STEP_ANGLE;
	record.readings.push_back( r );
	}
	i++;
}

// Now we have read a single line, converted into a list of fields, converted the fields
// from strings to doubles, and stored the results in the argument record, so
// we just return the argument stream as required for this kind of input overload function.
return ins;
}

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.
istream& operator >> ( istream& ins, data_t& data )
{
// make sure that the returned data only contains the CSV data we read here
data.clear();

// For every record we can read from the file, append it to our resulting data
record_t record;
while (ins >> record)
	{
	data.push_back( record );
	}

// Again, return the argument stream as required for this kind of input stream overload.
return ins;
}

//-----------------------------------------------------------------------------

int convert_from_rawseeds_CSV(const char *filename, const char *destination_dir)
// this function reads RAWSEEDS.ORG dataset (CSV file of laser scanner data)
// and converts it in a folder with individual scans (as used currently in sm_cli)
// TODO: I should do the opposite also ... and support different scan formats, scan pairs, sessions etc...
// the function returns the number of converted scans
{

data_t data;

ifstream infile( filename );
infile >> data;

// Complain if something went wrong.
if (!infile.eof())
{	cout << "Fooey!\n";
	return 1;
}

infile.close();

//for (unsigned n = 0; n < data.size(); n++)
for (unsigned n = 0; n < 500; n++)
{	string s;
		stringstream st; st << destination_dir << "scan" << n; // convert to string
		s=st.str(); cout << s << endl;
	data[ n ].write(s.c_str());
}

return data.size();
}


