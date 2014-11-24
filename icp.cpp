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
 * icp.cpp
 */


#ifndef _MATH_H
    #include <math.h>
#endif


#ifndef ICPH
    #include "icp.h"
#endif

#ifndef utilsH
    #include "utils.h"
#endif

#ifdef DRAW_PNG
    #include "mydraw.h"  // for drawing
    const int FILENAME_SIZE = 255;
#endif

#include "mydraw.h"
namespace ICP
{

//---------------------------------------------------------------------------

static std::list<corrPoints> correspondances;
std::list<corrPoints> myCorrespondances;

static char corrSearchMethod; // holds the currently selected correspondance search method (FindCorr...)
void setCorrSearchMethod(char i)
{
    corrSearchMethod=i;
}

double corr_angle_treshold; //     = 0.148352986;
double sensor_range; //            = 8.0;
double stepangle; //               = 0.015707963;
double maxfloat; //                = 0.2;

double corr_distance_treshold; //  = 0.3
unsigned max_iteration; //           = 50
double min_ddisp; //               = 0.002
double min_rdisp; //               = 0.001


void setParameters(
        double ucorr_angle_treshold,
        double usensor_range,
        double ustepangle,
        double umaxfloat,
        double ucorr_distance_treshold,
        unsigned umax_iteration,
        double umin_ddisp,
        double umin_rdisp)
{
    corr_angle_treshold=ucorr_angle_treshold;
    sensor_range= usensor_range;
    stepangle=ustepangle;
    maxfloat=umaxfloat;
    corr_distance_treshold=ucorr_distance_treshold;
    max_iteration=umax_iteration;
    min_ddisp=umin_ddisp;
    min_rdisp=umin_rdisp;
}

bool draw_iterations;
unsigned verbose_level;
std::string     map_filename;
double map_size_x;
double map_size_y;

void setDebugParameters(bool udraw_iterations, unsigned uverbose_level, double umap_size_x, double umap_size_y, std::string umap_filename)
{
	draw_iterations=udraw_iterations;
	verbose_level=uverbose_level;
	map_size_x=umap_size_x;
	map_size_y=umap_size_y;
	map_filename=umap_filename;
}



//---------------------------------------------------------------------------
//-------------------------- CORRESPONDENCES --------------------------------
//---------------------------------------------------------------------------



void FindCorrSmartInjective(scan& S1, scan& S2, double corr_distance_treshold)
{   SIT infReading, supReading;
    infReading=supReading=S2.readings.begin();
    const int MAX_CYCLES = 8;

    // inizializza ra con la distanza molto grande
    SIT i1=S1.readings.begin();
    while (i1!=S1.readings.end())
    {   i1->ra=sensor_range;
        i1->match=false;
        i1++;
    }
    i1=S2.readings.begin();
    while (i1!=S2.readings.end())
    {   i1->ra=sensor_range;
        i1->match=false;
        i1++;
    }

    int cycle=1;
    while (cycle<MAX_CYCLES)
    {
        i1=S1.readings.begin();
        bool empty=true;
        while (i1!=S1.readings.end())
        {   if (!i1->match)
            {
                readingsInAngleRange(S1,S2,i1->angle-corr_angle_treshold,i1->angle+corr_angle_treshold,infReading,supReading);
                if (supReading!=S2.readings.end() && infReading!=S2.readings.end())
                {
                        // scorri tutti gli elementi di S2 all-interno di range angolo
                    SIT i2;
                    SIT im=S2.readings.end();
                    double minDist=sensor_range;

                    i2=infReading;
                    while (true)
                    {   if (i2==S2.readings.end()) i2=S2.readings.begin();
                        //if (!(i2->match<cycle) || i2->match==0)       // equivalente a questa riga dopo
                        if (i2->match==0 || (i2->match==cycle))  // qui ho corretto la seconda condizione in == invece di >=. Infatti significa di considerare il candidato solo se non matchato oppure se matchato nello stesso ciclo
                        {   double dist=distPolar(i1->angle, i1->distance, i2->angle, i2->distance);
                            if (minDist>dist && dist<corr_distance_treshold)
                            {   minDist=dist;
                                im=i2;
                            }
                        }
                        if (i2==supReading) break;
                        i2++;
                    }
                    if (im!=S2.readings.end())                      // se ho trovato un candidato
                    {   if (!im->match)
                        {   im->match=cycle;                                    //  segna il candidato come matched in questo cycle
                            i1->match=cycle;
                            im->ra=minDist;
                            corrPoints cp;
                            cp.pref=i1;
                            cp.pnew=im;
                            correspondances.push_back(cp);
                            empty=false;
                        } else
                        {   if (im->match==cycle && im->ra>minDist)
                            { // cerca tra i correspondances dalla fine per cp.pnew==im find
                                CPIT ci=correspondances.end();
                                while (ci!=correspondances.begin() && ci->pnew!=im) ci--;
                                if (ci->pnew==im)
                                {   ci->pref->match=false;                         // liberiamo quello precedente
                                    ci->pref=i1;
                                    i1->match=cycle;
                                    im->match=cycle;
                                    im->ra=minDist;
                                }
                            }
                        }
                    }
                }
            }
            i1++;
        }
        cycle++;
        if (empty) break;
    }
    // a questo punto cycle contiene il numero di cicli necessario per fare tutte le corrispondenze oppure e ugule a MAX_CYCLES

    // per eliminare gli incroci:
        // fai sort dei corripondenti in base all-angolo di una delle scansioni (diciamo S1)
        // cerca la prima coppia di corrispondenti con i1->match=cycle (cioe cycle massimo raggiunto)
        // scorri in su al prossimo elemento di S1 che ha corrispondenti e trova il suo corrispondente in S2
        // se sono di cycle minore controllo se fanno incrocio altrimenti la coppia diventa quella corrente
        // controllo incrocio: elemeto di S2 ha angolo minore mentre S1 ha angolo maggiore e viceversa? se si abbiamo trovato l'incrocio e possiamo cancellarlo!!!
        // finito il giro (in base all'angolo iniziale?) cala su Cycle minore e ripeti fino a Cycle=2;

    if (cycle<3 || correspondances.size()<3) return;
    correspondances.sort(compareByS1Angle);

    for (int i=0; i<5; i++)
    {
        CPIT ci=correspondances.begin();
        CPIT cj=ci; cj++;

        while (cj!=correspondances.end())
        {   // unwrap della fase (per confrontare gli angoli facciamo in modo che non distino piu di PI)
            double a1=ci->pref->angle;
            double a2=ci->pnew->angle;
            double b1=cj->pref->angle;
            double b2=cj->pnew->angle;
            if (a1-b1>M_PI)
                a1-=2.0*M_PI;
            else if (a1-b1<-M_PI)
                a1+=2.0*M_PI;
            if (a2-b2>M_PI)
                a2-=2.0*M_PI;
            else if (a2-b2<-M_PI)
                a2+=2.0*M_PI;

            if ((a1-b1>0.0 && a2-b2<0.0) || (a1-b1<0.0 && a2-b2>0.0))
            {    if (ci->pref->match>cj->pref->match)
                {   ci=correspondances.erase(ci);
                    cj++;
                } else
                if (cj->pref->match>ci->pref->match)
                {    cj=correspondances.erase(cj);
                    if (cj!=correspondances.end()) cj++;
                    ci++;
                } else
                {   cj++;
                    ci++;
                }
            } else
            {   cj++;
                ci++;
            }
        }
     }
}
//---------------------------------------------------------------------------

void FindCorrSmartInjective2(scan& S1, scan& S2, double corr_distance_treshold)
{   SIT infReading, supReading;
    infReading=supReading=S2.readings.begin();
    const int MAX_CYCLES = 8;

    // inizializza ra con la distanza molto grande
    SIT i1=S1.readings.begin();
    while (i1!=S1.readings.end())
    {   i1->ra=sensor_range;
        i1->match=false;
        i1++;
    }
    i1=S2.readings.begin();
    while (i1!=S2.readings.end())
    {   i1->ra=sensor_range;
        i1->match=false;
        i1++;
    }

    int cycle=1;
    while (cycle<MAX_CYCLES)
    {
        i1=S1.readings.begin();
        bool empty=true;
        while (i1!=S1.readings.end())
        {   if (!i1->match)
            {
                readingsInAngleRange(S1,S2,i1->angle-corr_angle_treshold,i1->angle+corr_angle_treshold,infReading,supReading);
                if (supReading!=S2.readings.end() && infReading!=S2.readings.end())
                {   // scorri tutti gli elementi di S2 all-interno di range angolo
                    SIT i2;
                    SIT im=S2.readings.end();
                    double minDist=sensor_range;

                    i2=infReading;
                    while (true)
                    {   if (i2==S2.readings.end()) i2=S2.readings.begin();
                        //if (!(i2->match<cycle) || i2->match==0)       // equivalente a questa riga dopo
                        if (i2->match==0 || (i2->match==cycle))  // qui ho corretto la seconda condizione in == invece di >=. Infatti significa di considerare il candidato solo se non matchato oppure se matchato nello stesso ciclo
                        {   double dist=distPolar(i1->angle, i1->distance, i2->angle, i2->distance);
                            if (minDist>dist && dist<corr_distance_treshold)
                            {   minDist=dist;
                                im=i2;
                            }
                        }
                        if (i2==supReading) break;
                        i2++;
                    }
                    if (im!=S2.readings.end())                      // se ho trovato un candidato
                    {   if (!im->match)
                        {   im->match=cycle;                                    //  segna il candidato come matched in questo cycle
                            i1->match=cycle;
                            im->ra=minDist;
                            corrPoints cp;
                            cp.pref=i1;
                            cp.pnew=im;
                            correspondances.push_back(cp);
                            empty=false;
                        } else
                        {   if (im->match==cycle && im->ra>minDist)
                            { // cerca tra i correspondances dalla fine per cp.pnew==im find
                                CPIT ci=correspondances.end();
                                while (ci!=correspondances.begin() && ci->pnew!=im) ci--;
                                if (ci->pnew==im)
                                {   ci->pref->match=false;                         // liberiamo quello precedente
                                    ci->pref=i1;
                                    i1->match=cycle;
                                    im->match=cycle;
                                    im->ra=minDist;
                                }
                            }
                        }
                    }
                }
            }
            i1++;
        }
        cycle++;
        if (empty) break;
    }
}
//---------------------------------------------------------------------------

void FindCorrFull(scan& S1, scan& S2, double corr_distance_treshold)
{   corrPoints cp;
    SIT i1=S1.readings.begin();
    while (i1!=S1.readings.end())
    {   double minDist=maxfloat;
        cp.pnew=S2.readings.end();
        SIT i2=S2.readings.begin();
        while (i2!=S2.readings.end())
        {   double dist=distPolar(i1->angle, i1->distance, i2->angle, i2->distance);
            if (minDist>dist && dist<corr_distance_treshold)
            {   minDist=dist;
                cp.pnew=i2;
            }
            i2++;
        }
        if (cp.pnew!=S2.readings.end()) // ok trovata corr
        {   cp.pref=i1;
            correspondances.push_back(cp);
        }
        i1++;
    }
}
//---------------------------------------------------------------------------

void FindCorrFullInjective(scan& S1, scan& S2, double corr_distance_treshold)
{

}

//---------------------------------------------------------------------------

void FindCorrSmart(scan& S1, scan& S2, double corr_distance_treshold)
{   SIT infReading, supReading;
    infReading=supReading=S2.readings.begin();
    SIT i1=S1.readings.begin();
    while (i1!=S1.readings.end())
    {   readingsInAngleRange(S1,S2,i1->angle-corr_angle_treshold,i1->angle+corr_angle_treshold,infReading,supReading);
        if (supReading!=S2.readings.end() && infReading!=S2.readings.end())
        {   SIT i2=infReading;
            SIT im=S2.readings.end();
            double minDist=sensor_range;
            while (true)
            {   if (i2==S2.readings.end()) i2=S2.readings.begin();
                double dist=distPolar(i1->angle, i1->distance, i2->angle, i2->distance);
                if (minDist>dist && dist<corr_distance_treshold)
                {   minDist=dist;
                    im=i2;
                }
                if (i2==supReading) break;
                i2++;
            }
            if (im!=S2.readings.end())                      // se ho trovato un candidato
            {   corrPoints cp;
                cp.pref=i1;
                cp.pnew=im;
                correspondances.push_back(cp);
            }
        }
        i1++;
    }
}
//---------------------------------------------------------------------------

void FindCorrLookup(scan& S1, scan& S2, double corr_distance_treshold)
{
//    SIT i2=S2.readings.begin();
//    while (i2!=S2.readings.end())
//    {   // coordinate X,Y gobali del reading di S2
//        double X=S2.pos.x+i2->distance*cos(S2.pos.rot+i2->angle);
//        double Y=S2.pos.y+i2->distance*sin(S2.pos.rot+i2->angle);
//
//        // discretizzazione
//        unsigned int ix=floor((X-LOOKUP_OFFSET_X)/LOOKUP_STEP_X+0.5);          // discretizzo X+LOOKUP_OFFSET_X cioč X espresso in coord locali nel frame
//        unsigned int iy=floor((Y-LOOKUP_OFFSET_Y)/LOOKUP_STEP_Y+0.5);          // -
//
//        // lookup
//        if (closestReading[ix][iy]!=NULLMARKER)    // non esiste null per gli iterator. Probabilmente questo confronto gli fa perdere tempo. trovare una soluzione
//        {   corrPoints cp;
//            cp.pref=closestReading[ix][iy];
//            cp.pnew=i2;
//            correspondances.push_back(cp);
//        }
//        i2++;
//    }
}
//---------------------------------------------------------------------------

void FindCorrPolar(scan& S1, scan& S2, double corr_distance_treshold)
{   SIT i1=S1.readings.begin();
    SIT i2=S2.readings.begin();
    while (i1!=S1.readings.end())
    {   double inf=i1->angle-stepangle/2.0;
        while (i2->angle<inf && i2!=S2.readings.end()) i2++;
        if (i2==S2.readings.end()) break;                                       // finiamo la ricerca, S1 e S2 ordinate=> non possono esserci altre corrispondenze, abbiamo esaurito S2
        if (i2->angle<inf+stepangle)                                            // ok! ci siamo, trovato corrispondenza
        {   double dist=distPolar(i1->angle, i1->distance, i2->angle, i2->distance);
            if (dist<corr_distance_treshold)
            {   corrPoints cp;
                cp.pref=i1;
                cp.pnew=i2;
                correspondances.push_back(cp);
            }
        }
        i1++;
    }
}
//---------------------------------------------------------------------------
position ICPMatchingStep(void)
{
    position p={0.0,0.0,0.0};
    int numCorr=myCorrespondances.size();
    if (numCorr>2) {

        double SigmaXr=0.0;
        double SigmaYr=0.0;
        double SigmaXn=0.0;
        double SigmaYn=0.0;
        double SigmaXrXn=0.0;
        double SigmaXrYn=0.0;
        double SigmaYrXn=0.0;
        double SigmaYrYn=0.0;


        // calcolo sommatorie della formula 13 dell-articolo di martinez
        CPIT C=myCorrespondances.begin();
        while (C!=myCorrespondances.end())
        {
            double xn=C->pnew->distance*cos(C->pnew->angle);        // xn,yn riferiti al sistema di riferimento di snew (dunque locale)
            double yn=C->pnew->distance*sin(C->pnew->angle);

            double xr=C->pref->distance*cos(C->pref->angle);        // xr,yr riferiti al sistema di riferimento di sref (dunque locale)
            double yr=C->pref->distance*sin(C->pref->angle);

            SigmaXr+=xr;
            SigmaYr+=yr;
            SigmaXn+=xn;
            SigmaYn+=yn;
            SigmaXrXn+=xr*xn;
            SigmaXrYn+=xr*yn;
            SigmaYrXn+=yr*xn;
            SigmaYrYn+=yr*yn;;

            C++;
        }

        double den=numCorr*SigmaXrXn+numCorr*SigmaYrYn-SigmaXr*SigmaXn-SigmaYr*SigmaYn;
        if (den!=0.0) {
            p.rot=atan((SigmaXr*SigmaYn+numCorr*SigmaYrXn-numCorr*SigmaXrYn-SigmaXn*SigmaYr)/den);
        }
        p.rot=angleUnwrap(p.rot);
        p.x=(SigmaXr-cos(p.rot)*SigmaXn+sin(p.rot)*SigmaYn)/numCorr;
        p.y=(SigmaYr-sin(p.rot)*SigmaXn-cos(p.rot)*SigmaYn)/numCorr;
    }
    return p;
}
//---------------------------------------------------------------------------

std::list<corrPoints> &  findCorrPoints(scan& S1, scan& S2, double corr_distance_treshold)
{
    correspondances.clear();
    if (S1.readings.size()>1 && S2.readings.size()>1)
    {   position p=S1.pos;

        scanChangeRefFrame(S1, S2.pos);                                          // put both scans in the same reference frame (i.e. put S1 in ref frame of S2)
        S1.readings.sort(readingsAngleCompare);

        switch (corrSearchMethod)
        {   case 0 : FindCorrFull(S1,S2,corr_distance_treshold); break;  // forma tutte le possibili coppie. Punti multipli della scansione S2 nei corrispondenti (2pti di S1 possono avere stesso punto di S2). Usato da Martinez et al.
            case 1 : FindCorrFullInjective(S1,S2,corr_distance_treshold); break;    // forma tutte le possibili coppie. No punti multipli (associazione esclusiva 1-1 tra punti). Ordina per distanze crescenti.
            case 2 : FindCorrSmart(S1,S2,corr_distance_treshold); break;      // versione subottimale ottimizzata di associazione esclusiva 1-1
            case 3 : FindCorrSmartInjective(S1,S2,corr_distance_treshold); break;     // versione semplificata che non elimina incroci (per fitness function)
            case 4 : FindCorrSmartInjective2(S1,S2,corr_distance_treshold); break;     // versione semplificata che non elimina incroci (per fitness function)
            case 5 : FindCorrLookup(S1,S2,corr_distance_treshold); break;  // versione adatta per scan matching globale. Precalcola aree discretizzate nella tabella per un lookup veloce durante la ricerca.
            case 6 : FindCorrPolar(S1,S2,corr_distance_treshold); break;   // versione supersemplice usata nel genetico da Martinez et al. L-unica diff. e' implementativa in quanto qua scorro le liste mentre loro ottnegono direttamente l'indice per l'array. Questa versione associa letture con stesso angolo discretizzato (e' come avere anglerange ridotto ad un unico stepangle).
        };

        scanChangeRefFrame(S1, p);                                      // torna
        S1.readings.sort(readingsAngleCompare);
    }
    return correspondances;
}

position icp(scan& sref, scan& snew, unsigned& numiter)
//  input parameters:
//    reference and new scan (the provided position of the new scan is used as initial position estimate)
//     corresponding point pairs maximum distance
//     convergence criteria parameters (tell the algorithm when to exit from iterations)
{

	//std::cout << "\tENTER ICP:" << myCorrespondances.size() << std::endl;
	position dp={0.0,0.0,0.0};                                                     // la posizione risultato di un iteration step
    position previous={0.0,0.0,0.0};                                            // il valore precedente per confronto
    int convergence=0;
    unsigned i=0;
    for(i=0; i<max_iteration; i++)
    {

        if (sref.readings.size()>1 && snew.readings.size()>1)
        {
            myCorrespondances=findCorrPoints(sref,snew,corr_distance_treshold);
//std::cout << "\tMY CORR:" << myCorrespondances.size() << std::endl;
        }
#ifdef DRAW_PNG
        if (draw_iterations)
        {
			// create image
        	// todo: put filename in icp.ini and pass it as option
			char filename[FILENAME_SIZE];
			snprintf(filename,FILENAME_SIZE,"./images/ite_%02d.png",i);
			drawing img(filename,map_size_x,map_size_y,-map_size_x/2.0,-map_size_y/2.0,map_filename.c_str());
			img.grid();

			// draw correspondances
			img.setpointsize(1);
			if (myCorrespondances.size()>1) {
			CPIT C=myCorrespondances.begin();
				while (C!=myCorrespondances.end())
				{   double X1=sref.pos.x+C->pref->distance*cos(sref.pos.rot+C->pref->angle);
					double Y1=sref.pos.y+C->pref->distance*sin(sref.pos.rot+C->pref->angle);
					double X2=snew.pos.x+C->pnew->distance*cos(snew.pos.rot+C->pnew->angle);
					double Y2=snew.pos.y+C->pnew->distance*sin(snew.pos.rot+C->pnew->angle);
					img.setcolor(0,0,0);
					img.line(X1,Y1,X2,Y2);
					img.setcolor(0,0,30000);
					img.readings(sref,sref.pos);
					img.setcolor(0,50000,0);
					img.readings(snew,snew.pos);
					C++;
				}
			}
			// close image
			img.write();
        }
#endif


		if (myCorrespondances.size()<3)
		{	position p;
			p.x=NAN;
			p.y=NAN;
			p.rot=NAN;
			numiter=i;
			return p;
        }


        previous=dp;            // memorizza per confronto
        dp=ICPMatchingStep();    // calcola posizione (1 iteration)

        // il risultato ottenuto da ICPMatchingStep() e' la posizione di snew
        // corretta riferita al sistema di coordinate locali di sref. Percio per
        // aggiornare snew devo prima trasformare in coord. globali (tranne la
        //rotazione che non serve)
        double d=sqrt(dp.x*dp.x+dp.y*dp.y);
        double a=pointAngle(dp.x,dp.y);
        double dx=cos(a+sref.pos.rot)*d;
        double dy=sin(a+sref.pos.rot)*d;
        snew.pos.x=sref.pos.x+dx;
        snew.pos.y=sref.pos.y+dy;
        snew.pos.rot=angleUnwrap(sref.pos.rot+dp.rot);
        // usciamo dal loop se la convergenza si e' rallentata troppo
        // (distance displacement e rotation displacement troppo piccoli
        // tra due step)
        double ddisp=sqrt((dp.x-previous.x)*(dp.x-previous.x)+(dp.y-previous.y)*(dp.y-previous.y));
        double tr1=angleUnwrap(dp.rot); if (tr1 > M_PI) tr1=tr1-2*M_PI;
        double tr2=angleUnwrap(previous.rot); if (tr2 > M_PI) tr2=tr2-2*M_PI;
        double rdisp=sqrt((tr1-tr2)*(tr1-tr2)); if (rdisp > M_PI) rdisp=rdisp-2*M_PI;
        //std::cout << "\tMY CORR: dist:" << ddisp << "    rot:" << rdisp << "   rot:" << fabs(rdisp-2*M_PI) << "  conv:" << convergence << std::endl;
        if (ddisp<min_ddisp && (rdisp<min_rdisp || fabs(rdisp-2*M_PI)<min_rdisp)) convergence++; else convergence=0;
        if (i>2 && convergence>2)
        {	convergence=0;
        	break;
        }
        // se 3 volte di fila lo spostamento e' minore della soglia e abbiamo
        // fatto almeno 2 iterazioni iniziali, abbiamo finito
    }
    numiter=i;
    std::cout << "\tICP numiter= "<< numiter << std::endl;
    return snew.pos;
}


} // namespace
