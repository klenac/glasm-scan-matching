/*  
    Copyright (C) 2014  Kristijan Lenac
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
//---------------------------------------------------------------------------

#ifndef utilsH
#define utilsH

//---------------------------------------------------------------------------

#include <string>
#include <list>
#include <vector>

// converte i radianti in gradi
//
#define rad_gra(r) ((r)*360)/(2*M_PI)

// converte i gradi in radianti
//
#define gra_rad(g) ((g)*2*M_PI)/360

struct point
{
	double x,y;
};

struct position
{
	double x;
	double y;
	double rot;
};

struct relPosition
{
	double direction;
	double distance;
};


struct wall
{
	std::string surface;
	double x1, y1, x2, y2;
};


typedef std::list<double>* list_ptr_type;
/*typedef std::map<string, list_ptr_type> surfaces_map;*/
typedef std::list<wall> wall_list;


struct reading
{
	double angle;
	double distance;
	int    match;
	double ra;
	double rb;
	double incidence;
	double fiterr;
};


typedef std::list<reading>::const_iterator CSIT;
typedef std::list<reading>::iterator SIT;


struct corrPoints
{
	SIT         pref;
	SIT         pnew;
	double      det_Pij;
	double      P_ij[2][2];
	double      inv_P_ij[2][2];
};


typedef std::list<corrPoints>::iterator CPIT;

class scan
{
public:
	position pos;

	std::list<reading> readings;

	scan(double x=0, double y=0, double a=0)
	{
		pos.x=x;
		pos.y=y;
		pos.rot=a;
		readings.clear();
	};

	void clear()
	{
		pos.x=0.0;
		pos.y=0.0;
		pos.rot=0.0;
		readings.clear();
	};

	void copy(const scan &pSx)
	{   // Copy of scan
		CSIT si;
		readings.clear();
		si=pSx.readings.begin();
		while (si!=pSx.readings.end())
			readings.push_back(*si++);
		pos=pSx.pos;
	};

	void setPos(double x, double y, double a)
	{
		pos.x=x;
		pos.y=y;
		pos.rot=a;
	};

	void read(const char *filename);

	void read(std::string &rawseeds_scan);

	void write(const char *filename);
};

class scan_map
{
public:

	scan scan_node;

	scan * scan_next, scan_prev;

	std::list<scan> near_scan;

	bool scaned; // All points are drawed

	scan_map(scan &pScanNode, scan &pScanPrevNode){
		scan_node=pScanNode;
		scan_prev=pScanPrevNode;
		near_scan.push_back(pScanPrevNode);
	}

};


struct matching_result
{
	position p;
	unsigned int numiter;
	int status;
};

typedef std::list<scan_map>::iterator SMI;

//---------------------------------------------------------------------------

void assign_pos(position &posTo, const position &posFrom);
void calc_RanRot(position pPosOn, position pPosTo, double& newPosRot, double& newPosRan);
void calc_xy(double newPosRot, double newPosRan, position& pPos);

double maxRange(void);
void   setMaxRange(double pMRange);

void    locToGlob(const position, double, double, double&, double&);
void    globToLoc(const position, double, double, double&, double&);
double  pointAngle(double, double);
bool    CCW(point a, point b, point c);
bool    intersect(point a, point b, point c, point d);
double  angleInternalDifference(double, double);
double  angleUnwrap(double alfa);
bool    angleInArc(double, double, double);
double  gauss(double var, double med);
double randDouble(double low, double high);
double  distPolar(double a1, double d1, double a2, double d2);
int     _random(unsigned lim);
void    scanChangeRefFrame(scan&, position);
void    readingsChangeRefFrame(std::list<reading>&, position, position);
void    scanMedianFilter(scan&, const unsigned N = 3, const double DIST = 0.06); //0.06
void    scanResample(scan&, const double DIST = 0.03); //0.03
void 	scanDecimate(scan&, const unsigned ORDER = 2, const unsigned OFFSET = 0);
bool    readingsAngleCompare(reading r1, reading r2);
bool    compareByS1Angle(corrPoints c1, corrPoints c2);
void readingsInAngleRange( scan& S1, scan& S2, double infAngle, double supAngle
						 , SIT& infReading, SIT& supReading);

// procedures used in outlier detection and removal
void outliers(scan&, scan&);
void POChangeAngle(scan&);
void POHiddenBySelf(scan&);
void POHiddenByOther(scan&, scan&);
void POTooFar(scan&);
void fuseScans(scan&, const scan&);


int convert_from_rawseeds_CSV(const char *, const char *);


/*
template<typename T>
inline bool isnan(T value)
{
return value != value;

}
*/

#endif
