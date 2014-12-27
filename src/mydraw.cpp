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
 * mydraw.cpp
 *
 *  Created on: Jul 7, 2013
 *      Author: klenac
 */




#include <time.h>
#include "mydraw.h"

drawing::drawing(  const char * filename,      // full path
                        const double w,
                        const double h)
{   envwidth=w; envheight=h;
    X0 = 0.0;
    Y0 = 0.0;
    cols = 512;
    rows = 512;
    pngwriter imageTTT(rows, cols, 1.0, filename);
    img=imageTTT;
    fx=cols/envwidth;
    fy=rows/envheight;
    pointsize=3;
    fontsize=10;
    red=green=blue=10000;
    arrowsize=5;      // size of head of arrow in pixels
    arrowangle=3.141;     // in radiant
    arrowlength=1.0;
}
// -----------------------------------------------------------------------------

drawing::drawing(   const char * filename,      // full path
                    const double w,
                    const double h,
                    const double x0,
                    const double y0,
                    const char * mapname)
// cols,rows will be deduced from the mapname image loaded from disk
{   envwidth=w; envheight=h; X0=x0; Y0=y0;
    pngwriter imageTTT(1000, 1000, 1.0, filename);
    img=imageTTT;
    img.readfromfile(mapname);
    rows=img.getheight();
    cols=img.getwidth();
    fx=cols/envwidth;
    fy=rows/envheight;
    pointsize=3;
    fontsize=10;
    red=green=blue=10000;
    arrowsize=5;      // size of head of arrow in pixels
    arrowangle=3.141;     // in radiant
    arrowlength=1.0;


    /*
    printf("\n init_png: filename  : %s",filename);
    printf("\n init_png: envwidth     : %f",envwidth);
    printf("\n init_png: envheight     : %f",envheight);
    printf("\n init_png: mapname   : %s",mapname);
    printf("\n init_png: rows: %d",rows);
    printf("\n init_png: columns : %d",cols);
    printf("\n init_png: fx: %f",fx);
    printf("\n init_png: fy : %f\n",fy);
    */
}
// -----------------------------------------------------------------------------

drawing::drawing(   const char * filename,      // full path
                    const double w,
                    const double h,
                    const double x0,
                    const double y0,
                    const unsigned c,
                    const unsigned r)
{   cols=c; rows=r; envwidth=w; envheight=h; X0=x0; Y0=y0;
    pngwriter imageTTT(rows, cols, 1.0, filename);
    fx=cols/envwidth;
    fy=rows/envheight;
    pointsize=3;
    fontsize=10;
    red=green=blue=10000;
    arrowsize=5;      // size of head of arrow in pixels
    arrowangle=3.141;     // in radiant
    arrowlength=1.0;
}




void drawing::readings(const scan & s, const position p)
{
    CSIT i=s.readings.begin();
    while ( i!=s.readings.end() )
    {
        // coordinates of reading (from position p)
        double x=p.x+i->distance*cos(p.rot+i->angle);
        double y=p.y+i->distance*sin(p.rot+i->angle);
        point(x,y);
        i++;
    };
}
// -----------------------------------------------------------------------------

void drawing::grid()
{   double X=ceil(X0);
    double Y=ceil(Y0);
    for (double a=X;a<X0+envwidth; a+=1.0)
        for (double b=Y;b<Y0+envheight;b+=1.0)
        {    int x=fx*(a-X0);
            int y=fy*(b-Y0);
            img.line(x,y,x+cols,y,10, 10, 10);
            img.line(x,y,x,y+rows,10, 10, 10);
        }
}
// ----------------------------------------------------------------------------

void drawing::arrow(const position p)
{   double H=sin(p.rot)*arrowlength;
    double W=cos(p.rot)*arrowlength;
    int x1=fx*(p.x-X0)-W/2*fx;
    int y1=fy*(p.y-Y0)-H/2*fy;
    int x2=fx*(p.x-X0)+W/2*fx;
    int y2=fy*(p.y-Y0)+H/2*fy;
    img.arrow(x1,y1,x2,y2,arrowsize,arrowangle,red,green,blue);
} // draw()
// -----------------------------------------------------------------------------
void drawing::init(  const char * filename,      // full path
                        const double w,
                        const double h)
{   envwidth=w; envheight=h;
    X0 = 0.0;
    Y0 = 0.0;
    cols = 512;
    rows = 512;
    pngwriter imageTTT(rows, cols, 1.0, filename);
    img=imageTTT;
    fx=cols/envwidth;
    fy=rows/envheight;
    pointsize=3;
    fontsize=10;
    red=green=blue=10000;
    arrowsize=5;      // size of head of arrow in pixels
    arrowangle=3.141;     // in radiant
    arrowlength=1.0;
}
// -----------------------------------------------------------------------------

void drawing::init(   const char * filename,      // full path
                    const double w,
                    const double h,
                    const double x0,
                    const double y0,
                    const char * mapname)
// cols,rows will be deduced from the mapname image loaded from disk
{   envwidth=w; envheight=h; X0=x0; Y0=y0;
    pngwriter imageTTT(1000, 1000, 1.0, filename);
    img=imageTTT;
    img.readfromfile(mapname);
    rows=img.getheight();
    cols=img.getwidth();
    fx=cols/envwidth;
    fy=rows/envheight;
    pointsize=3;
    fontsize=10;
    red=green=blue=10000;
    arrowsize=5;      // size of head of arrow in pixels
    arrowangle=3.141;     // in radiant
    arrowlength=1.0;


    /*
    printf("\n init_png: filename  : %s",filename);
    printf("\n init_png: envwidth     : %f",envwidth);
    printf("\n init_png: envheight     : %f",envheight);
    printf("\n init_png: mapname   : %s",mapname);
    printf("\n init_png: rows: %d",rows);
    printf("\n init_png: columns : %d",cols);
    printf("\n init_png: fx: %f",fx);
    printf("\n init_png: fy : %f\n",fy);
    */
}
// -----------------------------------------------------------------------------

void drawing::init(   const char * filename,      // full path
                    const double w,
                    const double h,
                    const double x0,
                    const double y0,
                    const unsigned c,
                    const unsigned r)
{   cols=c; rows=r; envwidth=w; envheight=h; X0=x0; Y0=y0;
    pngwriter imageTTT(rows, cols, 1.0, filename);
    fx=cols/envwidth;
    fy=rows/envheight;
    pointsize=3;
    fontsize=10;
    red=green=blue=10000;
    arrowsize=5;      // size of head of arrow in pixels
    arrowangle=3.141;     // in radiant
    arrowlength=1.0;
}

