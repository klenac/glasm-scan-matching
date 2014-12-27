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
#ifndef drawingH
#define drawingH
//---------------------------------------------------------------------------


#include <pngwriter.h>
	/* needs packages: libpngwriter0c2, libpngwriter0-dev */

#ifndef utilsH
    #include "utils.h"
#endif



/*
 * You draw in the image which represents the environment of size envwidth,
 * envheight (meters). Pixel (0,0) of lower left corner is mapped in some point
 * (X01,Y0) in meters of the environment. Upper right pixel of the image therefore
 * represents a point (X0+envwidth,Y0+envheight) a generic pixel (a,b) represents
 * a point (X0+envwidth/cols*a,Y0+envheight/rows*b).
 */


class drawing
{
private:
    double envwidth;            // width of the depicted environment (meters)
    double envheight;           // height

    // lower left corner of the rectangular portion of the environment that we want to depict in the image
    double X0;                  // coordinate X in meters in the environment of the leftmost pixel in the image
    double Y0;                  // coordinate Y in meters in the environment of the lowest pixel in the image

    unsigned cols;              // number of columns (image width in pixels)
    unsigned rows;              // number of rows    (image height in pixels)

    unsigned pointsize;
    unsigned fontsize;

    pngwriter img;
    int red, green, blue;       // current color

    unsigned    arrowsize;      // size of head of arrow in pixels
    double      arrowangle;     // in radiant
    double      arrowlength;    // in meters

    double fx;                  // ratio between environment and image size in x
    double fy;                  //


public:
    // constructors
    drawing(){};
    drawing(const char * filename, const double envwidth, const double envheight);
    drawing(const char * filename,
            const double envwidth,
            const double envheight,
            const double X0,
            const double Y0,
            const char * mapname);
    drawing(const char * filename,
            const double envwidth,
            const double envheight,
            const double X0,
            const double Y0,
            const unsigned cols,
            const unsigned rows);
    void init(const char * filename, const double envwidth, const double envheight);
    void init(const char * filename,
            const double envwidth,
            const double envheight,
            const double X0,
            const double Y0,
            const char * mapname);
    void init(const char * filename,
            const double envwidth,
            const double envheight,
            const double X0,
            const double Y0,
            const unsigned cols,
            const unsigned rows);



    // setters
    void setcolor(const int r, const int g, const int b) { red=r; green=g; blue=b; }
    void setpointsize(const unsigned ps) { pointsize=ps;}
    void setfontsize(const unsigned fs) { fontsize=fs;}
    void setenvsize(const double w, const double h) { envwidth=w, envheight=h;}
    void setimgsize(const unsigned c, const unsigned r) { cols=c, rows=r;}
    void setcorner(const double x, const double y) { X0=x; Y0=y;}
    void setarrow(const unsigned s, const double a, const double l) { arrowsize=s, arrowangle=a; arrowlength=l;}

    // drawing procedures
    void point(const double x, const double y)
    {   img.filledcircle(fx*(x-X0), fy*(y-Y0), pointsize, red,green,blue);
    }
    void pixel(const double x, const double y)
    {   img.plot(fx*(x-X0), fy*(y-Y0), red,green,blue);
    }
    void line(const double x1, const double y1, const double x2, const double y2)
    {   img.line(fx*(x1-X0),fy*(y1-Y0),fx*(x2-X0),fy*(y2-Y0),red,green,blue);
    }
    void text(const double x, const double y, char * mytext)
    {    img.plot_text("/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans.ttf", fontsize, x, y, 0, mytext, red, green, blue);
    }
    void readings(const scan & s, const position p);
    void grid();
    void arrow(const position p);

    //
    void write() { img.write_png(); }



};

#endif
