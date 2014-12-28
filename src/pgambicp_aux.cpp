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


#include "pgambicp_aux.h"

#ifdef DRAW_PNG
		#include "mydraw.h"
#endif

namespace PGAMBICP
{
	scan sr,sn;

	void init(const char* filename)
	{	POLARGA::init(filename);
		MBICP::init(filename);
	}

	void pre_match(const scan &sref, const scan &snew)
	{	POLARGA::pre_match(sref,snew);
		sr.copy(sref); // we need to store these two for the mbicp step
		sn.copy(snew);
	}


	matching_result match()
	{
		matching_result mres;

		// step 1 of hybrid algorithm: PolarGA
		mres=POLARGA::match();
		//std::cout<<"step1 estimate: "<<mres.p.x<<", "<<mres.p.y<<", "<<mres.p.rot <<std::endl;

#ifdef DRAW_PNG
//		drawing img1("images/pgambicp_step1.png", 25, 25,-25/2.0,-25/2.0,"./cfg/bitmaps/empty_1000x1000.png");
//		img1.grid();
//		img1.setpointsize(6);
//		img1.setcolor(0,0,30000);
//		img1.readings(sr,sr.pos);
//		img1.setpointsize(4);
//		img1.setcolor(0,0,50000);
//		img1.readings(sn,mres.p);
//		img1.write();
#endif
		readingsChangeRefFrame(sn.readings,mres.p,sn.pos);
		position pga=mres.p; // we need to store this because step1 estimate will be combined with step 2 refined estimate
		MBICP::pre_match(sr,sn);
		mres=MBICP::match();

#ifdef DRAW_PNG
//		drawing img3("images/pgambicp_step2.png", 25, 25,-25/2.0,-25/2.0,"./cfg/bitmaps/empty_1000x1000.png");
//		img3.grid();
//		img3.setpointsize(6);
//		img3.setcolor(0,0,30000);
//		img3.readings(sr,sr.pos);
//		img3.setpointsize(4);
//		img3.setcolor(0,0,50000);
//		img3.readings(sn,mres.p);
//		img3.write();
#endif
		mres.p.x=mres.p.x+pga.x;
		mres.p.y=mres.p.y+pga.y;
		mres.p.rot=angleUnwrap(mres.p.rot+pga.rot);
		return mres;
	}

	void post_match()
	{	MBICP::post_match();
		POLARGA::post_match();
	}

	void deinit()
	{	MBICP::deinit();
		POLARGA::deinit();
	}

} // namespace PGAMBICP
