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


#ifndef GCP_AUXH
#define GCP_AUXH

#include "simulation.h"
#include "utils.h"
#include "gcp.h"


namespace GCP
{

	void read_parameters(const char* filename);
	void init();                                        // one time initialization
	void pre_match(const scan &snew, const scan &sref); // before matching
	matching_result match();                            // matching
	void post_match();                                  // after matching
	void deinit();                                      // one time deinitialization

}

#endif
