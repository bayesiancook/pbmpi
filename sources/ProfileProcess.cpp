
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void ProfileProcess::Create(int innsite, int indim)	{
	if (nsite || dim)	{
		if (nsite != innsite)	{
			cerr << "error in phyloprocess creation: non matching number of sites\n";
			exit(1);
		}
		if (dim != indim)	{
			cerr << "error in phyloprocess creation: non matching number of states\n";
			cerr << dim << '\t' << indim << '\n';
			exit(1);
		}
	}
	else	{
		nsite = innsite;
		dim = indim;
	}
}

