
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GTRMixtureProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GTRMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GTRMixtureProfileProcess::Create(int innsite, int indim)	{
	GTRProfileProcess::Create(innsite,indim);
	MatrixMixtureProfileProcess::Create(innsite,indim);
}

void GTRMixtureProfileProcess::Delete() {
	MatrixMixtureProfileProcess::Delete();
	GTRProfileProcess::Delete();
}

void GTRMixtureProfileProcess::CreateMatrix(int k)	{
	if (matrixarray[k])	{
		cerr << "error in gen path suff stat gtr dp profile process: matrixarray is not 0\n";
		cerr << matrixarray[k]->GetNstate() << '\n';
		exit(1);
	}
	// matrixarray[k] = new LGSubMatrix(profile[k],false);
	matrixarray[k] = new GTRSubMatrix(GetDim(),rr,profile[k],false);
}

void GTRMixtureProfileProcess::UpdateMatrix(int k)	{
	/*
	if (! matrixarray[k])	{
		cerr << "error in GeneralPathSuffStatGTRMixtureProfileProcess::UpdateComponent: null matrix\n";
		exit(1);
	}
	*/
	if (matrixarray[k])	{
	GetGTRMatrix(k)->CorruptMatrix();
	}
	// GetGTRMatrix(k)->UpdateMatrix();
}
