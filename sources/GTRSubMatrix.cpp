
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GTRSubMatrix.h"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//		 GTRSubMatrix
// ---------------------------------------------------------------------------


GTRSubMatrix::GTRSubMatrix(int inNstate, const double* rr, const double* stat, bool innormalise) : SubMatrix(inNstate, innormalise)	{

	Nrr = Nstate * (Nstate-1) / 2;
	mRelativeRate = rr;
	ExternalStat = stat;
	CorruptMatrix();
}

void GTRSubMatrix::ComputeStationary()	{
	for (int k=0; k<Nstate; k++)	{
		mStationary[k] = ExternalStat[k];
		// mStationary[k] = instat[k];
	}
}

// ---------------------------------------------------------------------------
//		 ComputeArray
// ---------------------------------------------------------------------------

void	GTRSubMatrix::ComputeArray(int i)	{

	if (mRelativeRate)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				Q[i][j] = RelativeRate(i,j) * mStationary[j];
				total += Q[i][j];
			}
		}

		// should always ensure that the diagonal entry of the matrix Q[i][i] is such that 
		// the sum over all entries of the row is equal to 0
		Q[i][i] = - total;
	}
	else	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] = mStationary[j];
		}
		Q[i][i] -= 1;
	}
}

