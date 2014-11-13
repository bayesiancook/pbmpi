
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixOneProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MatrixOneProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void MatrixOneProfileProcess::Create(int innsite, int indim)	{
	if (! matrix)	{
		OneProfileProcess::Create(innsite,indim);
		// ??? CreateMatrix(); 
	}
}

void MatrixOneProfileProcess::Delete() {
	if (matrix)	{
		delete matrix;
		matrix = 0;
		OneProfileProcess::Delete();
	}
}

double MatrixOneProfileProcess::MoveProfile(double tuning, int n, int nrep)	{

	int naccepted = 0;
	double* bk = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bk[k] = profile[k];
	} 
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogProfilePrior() - ProfileSuffStatLogProb();
		double loghastings = ProfileProposeMove(profile,tuning,n,0);
		// double loghastings = ProfileProposeMove(profile[cat],tuning,n,0,g);
		UpdateProfile();
		deltalogprob += LogProfilePrior() + ProfileSuffStatLogProb();
		deltalogprob += loghastings;
		// int accepted = (g.RandU01() < exp(deltalogprob));
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted ++;
			for (int k=0; k<GetDim(); k++)	{
				bk[k] = profile[k];
			} 
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[k] = bk[k];
			} 
			UpdateProfile();
		}
	}
	delete[] bk;
	return naccepted / nrep;
}

