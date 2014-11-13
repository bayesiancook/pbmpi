
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/



#include "MatrixProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MatrixProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


double MatrixProfileProcess::ProfileProposeMove(double* profile, double tuning, int n, int K, double statmin)	{ // n==0dirichlet resampling, otherwise, vase communiquants
	if (! statmin)	{
		statmin = stateps;
	}
	if (! K)	{
		K = dim;
	}
	double ret = 0;
	if (!n)	{ // dirichlet
		cerr << "dirichlet move\n";
		exit(1);
		double* oldprofile = new double[K];
		for (int i=0; i<K; i++)	{
			oldprofile[i] = profile[i];
		}
		double total = 0;
		for (int i=0; i<K; i++)	{
			profile[i] = rnd::GetRandom().sGamma(tuning*oldprofile[i]);
			if (profile[i] == 0)	{
				cerr << "error in dirichlet resampling : 0 \n";
				exit(1);
			}
			total += profile[i];
		}

		double logHastings = 0;
		for (int i=0; i<K; i++)	{
			profile[i] /= total;

			logHastings += - rnd::GetRandom().logGamma(tuning*oldprofile[i]) + rnd::GetRandom().logGamma(tuning*profile[i])
						-  (tuning*profile[i] -1.0) * log(oldprofile[i]) + (tuning * oldprofile[i] -1.0) * log(profile[i]);
		}

		delete[] oldprofile;
		return logHastings;
	}
	else	{
		if (2*n > K)	{
			n = K / 2;
		}
		int* indices = new int[2*n];
		rnd::GetRandom().DrawFromUrn(indices,2*n,K);
		for (int i=0; i<n; i++)	{
			int i1 = indices[2*i];
			int i2 = indices[2*i+1];
			double tot = profile[i1] + profile[i2];
			double x = profile[i1];
			
			// double h = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double h = tot * tuning * (rnd::GetRandom().Uniform() - 0.5);
			/*
			int c = (int) (h / (2 * tot));
			h -= c*2*tot;
			*/
			x += h;
			while ((x<0) || (x>tot))	{
				if (x<0)	{
					x = -x;
				}
				if (x>tot)	{
					x = 2*tot - x;
				}
			}
			profile[i1] = x;
			profile[i2] = tot - x;
		}
		delete[] indices;
	}
	return ret;
}
