
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "Random.h"

int main(int argc, char* argv[])	{

	int N = atoi(argv[1]);
	double alpha = atof(argv[2]);

	double mean = 0;
	double var = 0;

	int nzero = 0;
	for (int i=0; i<N; i++)	{
		double tmp1 = rnd::GetRandom().sGamma(alpha+1);
		double tmp2 = rnd::GetRandom().Uniform();
		double tmp3 = exp(log(tmp2) / alpha);
		double tmp = tmp1 * tmp3;
		// double tmp = rnd::GetRandom().sGamma(alpha);
		if (! tmp3)	{
			cerr << log(tmp2) << '\t' << log(tmp2) / alpha << '\n';
			nzero++;
		}
		mean += tmp;
		var += tmp * tmp;
	}
	mean /= N;
	var /= N;
	var -= mean*mean;

	cerr << "mean : " << mean << '\n';
	cerr << "var  : " << var << '\n';
	cerr << "nulls: " << nzero << '\n';

}

