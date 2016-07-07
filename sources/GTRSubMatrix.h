
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GTRSUBMATRIX_H
#define GTRSUBMATRIX_H

#include "BiologicalSequences.h"
#include "SubMatrix.h"

class GTRSubMatrix : public virtual SubMatrix	{

	public:


				GTRSubMatrix(int nstate, const double* rr, const double* stat, bool innormalise = false);
				~GTRSubMatrix() {};

	int			GetNRelativeRate() {return Nrr;}
	double 			RelativeRate(int i, int j) {return mRelativeRate[rrindex(i,j,GetNstate())];}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	static void SetMinStat(double in)	{
		minstat = in;
	}

	// make a copy of the entries (not of the pointer)
	// void 			CopyStationary(const double* instat);

	// copy the pointer
	void			SetRelativeRate(const double* inrelrate) {mRelativeRate = inrelrate;}

	protected:

	void 			ComputeArray(int state);
	void 			ComputeStationary();

	// data members
	const double* mRelativeRate;
	const double* ExternalStat;
	int Nrr;

	static int minstat;

};



class LGSubMatrix : public GTRSubMatrix	{

	public:

				LGSubMatrix(const double* stat, bool innormalise = false) : SubMatrix(Naa,innormalise), GTRSubMatrix(Naa,0,stat,innormalise)	{
					mRelativeRate = LG_RR;
				}
};

#endif 
