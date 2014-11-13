
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatGTRProfileProcess.h"
	
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GeneralPathSuffStatGTRProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


double  GeneralPathSuffStatGTRProfileProcess::MoveRR(double tuning)	{
	GlobalUpdateRRSuffStat();
	int naccepted = 0;
	for (int i=0; i<GetNrr(); i++)	{
		double deltalogratio = - LogRRPrior() - ProfileSuffStatLogProb();
		double bk = rr[i];
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		rr[i] *= e;
		UpdateMatrices();
		deltalogratio += LogRRPrior() + ProfileSuffStatLogProb();
		deltalogratio += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);
		if (accepted)	{
			naccepted++;
		}
		else	{
			rr[i] = bk;
			UpdateMatrices();
		}
	}
	return ((double) naccepted) / GetNrr();
}

