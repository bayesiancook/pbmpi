
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixPhyloProcess.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Matrix PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void MatrixPhyloProcess::Unfold()	{
	if (condflag)	{
		cerr << "error in PhyloProcess::Unfold\n";
		exit(1);
	}
	/*
	static bool first = true;
	if (first) {
		first = false;
	}
	else {
		DeleteSuffStat();
	}
	*/
	DeleteSuffStat();
	DeleteMappings();
	ActivateSumOverRateAllocations();

	CreateMatrices();

	// UpdateSubstitutionProcess();

	CreateCondSiteLogL();
	CreateConditionalLikelihoods();
	UpdateConditionalLikelihoods();
}

void MatrixPhyloProcess::UpdateConditionalLikelihoods()	{

	PostOrderPruning(GetRoot(),condlmap[0]);

	// not necessary
	MultiplyByStationaries(condlmap[0]);
	ComputeLikelihood(condlmap[0]);

	PreOrderPruning(GetRoot(),condlmap[0]);

	// CheckLikelihood();
}


void MatrixPhyloProcess::Collapse()	{

	if (! condflag)	{
		cerr << "error in PhyloProcess::Collapse\n";
		exit(1);
	}
	// UpdateConditionalLikelihoods();
	DrawAllocations();
	SampleNodeStates();
	if (! dataclamped)	{
		if (rateprior)	{
			DrawAllocationsFromPrior();
		}
		if (profileprior)	{
			DrawProfileFromPrior();
		}
		SimulateForward();
	}
	DeleteCondSiteLogL();
	DeleteConditionalLikelihoods();
	InactivateSumOverRateAllocations(ratealloc);
	FillMissingMap();
	SampleSubstitutionMappings(GetRoot());
	DeleteMatrices();
	CreateSuffStat();
}

/*
void MatrixPhyloProcess::UpdateSubstitutionProcess()	{
	ActivateSumOverRateAllocations();
	CreateMatrices();
	UpdateMatrices();
	// DiagonaliseMatrices();
}
*/

double GTRPhyloProcess::LengthRelRateMove(double tuning, int nrep)	{

	double naccept = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogratio = - LogRRPrior() - LogLengthPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		int nbranch = MoveAllBranches(e);
		for (int i=0; i<GetNrr(); i++)	{
			rr[i] /= e;
		}
		deltalogratio += LogRRPrior() + LogLengthPrior();
		deltalogratio += (nbranch-GetNrr()) * m;

		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

		if (accepted)	{
			naccept++;
		}
		else	{
			MoveAllBranches(1.0/e);
			for (int i=0; i<GetNrr(); i++)	{
				rr[i] *= e;
			}
		}	
	}
	return naccept / nrep;
}

