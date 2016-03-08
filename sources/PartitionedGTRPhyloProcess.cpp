
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionedGTRPhyloProcess.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Matrix PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PartitionedGTRPhyloProcess::Unfold()	{
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

void PartitionedGTRPhyloProcess::Collapse()	{

	if (! condflag)	{
		cerr << "error in PhyloProcess::Collapse\n";
		exit(1);
	}
	// UpdateConditionalLikelihoods();
	DrawAllocations();
	SampleNodeStates();
	DeleteCondSiteLogL();
	DeleteConditionalLikelihoods();
	InactivateSumOverRateAllocations(ratealloc);
	SampleSubstitutionMappings(GetRoot());
	DeleteMatrices();
	CreateSuffStat();
}

double PartitionedGTRPhyloProcess::LengthRelRateMove(double tuning, int nrep)	{

	double naccept = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogratio = - LogRRPrior() - LogLengthPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		int nbranch = MoveAllBranches(e);
		for(int p = 0; p < GetNpart(); p++)
		{
			for (int i=0; i<GetNrr(); i++)	{
				rr[p][i] /= e;
			}
		}
		deltalogratio += LogRRPrior() + LogLengthPrior();
		deltalogratio += (nbranch-GetNpart()*GetNrr()) * m;

		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

		if (accepted)	{
			naccept++;
		}
		else	{
			MoveAllBranches(1.0/e);
			for(int p = 0; p < GetNpart(); p++)
			{
				for (int i=0; i<GetNrr(); i++)	{
					rr[p][i] *= e;
				}
			}
		}	
	}
	return naccept / nrep;
}

double PartitionedGTRGammaPhyloProcess::LengthMultiplierMove(double tuning, int nrep)	{

	int nmult = PartitionedDGamRateProcess::GetNpart();

	double naccept = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogratio = - LogLengthPrior() - LogMultiplierPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		int nbranch = MoveAllBranches(e);
		for(int p = 0; p < nmult; p++)
		{
			ratemult[p] /= e;
		}
		deltalogratio += LogLengthPrior() + LogMultiplierPrior();
		deltalogratio += (nbranch-nmult) * m;

		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

		if (accepted)	{
			naccept++;
		}
		else	{
			MoveAllBranches(1.0/e);
			for(int p = 0; p < nmult; p++)
			{
				ratemult[p] *= e;
			}
		}
	}
	return naccept / nrep;
}

double PartitionedGTRGammaPhyloProcess::MultiplierRelRateMove(double tuning, int nrep)	{


	int nmult = PartitionedDGamRateProcess::GetNpart();
	int nrrpart = PartitionedGTRProfileProcess::GetNpart();

	map<int, vector<int> > gtrmults;

	// find each GTR partition and associated multipliers
	for(int p = 0; p < nmult; p++)
	{
		int site = PartitionedDGamRateProcess::GetPartSites(p).front();

		int rrpart = PartitionedGTRProfileProcess::GetSitePart(site);
		if(PartitionedGTRProfileProcess::GetPartType(rrpart) == "None")
		{
			gtrmults[rrpart].push_back(p);
		}
	}

	// iterate over each gtr partition and perform the move separately on each of them
	double naccept = 0;
	for (map<int, vector<int> >::iterator it = gtrmults.begin(); it != gtrmults.end(); it++)
	{
		int rrpart = it->first;
		vector<int> mults = it->second;

		for (int rep=0; rep<nrep; rep++)	{
			double deltalogratio = - LogRRPrior() - LogMultiplierPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);

			//perform the move
			for (int i=0; i<GetNrr(); i++)	{
				rr[rrpart][i] *= e;
			}
			for(int i = 0; i < mults.size(); i++)
			{
				ratemult[mults[i]] /= e;
			}

			deltalogratio += LogRRPrior() + LogMultiplierPrior();
			deltalogratio += (GetNrr()-mults.size()) * m;

			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

			if (accepted)	{
				naccept++;
			}
			else	{
				//undo the move
				for (int i=0; i<GetNrr(); i++)	{
					rr[rrpart][i] /= e;
				}
				for(int i = 0; i < mults.size(); i++)
				{
					ratemult[mults[i]] *= e;
				}
			}
		}
	}
	return naccept / nrep;
}

