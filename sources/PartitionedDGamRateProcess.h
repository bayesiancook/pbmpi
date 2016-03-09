
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTDGAMRATE_H
#define PARTDGAMRATE_H

#include "RateProcess.h"
#include "Partition.h"

class PartitionedDGamRateProcess : public virtual RateProcess, public PartitionProcess {

	public:

	PartitionedDGamRateProcess() : Ncat(0), rate(0) {}
	virtual ~PartitionedDGamRateProcess() {}

	double GetAlpha(int part) {return alpha[part];}

	double GetRateMultiplier(int part) {return ratemult[part];}

	double GetAlpha() {
		double total = 0.0;
		for(int i = 0; i < GetNpart(); i++)
		{
			total += alpha[i];
		}

		return total / GetNpart();
	}

	double GetMultiplierEntropy();

	int GetNrate(int site)	{
		if (SumOverRateAllocations())	{
			return Ncat;
		}
		return 1;
	}

	int GetNcat() {return Ncat;}

	double GetRate(int site, int cat = 0)	{
		// cat should be == 0
		if (SumOverRateAllocations())	{
			return rate[GetSitePart(site)][cat]*ratemult[GetSitePart(site)];
		}
		return rate[GetSitePart(site)][alloc[site]]*ratemult[GetSitePart(site)];
	}

	double GetRateWeight(int site, int cat)	{
		if (SumOverRateAllocations())	{
			return 1.0/Ncat;
		}
		return 1.0;
	}

	double GetMultHyper()	{
		return multHyper;
	}

	void ActivateSumOverRateAllocations() {
		sumflag = true;
	}

	void InactivateSumOverRateAllocations(int* ratealloc) {
		for (int i=0; i<GetNsite(); i++)	{
			alloc[i] = ratealloc[i];
		}
		sumflag = false;
	}

	double GetPriorMeanRate()	{
		double total = 0;
		for (int p=0; p<GetNpart(); p++)	{
			double part = 0.0;
			for (int k=0; k<GetNcat(); k++)	{
				part += rate[p][k];
			}
			part /= GetNcat();

			total += part * GetPartNsite(p) * ratemult[p];
		}
		return total / GetNsite();
	}

	// uses suffisicent stats
	double Move(double tuning, int nrep);
	double NonMPIMove(double tuning, int nrep);

	double MoveHyper(double tuning, int nrep);

	// uses suffisicent stats
	int MoveAlphas(double tuning);
	void MoveAlphaHyper();
	
	// uses suffisicent stats
	void MoveMultipliers();
	int MoveMultiplierHyper(double tuning);

	void SetAlpha(int inpart, double inalpha)	{
		alpha[inpart] = inalpha;
		UpdateDiscreteCategories(inpart);
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	void Create(int incat, PartitionScheme inscheme);
	void Delete();

	void SampleRate();

	double LogRatePrior();
	double LogAlphaPrior(int inpart);
	double LogRateLikelihood(int inpart);

	double LogMultiplierPrior();

	void GlobalUpdateRateSuffStat();
	void SlaveUpdateRateSuffStat();
	void UpdateRateSuffStat();

	void UpdateDiscreteCategories(int inpart);

	int Ncat;
	double** rate;
	int* alloc;
	int** ratesuffstatcount;
	double** ratesuffstatbeta;

	double* alpha;
	double alphaHyper;

	double* ratemult;
	double multHyper;
};

#endif

