
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RATE_H
#define RATE_H

#include <iostream>
using namespace std;

#include "Chrono.h"

class RateProcess {

	public:

	RateProcess() : nsite(0), ratefrac(1.0) {}
	virtual ~RateProcess() {}

	virtual string GetVersion() = 0;

	int GetNsite() {return nsite;}

	virtual int GetNrate(int site) = 0;
    virtual int GetMaxNrate() = 0;
	virtual double GetRate(int site, int cat = 0) = 0;
	virtual double GetRateWeight(int site, int cat) = 0;
	double GetMeanRate();
	virtual double GetPriorMeanRate() = 0;
	virtual double GetAlpha() {return 1;}

	virtual void SiteActivateSumOverRateAllocation(int site)	{
		cerr << "in RateProcess::SiteActivateSumOverRateAlloc\n";
		exit(1);
	}
	virtual void SiteInactivateSumOverRateAllocation(int site)	{
		cerr << "in RateProcess::SiteINactivateSumOverRateAlloc\n";
		exit(1);
	}
	virtual void ActivateSumOverRateAllocations() = 0;
	virtual void InactivateSumOverRateAllocations(int* ratealloc) = 0;
	bool SumOverRateAllocations() {return sumflag;}

	virtual double LogRatePrior() = 0;
	virtual void SampleRate() = 0;
    virtual void PriorSampleRate() = 0;

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;

	protected:

    void SetRateFrac(double infrac) {
        ratefrac = infrac;
    }

	// abstract classes will be implemented in phyloprocess
	virtual void GlobalUpdateSiteRateSuffStat() = 0;
	virtual void SlaveUpdateSiteRateSuffStat() = 0;

	virtual void UpdateSiteRateSuffStat() = 0;
	virtual double GetSiteRateSuffStatBeta(int site) = 0;
	virtual int GetSiteRateSuffStatCount(int site) = 0;

	void Create(int innsite)	{
		nsite = innsite;
	}
	void Delete() {}

	virtual int GetNprocs() = 0;
	virtual int GetMyid() = 0;
	virtual int GetSiteMin() = 0;
	virtual int GetSiteMax() = 0;

	bool sumflag;
	int nsite;

    double ratefrac;

	Chrono chronorate;
};

#endif

