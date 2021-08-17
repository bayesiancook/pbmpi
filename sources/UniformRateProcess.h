
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef UNIRATE_H
#define UNIRATE_H

#include "RateProcess.h"

class UniformRateProcess : public virtual RateProcess {

	public:

	UniformRateProcess() {}
	~UniformRateProcess() {}

    int GetMaxNrate() {return 1;}
	int GetNrate(int site) {return 1;}
	double GetRate(int site, int cat = 0)	{return 1.0;}
	double GetRateWeight(int site, int cat) {return 1.0;}

	double GetPriorMeanRate() {return 1;}

	void ActivateSumOverRateAllocations() {
		sumflag = true;
	}

	void InactivateSumOverRateAllocations(int* ratealloc) {
		sumflag = false;
	}

	void SiteActivateSumOverRateAllocation(int site) {
		sumflag = true;
	}

	void SiteInactivateSumOverRateAllocation(int site)	{
		sumflag = false;
	}

	virtual double LogRatePrior() {return 0;}
	virtual void SampleRate() {}
    virtual void PriorSampleRate() {}

	void ToStream(ostream& os) {}
	void FromStream(istream& is) {}

};

#endif

