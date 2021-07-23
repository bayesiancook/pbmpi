
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef SBDPPROFILE_H
#define SBDPPROFILE_H

#include <cmath>
#include "DPProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class SBDPProfileProcess: public virtual DPProfileProcess	{

	using MixtureProfileProcess::LogStatPrior;

	public:

	SBDPProfileProcess() : DPProfileProcess(), nmodemax(refnmodemax), V(0), maxweighterror(0) {}
	virtual ~SBDPProfileProcess(){}

	protected:

	virtual void DrawProfileFromPrior();

	double GetMaxWeightError() {return maxweighterror;}
	void ResetMaxWeightError() {maxweighterror = 0;}

	virtual void Create(int innsite, int indim);
	virtual void Delete();

	virtual int GetNmodeMax() {return GetNsite() > nmodemax ? nmodemax : GetNsite();}
	virtual void SetNmodeMax(int n) {nmodemax = n;}

	virtual double IncrementalDPMove(int nrep, double epsilon) = 0;

	double IncrementalDPMove(int nrep)	{
		cerr << "inc move deactivated\n";
		exit(1);
	}

	double GetWeightEnt()	{
		double tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (weight[k] > 1e-8)	{
				tot -= weight[k] * log(weight[k]);
			}
		}
		return tot;
	}

	int GetNDisplayedComponent()	{
		return GetNOccupiedComponent();
	}

	int GetLastOccupiedComponent()	{
		int kmax = 0;
		for (int i=0; i<GetNsite(); i++)	{
			if (kmax < alloc[i])	{
				kmax = alloc[i];
			}
		}
		return kmax;
	}

	int GetNCutoff(double cutoff)	{
		int n = (int) (GetNOccupiedComponent() * (1 - cutoff));
		int k = GetLastOccupiedComponent();
		int tot = occupancy[k];
		while (k && (tot < n))	{
			k--;
			if (k)	{
				tot += occupancy[k];
			}
		}
		return k;
	}
		
	virtual void SwapComponents(int cat1, int cat2);

	// double LogAllocPrior();
	double LogIntegratedAllocProb();
	double MoveKappa(double tuning, int nrep);

	// void ShedTail();

	// redefined
	void SampleAlloc();

	void IncrementalSampleAlloc();

	void SampleWeights();
	void ResampleWeights();
	// void ResampleLastWeight();
	double MoveOccupiedCompAlloc(int nrep = 1);
	double MoveAdjacentCompAlloc(int nrep = 1);

	double LogStatPrior();

	/*
	double totweight;
	double cumulProduct;
	*/

	int nmodemax;
	double* V;
	double* weight;

	double maxweighterror;
};

#endif

