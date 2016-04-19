
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef DPPROFILE_H
#define DPPROFILE_H

#include <cmath>
#include "MixtureProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class DPProfileProcess: public virtual MixtureProfileProcess	{

	public:

	DPProfileProcess() : kappa(1), movekappa(true), kappaprior(0) {}
	virtual ~DPProfileProcess(){}

	protected:

	virtual double IncrementalDPMove(int nrep) = 0;
	double MoveHyper(double tuning, int nrep);
	virtual double MoveKappa(double tuning, int nrep);


	// static allocation of many component-specific variables
	// such as: profiles, occupancy number
	// basically everything except substitution matrices

	// called at the beginning and end of the run (see PhyloProcess)
	/*
	virtual void Create(int innsite, int indim);
	virtual void Delete();
	*/

	// multinomial 
	virtual double LogProxy(int site, int cat);
	virtual void SampleAlloc();
	void SampleHyper();

	// kappa has an exponential prior of mean 10
	double LogHyperPrior();
	virtual double LogAllocPrior();

	double kappa;
	bool movekappa;
	int kappaprior;
	// 0 : exponential of mean 20
	// 1 : jeffreys prior 
	int dirweightprior;
	// 0 : flexible
	// 1 : rigid
};

#endif

