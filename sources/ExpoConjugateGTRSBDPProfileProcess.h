
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef EXPCONGTRSBDPPROFILE_H
#define EXPCONGTRSBDPPROFILE_H

#include "MatrixSBDPProfileProcess.h"
#include "ExpoConjugateGTRMixtureProfileProcess.h"

// Exponential conjugate GTR models
class ExpoConjugateGTRSBDPProfileProcess : public virtual MatrixSBDPProfileProcess, public virtual ExpoConjugateGTRMixtureProfileProcess {

	public:

	ExpoConjugateGTRSBDPProfileProcess() : InitIncremental(0) {}
	virtual ~ExpoConjugateGTRSBDPProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();
		// GlobalUpdateParameters();
		// GlobalUpdateSiteProfileSuffStat();
		// UpdateModeProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{

			// relative rates
			if (! fixrr)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				MoveRR();
			}

			incchrono.Start();
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			if ((!rep) && InitIncremental)	{
				cerr << "init incremental\n";
				InitIncremental--;
				IncrementalSampleAlloc();
				UpdateModeProfileSuffStat();
			}

			GlobalMixMove(5,1,0.001,40);
			// MixMove(5,1,0.001,40);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);
			incchrono.Stop();

			/*
			// allocations
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			// no parallel version of the allocation move for the moment
			// IncrementalDPMove(5,2);
			incchrono.Start();
			// ResampleEmptyProfiles();
			// GlobalIncrementalDPMove(25,3);
			// GlobalIncrementalDPMove(5,0.5);
			GlobalIncrementalDPMove(5,0.005);
			// IncrementalDPMove(5,3);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);
			incchrono.Stop();

			// GlobalIncrementalSBDPMove(5);

			// profiles
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			profilechrono.Start();
			GlobalMoveProfile(1,1,100);
			GlobalMoveProfile(1,3,100);
			GlobalMoveProfile(0.1,3,100);
			profilechrono.Stop();
			// MoveProfile(1,1,100);
			// MoveProfile(1,3,100);
			// MoveProfile(0.1,3,100);
			*/

			// hyperparameters
			// globalupdates useless as long as not relying on sufficient statistics
			// themselves dependent on new parameter values
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			MoveHyper(tuning,10);
		}
		totchrono.Stop();
		return 1;
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual double LogProxy(int site, int cat)	{
		return PoissonDiffLogSampling(cat,site);
	}


	virtual void SwapComponents(int cat1, int cat2)	{
		MatrixSBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	virtual void Create(int innsite, int indim)	{
		MatrixSBDPProfileProcess::Create(innsite,indim);
		ExpoConjugateGTRMixtureProfileProcess::Create(innsite,indim);
	}

	virtual void Delete()	{
		MatrixSBDPProfileProcess::Delete();
		ExpoConjugateGTRMixtureProfileProcess::Delete();
	}

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;

	int InitIncremental;
};

#endif

