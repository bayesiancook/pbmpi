
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef PARTEXPCONGTRSBDPPROFILE_H
#define PARTEXPCONGTRSBDPPROFILE_H

#include "PartitionedGTRSBDPProfileProcess.h"
#include "PartitionedExpoConjugateGTRMixtureProfileProcess.h"

// Exponential conjugate GTR models
class PartitionedExpoConjugateGTRSBDPProfileProcess : public virtual PartitionedGTRSBDPProfileProcess, public virtual PartitionedExpoConjugateGTRMixtureProfileProcess {

	public:

	PartitionedExpoConjugateGTRSBDPProfileProcess() {}
	virtual ~PartitionedExpoConjugateGTRSBDPProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();
		// GlobalUpdateParameters();
		// GlobalUpdateSiteProfileSuffStat();
		// UpdateModeProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{

			// relative rates
			if (nfreerr > 0)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				MoveRR();
			}

			incchrono.Start();
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			GlobalMixMove(5,1,0.001,40);
			// MixMove(5,1,0.001,40);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);
			incchrono.Stop();

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

	virtual void SwapComponents(int cat1, int cat2)	{
		PartitionedGTRSBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	virtual void Create(int indim, PartitionScheme rrscheme)	{
		PartitionedGTRSBDPProfileProcess::Create(indim, rrscheme);
		PartitionedExpoConjugateGTRMixtureProfileProcess::Create(indim, rrscheme);
	}

	virtual void Delete()	{
		PartitionedGTRSBDPProfileProcess::Delete();
		PartitionedExpoConjugateGTRMixtureProfileProcess::Delete();
	}

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;
};

#endif

