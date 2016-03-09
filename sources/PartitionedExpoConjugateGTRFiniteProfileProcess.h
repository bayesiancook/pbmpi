
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef PARTEXPCONGTRFINITEPROFILE_H
#define PARTEXPCONGTRFINITEPROFILE_H

#include "PartitionedGTRFiniteProfileProcess.h"
#include "PartitionedExpoConjugateGTRMixtureProfileProcess.h"

// Exponential conjugate GTR models
class PartitionedExpoConjugateGTRFiniteProfileProcess : public virtual PartitionedGTRFiniteProfileProcess, public virtual PartitionedExpoConjugateGTRMixtureProfileProcess {

	public:

	PartitionedExpoConjugateGTRFiniteProfileProcess() {}
	virtual ~PartitionedExpoConjugateGTRFiniteProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();

		for (int rep=0; rep<nrep; rep++)	{

			// relative rates
			if (nfreerr > 0)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				// useless in an expo suff stat context
				// UpdateModeProfileSuffStat();
				MoveRR();
			}

			// allocations
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			// no parallel version of the allocation move for the moment
			incchrono.Start();
			// IncrementalFiniteMove(1);
			GlobalIncrementalFiniteMove(1);
			incchrono.Stop();

			if (! empmix)	{
				// profiles
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				UpdateModeProfileSuffStat();
				profilechrono.Start();
				GlobalMoveProfile(1,1,100);
				GlobalMoveProfile(1,3,100);
				GlobalMoveProfile(0.1,3,100);
				/*
				MoveProfile(1,1,100);
				MoveProfile(1,3,100);
				MoveProfile(0.1,3,100);
				*/
				profilechrono.Stop();
				MoveHyper(tuning,10);
			}

			if (! fixncomp)	{
				MoveNcomponent(10);
				MoveWeightAlpha(tuning,10);
			}
		}
		totchrono.Stop();
		return 1;
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual void Create(int innsite, int indim)	{
		cerr << "in create 2 arguments\n";
		exit(1);
	}

	virtual void Create(int indim, PartitionScheme rrscheme, int ncat, int infixncomp, int inempmix, string inmixtype)	{
		PartitionedGTRFiniteProfileProcess::Create(indim,rrscheme,ncat,infixncomp,inempmix,inmixtype);
		PartitionedExpoConjugateGTRMixtureProfileProcess::Create(indim,rrscheme);
	}

	virtual void Delete()	{
		PartitionedGTRFiniteProfileProcess::Delete();
		PartitionedExpoConjugateGTRMixtureProfileProcess::Delete();
	}

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;
};

#endif

