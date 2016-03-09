
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRFINITEPROFILE_H
#define PARTGTRFINITEPROFILE_H

#include "PartitionedGTRMixtureProfileProcess.h"
#include "FiniteProfileProcess.h"

class PartitionedGTRFiniteProfileProcess : public virtual PartitionedGTRMixtureProfileProcess, public virtual FiniteProfileProcess {

	public:

	PartitionedGTRFiniteProfileProcess() {}
	virtual ~PartitionedGTRFiniteProfileProcess() {}

	protected:

	virtual void Create(int innsite, int indim)	{
		cerr << "in create 2 arguments\n";
		exit(1);
	}

	virtual void Create(int indim, PartitionScheme rrscheme, int ncat, int infixncomp = 0, int inempmix = 0, string inmixtype = "None")	{
		FiniteProfileProcess::Create(rrscheme.GetNsite(),indim,ncat,infixncomp,inempmix,inmixtype);
		PartitionedGTRMixtureProfileProcess::Create(indim,rrscheme);
	}

	virtual void Delete()	{
		FiniteProfileProcess::Delete();
		PartitionedGTRMixtureProfileProcess::Delete();
	}

	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	double GlobalIncrementalFiniteMove(int nrep);
	double SlaveIncrementalFiniteMove();
	double IncrementalFiniteMove(int nrep);

	// SubMatrix** matrixarray;
};

#endif
