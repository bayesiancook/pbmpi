
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRPHYLO_H
#define PARTGTRPHYLO_H


#include "PhyloProcess.h"
#include "PartitionedGTRSubstitutionProcess.h"


class PartitionedGTRPhyloProcess : public virtual PhyloProcess, public virtual PartitionedGTRSubstitutionProcess	{


	public:

	PartitionedGTRPhyloProcess() {}
	virtual ~PartitionedGTRPhyloProcess() {}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata,int indim, PartitionScheme rrscheme, int insitemin,int insitemax)	{
		PhyloProcess::Create(intree,indata,indata->GetNstate());
		PartitionedGTRSubstitutionProcess::Create(indata->GetNstate(),rrscheme, insitemin,insitemax);
	}

	void Delete() {
		PartitionedGTRSubstitutionProcess::Delete();
		PhyloProcess::Delete();
	}

	// override PhyloProcess functions because should also create and delete all the matrices
	virtual void Unfold();
	virtual void Collapse();


	double LengthRelRateMove(double tuning, int nrep);
};

#endif

