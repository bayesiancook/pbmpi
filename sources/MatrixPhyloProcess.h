
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXPHYLO_H
#define MATRIXPHYLO_H


#include "PhyloProcess.h"
#include "MatrixSubstitutionProcess.h"
#include "GTRSubstitutionProcess.h"

class MatrixPhyloProcess : public virtual PhyloProcess, public virtual MatrixSubstitutionProcess	{

	public:

	MatrixPhyloProcess() {}
	virtual ~MatrixPhyloProcess() {}

	// protected:

	// override PhyloProcess functions because should also create and delete all the matrices
	virtual void Unfold();
	virtual void Collapse();

	virtual void UpdateConditionalLikelihoods();

	// virtual void UpdateSubstitutionProcess();

};


class GTRPhyloProcess : public virtual MatrixPhyloProcess, public virtual GTRSubstitutionProcess	{


	public:

	GTRPhyloProcess() {}
	virtual ~GTRPhyloProcess() {}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata,int indim,int insitemin,int insitemax)	{
		MatrixPhyloProcess::Create(intree,indata,indata->GetNstate());
		GTRSubstitutionProcess::Create(indata->GetNsite(),indata->GetNstate(),insitemin,insitemax);
	}

	void Delete() {
		GTRSubstitutionProcess::Delete();
		MatrixPhyloProcess::Delete();
	}


	double LengthRelRateMove(double tuning, int nrep);
};

#endif

