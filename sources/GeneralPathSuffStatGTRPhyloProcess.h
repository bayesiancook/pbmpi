
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSGTRPHYLO_H
#define GENPATHSSGTRPHYLO_H

#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GeneralPathSuffStatGTRSubstitutionProcess.h"
#include "RateProcess.h"

// this is the General Suff Stat version
class GeneralPathSuffStatGTRPhyloProcess : public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GTRPhyloProcess, public virtual GeneralPathSuffStatGTRSubstitutionProcess, public virtual RateProcess {

	public:

	GeneralPathSuffStatGTRPhyloProcess() {}
	virtual ~GeneralPathSuffStatGTRPhyloProcess() {}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int indim, int sitemin, int sitemax)	{
		RateProcess::Create(indata->GetNsite());
		GeneralPathSuffStatGTRSubstitutionProcess::Create(indata->GetNsite(),indim,sitemin,sitemax);
		GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,indim,sitemin,sitemax);
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		GeneralPathSuffStatGTRSubstitutionProcess::Delete();
		RateProcess::Delete();
	}

	void UpdateRRSuffStat() {}
};

#endif

