
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRSUB_H
#define PARTGTRSUB_H

#include "MatrixSubstitutionProcess.h"
#include "PartitionedGTRProfileProcess.h"
#include "RateProcess.h"

class PartitionedGTRSubstitutionProcess : public virtual MatrixSubstitutionProcess, public virtual PartitionedGTRProfileProcess, public virtual RateProcess	{

	public:

	PartitionedGTRSubstitutionProcess() {}
	virtual ~PartitionedGTRSubstitutionProcess() {}

	protected:

	void Create(int indim, PartitionScheme rrscheme, int insitemin,int insitemax)	{
		PartitionedGTRProfileProcess::Create(indim, rrscheme);
		SubstitutionProcess::Create(rrscheme.GetNsite(),indim,insitemin,insitemax);
	}

	void Delete() {
		SubstitutionProcess::Delete();
		PartitionedGTRProfileProcess::Delete();
	}

	int GetNstate(int site) {return GetDim();}
	int GetNstate() {return GetDim();}

	/*
	double* GetStationary(int site)	{
		return GetProfile(site);
	}
	*/
};

#endif

