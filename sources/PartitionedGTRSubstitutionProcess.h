
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
		sitemask = vector<bool>(insitemax - insitemin, false);
	}

	void Delete() {
		SubstitutionProcess::Delete();
		PartitionedGTRProfileProcess::Delete();
	}

	int GetNstate(int site) {return GetDim();}
	int GetNstate() {return GetDim();}

	// CPU : level 1
	void Reset(double*** condl, bool condalloc = false);
	void Multiply(double*** from, double*** to, bool condalloc = false);
	void MultiplyByStationaries(double*** from, bool condalloc = false);
	void Offset(double*** condl, bool condalloc = false);
	void Initialize(double*** condl, const int* leafstates, bool condalloc = false);

	// CPU : level 2
	double ComputeLikelihood(double*** aux, bool condalloc = false);

	// CPU : level 3
	// implemented in GTR or POisson Substitution process
	void Propagate(double*** from, double*** to, double time, bool condalloc = false);

	/*
	double* GetStationary(int site)	{
		return GetProfile(site);
	}
	*/

	vector<bool> sitemask;
};

#endif

