
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTEXPCONGTRPROFILE_H
#define PARTEXPCONGTRPROFILE_H

#include "PartitionedGTRProfileProcess.h"

// superclass for GTR-like models using the exponential conjugate relation between relative rates and substitution processes
class PartitionedExpoConjugateGTRProfileProcess : public virtual PartitionedGTRProfileProcess {

	public:

	PartitionedExpoConjugateGTRProfileProcess() : rrsuffstatcount(0), rrsuffstatbeta(0) {}
	virtual ~PartitionedExpoConjugateGTRProfileProcess() {}

	// protected:

	virtual void Create(int indim, PartitionScheme inscheme);
	virtual void Delete();

	// profiles
	virtual const int* GetSiteProfileSuffStatCount(int site) = 0;
	virtual const double* GetSiteProfileSuffStatBeta(int site) = 0;

	// update of relative rates
	// conjugate Gibbs resampling
	void MoveRR();

	int** rrsuffstatcount;
	double** rrsuffstatbeta;
};

#endif

