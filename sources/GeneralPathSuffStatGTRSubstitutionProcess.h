
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSGTRSUB_H
#define GENPATHSSGTRSUB_H

#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"
#include "GTRSubstitutionProcess.h"
#include "GeneralPathSuffStatGTRProfileProcess.h"
#include "RateProcess.h"

class GeneralPathSuffStatGTRSubstitutionProcess : public virtual GeneralPathSuffStatMatrixSubstitutionProcess, public virtual GTRSubstitutionProcess, public virtual GeneralPathSuffStatGTRProfileProcess, public virtual RateProcess	{

	public:

	GeneralPathSuffStatGTRSubstitutionProcess() {}
	virtual ~GeneralPathSuffStatGTRSubstitutionProcess() {}

	protected:

	void Create(int innsite, int indim,int insitemin,int insitemax)	{
		GeneralPathSuffStatGTRProfileProcess::Create(innsite,indim);
		SubstitutionProcess::Create(innsite,indim,insitemin,insitemax);
	}

	void Delete() {
		SubstitutionProcess::Delete();
		GeneralPathSuffStatGTRProfileProcess::Delete();
	}
};

#endif

