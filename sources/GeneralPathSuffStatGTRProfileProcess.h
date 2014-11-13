
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSGTRPROFILE_H
#define GENPATHSSGTRPROFILE_H

#include "GTRProfileProcess.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"

// superclass for GTR-like models using the general sufficient statistics
class GeneralPathSuffStatGTRProfileProcess : public virtual GTRProfileProcess, public virtual GeneralPathSuffStatMatrixProfileProcess {

	public:

	GeneralPathSuffStatGTRProfileProcess() {}
	virtual ~GeneralPathSuffStatGTRProfileProcess() {}

	protected:

	// update of relative rates
	// Metropolis Hastings using global suff stats
	void MoveRR()	{
		MoveRR(1.0);
		MoveRR(0.1);
	}

	double MoveRR(double tuning);
};

#endif

