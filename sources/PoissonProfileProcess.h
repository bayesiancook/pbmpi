
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef POISSONPROFILE_H
#define POISSONPROFILE_H

#include "ProfileProcess.h"

// Poisson (F81) implementation
class PoissonProfileProcess : public virtual ProfileProcess {

	public:

	PoissonProfileProcess() {}
	virtual ~PoissonProfileProcess() {}

	protected:

	virtual void UpdateZip(int site) = 0;

	// implemented in specialized phyloprocess subclasses
	virtual const int* GetSiteProfileSuffStatCount(int site) = 0;

};

#endif

