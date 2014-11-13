
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSELFINITESUB_H
#define CODONMUTSELFINITESUB_H

#include "CodonMutSelFiniteProfileProcess.h"
#include "UniformRateProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"


class CodonMutSelFiniteSubstitutionProcess : public virtual CodonMutSelFiniteProfileProcess, public virtual UniformRateProcess, public virtual GeneralPathSuffStatMatrixSubstitutionProcess {

	// s'inspirer de GeneralPathSuffStatGTRSubstitutionProcess
	// et GeneralPathSuffStatRASCATGTRSubstitutionProcess

	public:

	CodonMutSelFiniteSubstitutionProcess() {}
	virtual ~CodonMutSelFiniteSubstitutionProcess() {}

	protected:

	void Create(int insite, int indim, int sitemin, int sitemax)	{
		cerr << "In four-argument Create of CodonMutSelFiniteSubstitutionProcess. Should not be here.\n";
		exit(1);
	}

	void Create(int innsite, int indim, int ncat, int infixncomp, int inempmix, string inmixtype, int sitemin, int sitemax, CodonStateSpace* instatespace)	{
		if (ncat == -1)	{
			ncat = innsite;
		}
	
		CodonMutSelFiniteProfileProcess::Create(innsite,indim,ncat,infixncomp,inempmix,inmixtype,instatespace);
		UniformRateProcess::Create(innsite);
		GeneralPathSuffStatMatrixSubstitutionProcess::Create(innsite,indim,sitemin,sitemax);
	}

	void Delete()	{
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
		UniformRateProcess::Delete();
		CodonMutSelFiniteProfileProcess::Delete();
	}

};

#endif

