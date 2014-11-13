
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSBDPSUB_H
#define AACODONMUTSELSBDPSUB_H

#include "AACodonMutSelGeneProfileProcess.h"
#include "UniformRateProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"


class AACodonMutSelGeneSubstitutionProcess : public virtual AACodonMutSelGeneProfileProcess, public virtual UniformRateProcess, public virtual GeneralPathSuffStatMatrixSubstitutionProcess {

	// s'inspirer de GeneralPathSuffStatGTRSubstitutionProcess
	// et GeneralPathSuffStatRASCATGTRSubstitutionProcess

	public:

	AACodonMutSelGeneSubstitutionProcess() {}
	virtual ~AACodonMutSelGeneSubstitutionProcess() {}

	protected:

	void Create(int insite, int indim, int sitemin, int sitemax)	{
		cerr << "In four-argument Create of AACodonMutSelGeneSubstitutionProcess. Should not be here.\n";
		exit(1);
	}

	void Create(int innsite, int indim, CodonStateSpace* instatespace, int infixbl, int infixcodonprofile, int infixomega, AACodonMutSelProfileProcess* hub, int offset)	{
		AACodonMutSelGeneProfileProcess::Create(innsite,indim,instatespace,infixbl,infixcodonprofile,infixomega,hub,offset);
		UniformRateProcess::Create(innsite);
		GeneralPathSuffStatMatrixSubstitutionProcess::Create(innsite,indim,0,innsite);
	}

	void Delete()	{
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
		UniformRateProcess::Delete();
		AACodonMutSelGeneProfileProcess::Delete();
	}

};

#endif

