
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GTRSUB_H
#define GTRSUB_H

#include "MatrixSubstitutionProcess.h"
#include "GTRProfileProcess.h"
#include "RateProcess.h"

class GTRSubstitutionProcess : public virtual MatrixSubstitutionProcess, public virtual GTRProfileProcess, public virtual RateProcess	{

	public:

	GTRSubstitutionProcess() {}
	virtual ~GTRSubstitutionProcess() {}

	protected:

	void Create(int innsite, int indim,int insitemin,int insitemax)	{
		GTRProfileProcess::Create(innsite,indim);
		SubstitutionProcess::Create(innsite,indim,insitemin,insitemax);
	}

	void Delete() {
		SubstitutionProcess::Delete();
		GTRProfileProcess::Delete();
	}

	// int GetNstate(int site) {return GetDim();}
	// int GetNstate() {return GetDim();}

	/*
	double* GetStationary(int site)	{
		return GetProfile(site);
	}
	*/
};

#endif

