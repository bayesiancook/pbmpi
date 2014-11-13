
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef POISSONDPPROFILE_H
#define POISSONDPPROFILE_H

#include "PoissonMixtureProfileProcess.h"
#include "DPProfileProcess.h"

// superclass for Poisson (F81) implementations
class PoissonDPProfileProcess: public virtual PoissonMixtureProfileProcess, public virtual DPProfileProcess	{

	public:

	PoissonDPProfileProcess() {}
	virtual ~PoissonDPProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			IncrementalDPMove(2);
			MoveProfile();
			MoveHyper(tuning,10);
			MoveHyper(0.1 * tuning,10);
		}
		return 1;
	/*
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			IncrementalDPMove(nrep);
			MoveProfile();
			MoveHyper(tuning,nrep);
	*/
	}

	virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is);

	protected:

	/*
	virtual void Create(int innsite, int indim)	{
		DPProfileProcess::Create(innsite,indim);
		PoissonMixtureProfileProcess::Create(innsite,indim);
	}

	virtual void Delete()	{
		DPProfileProcess::Delete();
		PoissonMixtureProfileProcess::Delete();
	}
	*/

	double IncrementalDPMove(int nrep);

};

#endif

