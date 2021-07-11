
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef POISSONFINITEPROFILE_H
#define POISSONFINITEPROFILE_H

#include "PoissonMixtureProfileProcess.h"
#include "FiniteProfileProcess.h"

// superclass for Poisson (F81) implementations
class PoissonFiniteProfileProcess: public virtual PoissonMixtureProfileProcess, public virtual FiniteProfileProcess	{

	public:

	PoissonFiniteProfileProcess() {}
	virtual ~PoissonFiniteProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			GlobalIncrementalFiniteMove(1);

			if (! empmix)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				UpdateModeProfileSuffStat();
				MoveProfile();
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				MoveHyper(tuning,10);
			}

			if (! fixncomp)	{
				MoveNcomponent(100);
            }
		}
		return 1;
	}

	virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is);

	protected:

	virtual void Create(int innsite, int indim)	{
		cerr << "in create 2 arguments\n";
		exit(1);
	}

	virtual void Create(int innsite, int indim, int ncat, int infixncomp, int inempmix, string inmixtype)	{
		FiniteProfileProcess::Create(innsite,indim,ncat,infixncomp,inempmix,inmixtype);
		PoissonMixtureProfileProcess::Create(innsite,indim);
	}

	virtual void Delete()	{
		PoissonMixtureProfileProcess::Delete();
		FiniteProfileProcess::Delete();
	}

	double IncrementalFiniteMove(int nrep)	{
		cerr << "error : in Poisson Finite Profile Process inc move\n";
		exit(1);
	}

	double GlobalIncrementalFiniteMove(int nrep);
	double SlaveIncrementalFiniteMove();


};

#endif

