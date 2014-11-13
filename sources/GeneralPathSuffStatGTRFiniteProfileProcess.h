
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSGTRFINITEPROFILE_H
#define GENPATHSSGTRFINITEPROFILE_H

#include "MatrixFiniteProfileProcess.h"
#include "GTRMixtureProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

// Exponential conjugate GTR models
// class GeneralPathSuffStatGTRFiniteProfileProcess : public virtual MatrixFiniteProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess {
class GeneralPathSuffStatGTRFiniteProfileProcess : public virtual MatrixFiniteProfileProcess, public virtual GTRMixtureProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess {

	public:

	GeneralPathSuffStatGTRFiniteProfileProcess() {}
	virtual ~GeneralPathSuffStatGTRFiniteProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();

		for (int rep=0; rep<nrep; rep++)	{

			// UpdateMatrices();
			if (! fixrr)	{
				// relative rates
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				UpdateModeProfileSuffStat();
				MoveRR();
			}

			// allocations
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			// no parallel version of the allocation move for the moment
			incchrono.Start();
			// IncrementalFiniteMove(1);
			GlobalIncrementalFiniteMove(1);
			// MoveNcomponent(10);
			incchrono.Stop();

			if (! empmix)	{
			// profiles
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				// UpdateModeProfileSuffStat();
				profilechrono.Start();
				GlobalMoveProfile(1,1,100);
				GlobalMoveProfile(1,3,100);
				GlobalMoveProfile(0.1,3,100);
				/*
				MoveProfile(1,1,100);
				MoveProfile(1,3,100);
				MoveProfile(0.1,3,100);
				*/
				profilechrono.Stop();
			}

			// hyperparameters
			// globalupdates useless as long as not relying on sufficient statistics
			// themselves dependent on new parameter values
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			MoveHyper(tuning,10);

			if (! fixncomp)	{
				MoveNcomponent(10);
			}
			MoveWeightAlpha(tuning,10);

			// important since matrices have changed (through profile and alloc moves)
			// and matrices allocated by MASTER are used in MoveRR
			UpdateMatrices();

		}
		totchrono.Stop();
		return 1;
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual void Create(int innsite, int indim)	{
		cerr << "in create 2 arguments\n";
		exit(1);
	}

	virtual void Create(int innsite, int indim, int ncat, int infixncomp, int inempmix, string inmixtype,string rrtype)	{
		MatrixFiniteProfileProcess::Create(innsite,indim,ncat,infixncomp,inempmix,inmixtype);
		GTRMixtureProfileProcess::Create(innsite,indim);
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create(innsite,indim);
		SetRR(rrtype);
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		MatrixFiniteProfileProcess::Delete();
		GTRMixtureProfileProcess::Delete();
	}

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;
};


#endif


