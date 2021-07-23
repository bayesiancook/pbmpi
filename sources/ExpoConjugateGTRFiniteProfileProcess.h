
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef EXPCONGTRFINITEPROFILE_H
#define EXPCONGTRFINITEPROFILE_H

#include "MatrixFiniteProfileProcess.h"
#include "ExpoConjugateGTRMixtureProfileProcess.h"

// Exponential conjugate GTR models
class ExpoConjugateGTRFiniteProfileProcess : public virtual MatrixFiniteProfileProcess, public virtual ExpoConjugateGTRMixtureProfileProcess {

	public:

	ExpoConjugateGTRFiniteProfileProcess() {}
	virtual ~ExpoConjugateGTRFiniteProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();

		for (int rep=0; rep<nrep; rep++)	{

			// relative rates
			if (! fixrr)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				// useless in an expo suff stat context
				// UpdateModeProfileSuffStat();
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
			incchrono.Stop();

			if (! empmix)	{
				// profiles
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				UpdateModeProfileSuffStat();
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
				MoveHyper(tuning,10);
			}

			// hyperparameters
			// globalupdates useless as long as not relying on sufficient statistics
			// themselves dependent on new parameter values
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();

			if (! fixncomp)	{
				MoveNcomponent(100);
				// MoveWeightAlpha(tuning,10);
			}
			// MoveWeightAlpha(tuning,10);
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
		ExpoConjugateGTRMixtureProfileProcess::Create(innsite,indim);
		SetRR(rrtype);
	}

	virtual void Delete()	{
		MatrixFiniteProfileProcess::Delete();
		ExpoConjugateGTRMixtureProfileProcess::Delete();
	}

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;
};

#endif

