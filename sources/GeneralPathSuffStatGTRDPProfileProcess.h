
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSGTRDPPROFILE_H
#define GENPATHSSGTRDPPROFILE_H

#include "MatrixDPProfileProcess.h"
#include "GTRMixtureProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

// Exponential conjugate GTR models
// class GeneralPathSuffStatGTRDPProfileProcess : public virtual MatrixDPProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess {
class GeneralPathSuffStatGTRDPProfileProcess : public virtual MatrixDPProfileProcess, public virtual GTRMixtureProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess {

	public:

	GeneralPathSuffStatGTRDPProfileProcess() {}
	virtual ~GeneralPathSuffStatGTRDPProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();

		for (int rep=0; rep<nrep; rep++)	{

			// relative rates
			if (! fixrr)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				UpdateModeProfileSuffStat();
				MoveRR();
			}

			/*
			incchrono.Start();
			GlobalMixMove(1,100);
			incchrono.Stop();
			*/

			// allocations
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			// no parallel version of the allocation move for the moment
			incchrono.Start();
			IncrementalDPMove(5);
			incchrono.Stop();
			// profiles
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			profilechrono.Start();
			GlobalMoveProfile(1,1,100);
			GlobalMoveProfile(1,3,100);
			GlobalMoveProfile(0.1,3,100);
			// MoveProfile(1,1,100);
			// MoveProfile(1,3,100);
			// MoveProfile(0.1,3,100);
			profilechrono.Stop();

			// CreateMatrices();
			UpdateMatrices();
			// Phyperparameters
			// globalupdates useless as long as not relying on sufficient statistics
			// themselves dependent on new parameter values
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			// UpdateModeProfileSuffStat();
			MoveHyper(tuning,10);

		}
		totchrono.Stop();
		return 1;
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	void Create(int innsite, int indim)	{
		MatrixDPProfileProcess::Create(innsite,indim);
		GTRMixtureProfileProcess::Create(innsite,indim);
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create(innsite,indim);
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		MatrixDPProfileProcess::Delete();
		GTRMixtureProfileProcess::Delete();
	}

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;
};


#endif


