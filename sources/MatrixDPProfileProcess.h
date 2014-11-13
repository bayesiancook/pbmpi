
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXDPPROFILE_H
#define MATRIXDPPROFILE_H

#include "MatrixMixtureProfileProcess.h"
#include "DPProfileProcess.h"

class MatrixDPProfileProcess : public virtual MatrixMixtureProfileProcess, public virtual DPProfileProcess {

	public:

	MatrixDPProfileProcess():  Nadd(30) , Ninc(3) {}
	virtual ~MatrixDPProfileProcess() {}

	protected:

	virtual void Create(int innsite, int indim)	{
		DPProfileProcess::Create(innsite,indim);
		MatrixMixtureProfileProcess::Create(innsite,indim);
	}

	virtual void Delete()	{
		MatrixMixtureProfileProcess::Delete();
		DPProfileProcess::Delete();
	}

	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	// double GlobalIncrementalDPMove(int nrep);
	// double SlaveIncrementalDPMove();
	double IncrementalDPMove(int nrep);
	double GlobalMixMove(int Nrep, int Nprofile);
	void SlaveMixMove();

	// SubMatrix** matrixarray;
	// in DP incremental move: Nadd additional components are proposed
	int Nadd;
	int Ninc;
};

#endif
