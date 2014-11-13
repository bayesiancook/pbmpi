
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXSBDPPROFILE_H
#define MATRIXSBDPPROFILE_H

#include "MatrixMixtureProfileProcess.h"
#include "SBDPProfileProcess.h"

class MatrixSBDPProfileProcess : public virtual MatrixMixtureProfileProcess, public virtual SBDPProfileProcess {

	public:

	MatrixSBDPProfileProcess() {}
	virtual ~MatrixSBDPProfileProcess() {}

	protected:

	virtual void SwapComponents(int cat1, int cat2);

	virtual void Create(int innsite, int indim)	{
		MatrixMixtureProfileProcess::Create(innsite,indim);
		SBDPProfileProcess::Create(innsite,indim);
	}

	virtual void Delete()	{
		SBDPProfileProcess::Delete();
		MatrixMixtureProfileProcess::Delete();
	}

	double MixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual double GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual void SlaveMixMove();

	double GlobalIncrementalDPMove(int nrep, double epsilon);
	void SlaveIncrementalDPMove();
	double IncrementalDPMove(int nrep, double epsilon);

};

#endif
