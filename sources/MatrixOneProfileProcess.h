
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXONEPROFILE_H
#define MATRIXONEPROFILE_H

#include "MatrixProfileProcess.h"
#include "OneProfileProcess.h"

class MatrixOneProfileProcess : public virtual MatrixProfileProcess, public virtual OneProfileProcess {

	public:

	MatrixOneProfileProcess() : matrix(0) {}
	virtual ~MatrixOneProfileProcess() {}

	SubMatrix* GetMatrix(int site)	{
		return matrix;
	}

	protected:

	// called at the beginning and the end of the run
	virtual void Create(int innsite, int indim);
	virtual void Delete();

	double MoveProfile(double tuning = 1, int n = 1, int nrep = 1);

	virtual void UpdateProfileSuffStat() = 0;

	// should be called each time global parameters are modified
	void UpdateProfile()	{
		UpdateMatrix();
	}

	virtual void CreateMatrix() = 0;

	virtual void DeleteMatrix()	{
		delete matrix;
		matrix = 0;
	}

	virtual void UpdateMatrix() = 0;

	SubMatrix* matrix;
};

#endif
