
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXPROFILE_H
#define MATRIXPROFILE_H

#include "SubMatrix.h"
#include "ProfileProcess.h"

// superclass for all Matrix implementations
class MatrixProfileProcess : public virtual ProfileProcess	{

	public:

	MatrixProfileProcess() {}
	virtual ~MatrixProfileProcess() {}

	// access to matrices
	virtual SubMatrix* GetMatrix(int site) = 0;

	protected:

	// create/delete all matrices
	// Create called when deactivating sufficient statistics and activating pruning-based computation (Unfold  in PhyloProcess)
	virtual void CreateMatrices() = 0;
	// Delete called when deactivating pruning-based computation and activating sufficient statistics (Collapse in PhyloProcess)
	virtual void DeleteMatrices() = 0;

	// auxiliary: for Metropolis Moves on profiles
	double ProfileProposeMove(double* profile, double tuning, int n, int K=0, double statmin = 0);

	// updates all matrices
	// (should be called, e.g. when performing a Metropolis on relative exchangeabilities or global mutation parameters)
	virtual void UpdateMatrices() = 0;
	// virtual void DiagonaliseMatrices() = 0;
};

#endif

