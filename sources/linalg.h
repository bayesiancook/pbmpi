
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef LINALG_H
#define LINALG_H

class LinAlg {

	public:

	// diagonalize a reversible rate matrix
	// first transforms reversible matrix into a symmetric matrix
	// then calls DiagonalizeSymmetricMatrix
	static int DiagonalizeRateMatrix(double** u, double* pi, int dim, double* eigenval, double** eigenvect, double** inveigenvect, int nmax=1000, double epsilon = 1e-10);

	// diagonalize a symmetric matrix
	// first applying Householder transformation (tri-diagonal)
	// then using QR reduction

	// nmax : max number of iterations
	// epsilon : loops until non-diagonal elements are < epsilon in absolute value
	// returns number of iterations
	// eigenvect matrix is orthonormal (its transpose is its inverse)
	static int DiagonalizeSymmetricMatrix(double** u, int dim, int nmax, double epsilon, double* eigenval, double** eigenvect);

	// computes inverse of matrix given as an input (a)
	// by Gauss elimination
	// store inverse in invu
	// does not corrupt matrix a
	// if invu not specified, just returns the logdet
	static double Gauss(double** a, int dim, double** invu = 0);

	private:

	static void QR(double** u, int dim, double** ql, double** r);
	static void HouseHolder(double** u, int dim, double** a, double** ql);

};

#endif
