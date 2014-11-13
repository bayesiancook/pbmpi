
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


// SubMatrix:
// this class implements
// the instant rate matrix of a substitution process
// but only the mathematical aspects of it
//
// RandomSubMatrix:
// this class derives from SubMatrix and from Dnode
// therefore, it knows everything about substitution processes (as a SubMatrix)
// and at the same time, it can be inserted as a deterministic node in a probabilistic model (as a Dnode)
//
// If you need to define a new substitution process
// - derive a new class deriving (directly or indirectly) from SubMatrix
// in this class, implements the ComputeArray and ComputeStationary functions
// (these are the functions that construct the instant rates and the stationary probabilities of the matrix)
//
// - derive a new class deriving both from RandomSubMatrix, and from your new class
// in this class, implement SetParameters
// (this function is responsible for updating the internal parameters that the SubMatrix uses in ComputeArray And ComputeStationary,
// based on the values stored by the parent nodes)
//
//

#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <iostream>
#include <cmath>
using namespace std;

#include "Random.h"

class SubMatrix  	{


	protected:

	// give a range (a time, or a branch length in expected number of substitutions)
	// give a temporary Nstate*Nstate array (temp)
	// give another Nstat*Nstate array (expo)
	// return the exponential of the matrix in expo
	void			ComputeExponential(double range, double** expo, double** temp);

	// these 2 pure virtual functions are the most essential component of the SubMatrix class
	// see GTRSubMatrix.cpp and CodonSubMatrix.cpp for examples

	// ComputeArray(int state) is in charge of computing the row of the rate matrix
	// corresponding to all possible rates of substitution AWAY from state
	//
	virtual void 		ComputeArray(int state) = 0;

	// ComputeStationary() is in charge of computing the vector of stationary probabilities (equilibirum frequencies)
	// of the substitution process
	virtual void 		ComputeStationary() = 0;

	public:

	static const int	UniSubNmax = 500;

	static int		nuni;
	static int		nunimax;

	static int		nunisubcount;

	static int		GetUniSubCount() {return nunisubcount;}

	static double		GetMeanUni() {return ((double) nunimax) / nuni;}

				SubMatrix(int Nstate, bool innormalise = false);
	virtual 		~SubMatrix();

	void			Create();

	double 			operator()(int, int);
	const double* 		GetRow(int i);

	virtual const double* 	GetStationary();
	double 			Stationary(int i);

	int 			GetNstate() {return Nstate;}

	virtual double		GetRate();
	void 			ScalarMul(double e);

	bool			isNormalised() {return normalise;}
	void			Normalise();

	void 			CorruptMatrix();
	void 			UpdateMatrix();

	void			ActivatePowers();
	void			InactivatePowers();
	double 			Power(int n, int i, int j);
	double			GetUniformizationMu();

	double* 		GetEigenVal();
	double** 		GetEigenVect();
	double** 		GetInvEigenVect();


	// uniformization resampling methods
	// CPU level 1
	double 			GetFiniteTimeTransitionProb(int stateup, int statedown, double efflength);
	int 			DrawUniformizedTransition(int state, int statedown, int n);
	int 			DrawUniformizedSubstitutionNumber(int stateup, int statedown, double efflength);
	//

	// used by accept-reject resampling method
	// CPU level 1
	int 			DrawOneStep(int state);

	void			ToStream(ostream& os);
	void			CheckReversibility();

	int 			GetDiagStat() { return ndiagfailed;}


	protected:

	void 			UpdateRow(int state);
	void 			UpdateStationary();


	void 			ComputePowers(int n);
	void 			CreatePowers(int n);

	bool			ArrayUpdated();

	int 			Diagonalise();

	// data members
	
	bool powflag;
	bool diagflag;
	bool statflag;
	bool* flagarray;

	int Nstate;
	int npow;
	double UniMu;

	double*** mPow;
	
	// Q : the infinitesimal generator matrix
	double ** Q;

	// the stationary probabilities of the matrix
	double* mStationary;

	bool normalise;

	private:

	// v : eigenvalues
	// vi : imaginary part
	// u : the matrix of eigen vectors
	// invu : the inverse of u


	double ** u;
	double ** invu;
	double * v;
	double * vi;

	int ndiagfailed;
};


//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------

inline double	SubMatrix::operator()(int i, int j)	{
	if (! flagarray[i])	{
		UpdateRow(i);
	}
	return Q[i][j];
}

inline const double* SubMatrix::GetRow(int i) {
	if (! flagarray[i])	{
		UpdateRow(i);
	}
	return Q[i];
}

inline const double* SubMatrix::GetStationary() {
	if (! statflag)	{
		UpdateStationary();
	}
	return mStationary;
}

inline double SubMatrix::Stationary(int i)	{
	if (! statflag)	{
		UpdateStationary();
	}
	return mStationary[i];
}


inline void SubMatrix::CorruptMatrix()	{
	diagflag = false;
	statflag = false;
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = false;
	}
	InactivatePowers();
}

inline bool SubMatrix::ArrayUpdated()	{

	bool qflag = true;
	for (int k=0; k<Nstate; k++)	{
		qflag &= flagarray[k];
	}
	return qflag;
}

inline void SubMatrix::UpdateStationary()	{
	ComputeStationary();
	statflag = true;
}
	
inline void SubMatrix::UpdateRow(int state)	{

	if (isNormalised())	{
		UpdateMatrix();
	}
	else	{
		if (! statflag)	{
			UpdateStationary();
		}
		ComputeArray(state);
		flagarray[state] = true;
	}
}

inline double SubMatrix::GetFiniteTimeTransitionProb(int stateup, int statedown, double efflength)	{
	double** invp = GetInvEigenVect();
	double** p = GetEigenVect();
	double* l = GetEigenVal();
	double tot = 0;
	for (int i=0; i<GetNstate(); i++)	{
		tot += p[stateup][i] * exp(efflength * l[i]) * invp[i][statedown];
	}
	return tot;
}

inline int SubMatrix::DrawUniformizedTransition(int state, int statedown, int n)	{

	double* p = new double[GetNstate()];
	double tot = 0;
	for (int l=0; l<GetNstate(); l++)	{
		tot += Power(1,state,l) * Power(n,l,statedown);
		p[l] = tot;
	}

	double s = tot * rnd::GetRandom().Uniform();
	int k = 0;
	while ((k<GetNstate()) && (s > p[k]))	{
		k++;
	}
	if (k == GetNstate())	{
		cerr << "error in DrawUniformizedTransition: overflow\n";
		throw;
	}
	delete[] p;
	return k;
}


inline int SubMatrix::DrawUniformizedSubstitutionNumber(int stateup, int statedown, double efflength)	{

	double mu = GetUniformizationMu();
	double fact = exp(- efflength * mu);
	int m = 0;
	double total = (stateup==statedown) * fact;
	double Z = GetFiniteTimeTransitionProb(stateup,statedown,efflength);
	double q = rnd::GetRandom().Uniform() * Z;
	
	while ((m<UniSubNmax) && (total < q)) 	{
		m++;
		fact *= mu * efflength / m;
		total += Power(m,stateup,statedown) * fact;
		if ((total-Z)>1e-12)	{
			cerr << "error in DrawUniformizedSubstitutionNumber: normalising constant\n";
			cerr << total << '\t' << Z << '\n';
			cerr << mu << '\n';
			cerr << m << '\n';
			cerr << stateup << '\t' << statedown << '\n';

			// ToStream(cerr);
			// CheckReversibility();
			throw;
		}
	}
	if (m == UniSubNmax)	{
		nunisubcount++;
	}
	return m;
}

inline int SubMatrix::DrawOneStep(int state)	{

	const double* row = GetRow(state);
	double p = -row[state] * rnd::GetRandom().Uniform();
	int k = -1;
	double tot = 0;
	do	{
		k++;
		if (k != state)	{
			tot += row[k];
		}
	}
	while ((k<GetNstate()) && (tot < p));
	if (k == GetNstate())	{
		cerr << "error in DrawOneStep\n";
		cerr << GetNstate() << '\n';
		for (int k=0; k<GetNstate(); k++)	{
			cerr << row[k] << '\n';
		}
		// ToStream(cerr);	
		exit(1);
	}
	return k;
}

#endif // SUBMATRIX_H
