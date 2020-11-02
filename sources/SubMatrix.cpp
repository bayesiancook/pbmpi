
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "linalg.h"
#include "SubMatrix.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

int SubMatrix::nuni = 0;
int SubMatrix::nunimax = 0;
int SubMatrix::nunisubcount = 0;

// ---------------------------------------------------------------------------
//		 SubMatrix()
// ---------------------------------------------------------------------------
/*
SubMatrix::SubMatrix() : Nstate(0), normalise(false) {}


SubMatrix::SubMatrix(const SubMatrix& from) : Nstate(from.inNstate), normalise(from.innormalise)	{
	Create();
}

SubMatrix& SubMatrix::operator=(const SubMatrix& from)	{

	if (Nstate != from.Nstate)	{
		cerr << "error in SubMatrix::operator= : non matching dimensions : " << Nstate << " and " << from.Nstate << '\n';
		throw;
	}
	// copy?
	// Corrupt ?
	CorruptMatrix();
}
*/
SubMatrix::SubMatrix(int inNstate, bool innormalise) : Nstate(inNstate), normalise(innormalise)	{
	ndiagfailed = 0;
	Create();
}

void SubMatrix::Create()	{

	Q = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		Q[i] = new double[Nstate];
	}

	u = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		u[i] = new double[Nstate];
	}

	invu = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		invu[i] = new double[Nstate];
	}

	v = new double[Nstate];
	vi = new double[Nstate];

	mStationary = new double[Nstate];

	UniMu = 1;
	mPow = new double**[UniSubNmax];
	for (int n=0; n<UniSubNmax; n++)	{
		mPow[n] = 0;
		/*
		mPow[n] = new double*[Nstate];
		for (int i=0; i<Nstate; i++)	{
			mPow[n][i] = new double[Nstate];
		}
		*/
	}

	flagarray = new bool[Nstate];
	diagflag = false;
	statflag = false;
	for (int i=0; i<Nstate; i++)	{
		flagarray[i] = false;
	}		
	powflag = false;

}

// ---------------------------------------------------------------------------
//		 ~SubMatrix()
// ---------------------------------------------------------------------------


SubMatrix::~SubMatrix()	{

	for (int i=0; i<Nstate; i++)	{
		delete[] Q[i];
		delete[] u[i];
		delete[] invu[i];
	}
	delete[] Q;
	delete[] u;
	delete[] invu;

	if (mPow)	{
		for (int n=0; n<UniSubNmax; n++)	{
			if (mPow[n])	{
				for (int i=0; i<Nstate; i++)	{
					delete[] mPow[n][i];
				}
				delete[] mPow[n];
			}
		}
		delete[] mPow;
	}
	delete[] mStationary;
	delete[] flagarray;
	delete[] v;
	delete[] vi;


}

// ---------------------------------------------------------------------------
//		 void ScalarMul()
// ---------------------------------------------------------------------------

void	SubMatrix::ScalarMul(double e)	{

	for (int i=0; i< Nstate ; i++)	{
		v[i] *= e;
		vi[i] *= e;
	}
	UniMu *= e;
}


// ---------------------------------------------------------------------------
//		 Diagonalise()
// ---------------------------------------------------------------------------

int SubMatrix::Diagonalise()	{

	if (! ArrayUpdated())	{
		UpdateMatrix();
	}

	// CheckQ();

	int nmax = 1000;
	double epsilon = 1e-20;
	double temptoosmall = 1e-20;
	int reducedStateCount = Nstate;
	int stateNull[Nstate];
	int isnull = 0;
	for (int i=0; i<Nstate; i++)	{
		stateNull[i] = 0;
	}
	for (int i=0; i<Nstate; i++)	{
		if (mStationary[i] < temptoosmall)	{
		//if (mStationary[i] == 0)	{
			stateNull[i] = 1;
			reducedStateCount--;
			isnull = 1;
		}
		else	{
			stateNull[i] = 0;
		}
	}
	bool failed;
	int n;
	if (!isnull)	{
		n = LinAlg::DiagonalizeRateMatrix(Q,mStationary,Nstate,v,u,invu,nmax,epsilon);
		failed = (n == nmax);
	}
	else {	

		double** reducedQ = new double*[reducedStateCount];
		double** reducedu = new double*[reducedStateCount];
		double** reducedinvu = new double*[reducedStateCount];
		double* reducedv = new double[reducedStateCount];
		double* reducedPi = new double[reducedStateCount];

		for (int i=0; i<reducedStateCount; i++)	{
			reducedQ[i] = new double[reducedStateCount];
			reducedu[i] = new double[reducedStateCount];
			reducedinvu[i] = new double[reducedStateCount];
		}

		int counti = 0;
		int countj;

		double reducedPiTotal = 0;	
		for (int i=0; i<Nstate; i++)	{
			countj = 0;
			if (!stateNull[i])	{
				reducedPi[counti] = mStationary[i];
				reducedPiTotal += reducedPi[counti];
				for (int j=0; j<Nstate; j++)	{
					if (!stateNull[j])	{
						reducedQ[counti][countj] = Q[i][j];
						countj++;
					}
				}
				if (countj != reducedStateCount)	{
					cerr << "countj not equal to reducedStateCount\n";
					cerr << "countj = " << countj << ", and reducedStateCount = " << reducedStateCount << "\n";
				}
				counti++;
			}
		}

		for (int i=0; i<reducedStateCount; i++)	{
			reducedPi[i] /= reducedPiTotal;
		}

		if (counti != reducedStateCount)	{
			cerr << "counti not equal to reducedStateCount\n";
		}


		n = LinAlg::DiagonalizeRateMatrix(reducedQ,reducedPi,reducedStateCount,reducedv,reducedu,reducedinvu,nmax,epsilon);
		failed = (n == nmax);

		LinAlg::Gauss(reducedu, reducedStateCount, reducedinvu);

		for (int i=0; i<Nstate; i++)	{
			v[i] = 0;
			for (int j=0; j<Nstate; j++)	{
				if (i==j)	{
					u[i][j] = 1;
					invu[i][j] = 1;
				}
				else	{
					u[i][j] = 0;
					invu[i][j] = 0;
				}
			}
		}

		counti = 0;
		for (int i=0; i<Nstate; i++)	{
			countj = 0;
			if (!stateNull[i])	{
				v[i] = reducedv[counti];
				for (int j=0; j<Nstate; j++)	{
					if (!stateNull[j])	{
						u[i][j] = reducedu[counti][countj];
						invu[i][j] = reducedinvu[counti][countj];
						countj++;
					}	
				}
				counti++;
			}
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				if (std::isnan(u[i][j]))	{
					cerr << "nan in diag: u\n";
					exit(1);
				}
				if (std::isnan(invu[i][j]))	{
					cerr << "nan in diag: invu\n";
					exit(1);
				}
			}
		}

		//LinAlg::Gauss(u, Nstate, invu);

		for (int i=0; i<reducedStateCount; i++)	{
			delete[] reducedQ[i];
			delete[] reducedu[i];
			delete[] reducedinvu[i];
		}
		delete[] reducedQ;
		delete[] reducedu;
		delete[] reducedinvu;
		delete[] reducedv;
		delete[] reducedPi;
	}

	if (failed)	{
		cerr << "error in SubMatrix::Diagonalise\n";
		CheckReversibility();
		ToStream(cerr);
		exit(1);
	}

	diagflag = true;

	return failed;
}

void	SubMatrix::ComputeExponential(double range, double** expo, double** temp)	{

	Diagonalise();

	for (int i=0; i<Nstate; i++)	{
		double tot = 0;
		for (int j=0; j<Nstate; j++)	{
			double tmp = 0;
			for (int k=0; k<Nstate; k++)	{
				tmp += u[i][k] * invu[k][j] * exp(v[k] * range);
			}
			expo[i][j] = tmp;
			tot += tmp;
			/*
			if (fabs(tmp - expo[i][j])>1e-8)	{
				cerr << "error in compute exponential matrix: " << tmp << '\t' << expo[i][j] << '\n';
				exit(1);
			}
			*/
		}
		if (fabs(tot-1) > 1e-8)	{
			cerr << "error : row does not sum to 1\n";
			exit(1);
		}
	}
}

// ---------------------------------------------------------------------------
//		 ComputeRate()
// ---------------------------------------------------------------------------

double SubMatrix::GetRate()	{

	if (! ArrayUpdated())	{
		UpdateStationary();
		for (int k=0; k<Nstate; k++)	{
			ComputeArray(k);
		}
	}
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = true;
	}
	double norm = 0;
	for (int i=0; i<Nstate-1; i++)	{
		for (int j=i+1; j<Nstate; j++)	{
			norm += mStationary[i] * Q[i][j];
		}
	}
	return 2 * norm;
}


double* SubMatrix::GetEigenVal() {
	if (! diagflag)	{
		Diagonalise();
	}
	return v;
}

double** SubMatrix::GetEigenVect() {
	if (! diagflag)	{
		Diagonalise();
	}
	return u;
}

double** SubMatrix::GetInvEigenVect() {
	if (! diagflag) Diagonalise();
	return invu;
}

// ---------------------------------------------------------------------------
//		 Update()
// ---------------------------------------------------------------------------

void SubMatrix::UpdateMatrix()	{
	UpdateStationary();
	for (int k=0; k<Nstate; k++)	{
		ComputeArray(k);
	}
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = true;
	}
	if (isNormalised())	{
		Normalise();
	}
	// CheckReversibility();
}



// ---------------------------------------------------------------------------
//		 Normalise()
// ---------------------------------------------------------------------------

void SubMatrix::Normalise()	{

	double norm = GetRate();
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] /= norm;
		}
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Powers
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void SubMatrix::ActivatePowers()	{

	if (! powflag)	{
		if (! ArrayUpdated())	{
			UpdateMatrix();
		}

		UniMu = 0;
		for (int i=0; i<Nstate; i++)	{
			if (UniMu < fabs(Q[i][i]))	{
				UniMu = fabs(Q[i][i]);
			}
		}

		CreatePowers(0);
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				mPow[0][i][j] = 0;
			}
		}
		for (int i=0; i<Nstate; i++)	{
			mPow[0][i][i] = 1;
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				mPow[0][i][j] += Q[i][j] / UniMu;
				if (mPow[0][i][j] < 0)	{
					cerr << "error in SubMatrix::ComputePowers: negative prob : ";
					cerr << i << '\t' << j << '\t' << mPow[0][i][j] << '\n';
					cerr << "Nstate : " << Nstate << '\n';
					exit(1);
				}
			}
		}
		npow = 1;
		powflag = true;
	}
}

void SubMatrix::InactivatePowers()	{

	if (powflag)	{
		for (int n=0; n<UniSubNmax; n++)	{
			if (mPow[n])	{
				for (int i=0; i<Nstate; i++)	{
					delete[] mPow[n][i];
				}
				delete [] mPow[n];
				mPow[n] = 0;
			}
		}
		nunimax += npow;
		nuni++;

		npow = 0;
		powflag = false;
	}
}

void SubMatrix::CreatePowers(int n)	{

	if (! mPow[n])	{
		mPow[n] = new double*[Nstate];
		for (int i=0; i<Nstate; i++)	{
			mPow[n][i] = new double[Nstate];
		}
	}
}


double SubMatrix::GetUniformizationMu() {

	if (! powflag)	{
		ActivatePowers();
	}
	return UniMu;
}

double SubMatrix::Power(int n, int i, int j)	{

	if (! powflag)	{
		ActivatePowers();
	}
	if (!n)	{
		return (i == j);
	}
	if (n > UniSubNmax)	{
		return Stationary(j);
	}
	if (n > npow)	{
		ComputePowers(n);
	}
	return mPow[n-1][i][j];
}


void SubMatrix::ComputePowers(int N)	{

	if (! powflag)	{
		ActivatePowers();
	}
	if (N>npow)	{
		for (int n=npow; n<N; n++)	{
			CreatePowers(n);
			for (int i=0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					double& t = mPow[n][i][j];
					t = 0;
					for (int k=0; k<Nstate; k++)	{
						t += mPow[n-1][i][k] * mPow[0][k][j];
					}
				}
			}
		}
		npow = N;
	}
}

void SubMatrix::ToStream(ostream& os)	{

	os << GetNstate() << '\n';
	os << "stationaries: \n";
	for (int i=0; i<GetNstate(); i++)	{
		os << Stationary(i) << '\t';
	}
	os << '\n';
	
	os << "rate matrix\n";
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			os << Q[i][j] << '\t';
		}
		os << '\n';
	}

	os << '\n';
	for (int i=0; i<GetNstate(); i++)	{
		os << v[i] << '\t';
	}
	os << '\n';
	
}

void SubMatrix::CheckReversibility()	{

	double max = 0;
	int imax = 0;
	int jmax = 0;
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=i+1; j<GetNstate(); j++)	{
			double tmp = fabs(Stationary(i) * Q[i][j] - Stationary(j) * Q[j][i]);
			if (max < tmp)	{
				max = tmp;
				imax = i;
				jmax = j;
			}
		}
	}
	if (max > 1e-6)	{
		cerr << "max irreversibility: " << max << '\n';
		cerr << imax << '\t' << jmax << '\t' << Stationary(imax) << '\t' << Q[imax][jmax] << '\t' << Stationary(jmax) << '\t' << Q[jmax][imax] << '\n';
		exit(1);
	}	
}

