
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixSubstitutionProcess.h"
#include "Random.h"

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;


void MatrixSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{

	const int nstate = GetMatrix(sitemin)->GetNstate();
	double cumul[nstate];
	double expdiag[nstate];
	for(int i=sitemin; i<sitemax; i++)	{

		int up = stateup[i];

		SubMatrix* matrix = GetMatrix(i);
		double** eigenvect = matrix->GetEigenVect();
		double** inveigenvect = matrix->GetInvEigenVect();
		double* eigenval = matrix->GetEigenVal();

		int j = ratealloc[i];
		double length = time * GetRate(i,j);
		for (int k=0; k<nstate; k++)	{
			expdiag[k] = exp(length * eigenval[k]);
		}

		double totprob = 0;
		for (int k=0; k<nstate; k++)	{
			double tot = 0;
			for (int l=0; l<nstate; l++)	{
				tot += eigenvect[up][l] * expdiag[l] * inveigenvect[l][k];
			}
			totprob += tot;
			cumul[k] = totprob;
		}
		if (fabs(totprob - 1) > 1e-6)	{
			cerr << "error in MatrixSubstitutionProcess::SimuPropagate: tot prob is not 1\n";
			cerr << totprob << '\n';
			exit(1);
		}

		double u = rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<nstate) && (u>cumul[k]))	{
			k++;
		}
		if (k == nstate)	{
			cerr << "error in MatrixSubstitutionProcess::SimuPropagate: overflow\n";
			exit(1);
		}
		statedown[i] = k;
	}
}

//-------------------------------------------------------------------------
//	(CPU level 3)
//
//	* conditional likelihood propagation
//
//	(CPU level 3)
//-------------------------------------------------------------------------

void MatrixSubstitutionProcess::Propagate(double*** from, double*** to, double time, bool condalloc)	{
	PropagateImpl(from,to,time,condalloc,0);
}

void MatrixSubstitutionProcess::PropagateFromLeaf(const int* state, double*** from, double*** to, double time, bool condalloc)	{
	PropagateImpl(from,to,time,condalloc,state);
}

// substitution matrix Q = P L P^{-1} where L is diagonal (eigenvalues) and P is the eigenvector matrix
// we need to compute
// down = exp(length * Q) . up
// which we express as
// down = P ( exp(length * L) . (P^{-1} . up) )
//
// thus we successively do the following matrix.vector products
//
// P^{-1} . up  -> aux
// exp(length * L) . aux  -> aux 	(where exp(length*L) is diagonal, so this is linear)
// P . aux -> down
//
// Two things are done to save time, without changing the result:
// - exp(length * L) only depends on the matrix and on the rate: it is computed once
//   per (matrix,rate) pair and then shared by all sites having that matrix (mixture models)
// - if leafstate is given, 'up' is known to be the unit vector of the observed state
//   (for non-missing data), in which case P^{-1} . up is just a column of P^{-1}

void MatrixSubstitutionProcess::PropagateImpl(double*** from, double*** to, double time, bool condalloc, const int* leafstate)	{

	int i,j,k,l;
	const int nstate = GetMatrix(sitemin)->GetNstate();
	if ((int) propaux.size() < nstate)	{
		propaux.resize(nstate);
	}
	double* __restrict__ aux = &propaux[0];
	expindex.clear();
	expbuffer.clear();

	for(i=sitemin; i<sitemax; i++)	{
		if (ActiveSite(i))  {
			SubMatrix* matrix = GetMatrix(i);
			const double* __restrict__ eigenvect_t = matrix->GetTransposedEigenVect();
			const double* __restrict__ inveigenvect_t = matrix->GetTransposedInvEigenVect();
			double* eigenval = matrix->GetEigenVal();
			const int nrate = GetNrate(i);
			const int obs = leafstate ? leafstate[i] : -1;
			for(j=0; j<nrate; j++)	{
				if ((!condalloc) || (ratealloc[i] == j))	{
					const double* up = from[i][j];
					double* down = to[i][j];
					const double rate = GetRate(i,j);

					// exp(length * L)
					ExpKey key;
					key.matrix = matrix;
					key.rate = rate;
					std::pair<std::unordered_map<ExpKey,size_t,ExpKeyHash>::iterator,bool> ins = expindex.insert(std::make_pair(key,expbuffer.size()));
					if (ins.second)	{
						const double length = time * rate;
						expbuffer.resize(expbuffer.size() + nstate);
						double* e = &expbuffer[ins.first->second];
						for(k=0; k<nstate; k++)	{
							e[k] = exp(length * eigenval[k]);
						}
					}
					const double* expdiag = &expbuffer[ins.first->second];

					// The matrix x vector products are done column by column:
					// for each k, the terms are still added in the order l = 0..nstate-1 (same result
					// as a row by row computation), but the inner loop over k can be vectorized.

					// P^{-1} . up  -> aux
					if (obs >= 0)	{
						const double* col = inveigenvect_t + obs*nstate;
						for(k=0; k<nstate; k++)	{
							// 0.0 + x, so as to get exactly what the general case gives (including when x == -0.0)
							aux[k] = 0.0 + col[k];
						}
					}
					else	{
						for(k=0; k<nstate; k++)	{
							if (up[k] < 0.0)	{
								cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
								exit(1);
							}
						}
						for(k=0; k<nstate; k++)	{
							aux[k] = 0.0;
						}
						for(l=0; l<nstate; l++)	{
							const double x = up[l];
							const double* col = inveigenvect_t + l*nstate;
							for(k=0; k<nstate; k++)	{
								aux[k] += col[k] * x;
							}
						}
					}

					// exp(length * L) . aux  -> aux
					for(k=0; k<nstate; k++)	{
						aux[k] *= expdiag[k];
					}

					// P . aux -> down
					for(k=0; k<nstate; k++)	{
						down[k] = 0.0;
					}
					for(l=0; l<nstate; l++)	{
						const double x = aux[l];
						const double* col = eigenvect_t + l*nstate;
						for(k=0; k<nstate; k++)	{
							down[k] += col[k] * x;
						}
					}
					bool isnan = false;
					for(k=0; k<nstate; k++)	{
						if (std::isnan(down[k]))	{
							isnan = true;
						}
					}

					// exit in case of numerical errors
					if (isnan)	{
						cerr << "error in back prop\n";
						for(l=0; l<nstate; l++)	{
							cerr << up[l] << '\t' << down[l] << '\t' << matrix->Stationary(l) << '\n';
						}
						exit(1);
					}
					for(k=0; k<nstate; k++)	{
						if (down[k] < 0.0)	{
							infprobcount++;
							down[k] = 0.0;
						}
					}

					// this is the offset (in log)
					down[nstate] = up[nstate];
				}
			}
		}
	}
}

void MatrixSubstitutionProcess::SitePropagate(int i, double** from, double** to, double time, bool condalloc)	{

	// propchrono.Start();
	int j,k,l;
	double length,max,maxup;

	// should be dependent on site
	// const int nstate = GetMatrix(GetSiteMin())->GetNstate();

	double* aux = new double[GetNstate(i)];

	SubMatrix* matrix = GetMatrix(i);
	double** eigenvect = matrix->GetEigenVect();
	double** inveigenvect = matrix->GetInvEigenVect();
	double* eigenval = matrix->GetEigenVal();

	int nstate = matrix->GetNstate();

	for(j=0; j<GetNrate(i); j++)	{

		if ((!condalloc) || (ratealloc[i] == j))	{

			double* up = from[j];
			double* down = to[j];

			length = time * GetRate(i,j);

			// substitution matrix Q = P L P^{-1} where L is diagonal (eigenvalues) and P is the eigenvector matrix
			// we need to compute 
			// down = exp(length * Q) . up
			// which we express as 
			// down = P ( exp(length * L) . (P^{-1} . up) )  
	
			// thus we successively do the following matrix.vector products

			// P^{-1} . up  -> aux
			// exp(length * L) . aux  -> aux 	(where exp(length*L) is diagonal, so this is linear)
			// P . aux -> down

			// P^{-1} . up  -> aux

			for(k=0; k<nstate; k++)	{
				aux[k] = 0.0;
			}

			for(k=0; k<nstate; k++)	{
				for(l=0; l<nstate; l++)	{
					aux[k] += inveigenvect[k][l] * up[l];
				}
			}

			// exp(length * L) . aux  -> aux
			for(k=0; k<nstate; k++)	{
				aux[k] *= exp(length * eigenval[k]);
			}

			// P . aux -> down
			for(k=0; k<nstate; k++)	{
				down[k] = 0.0;
			}

			for(k=0; k<nstate; k++)	{
				for(l=0; l<nstate; l++)	{
					down[k] += eigenvect[k][l] * aux[l]; 
				}
			}

			// exit in case of numerical errors
			for(k=0; k<nstate; k++)	{
				if (isnan(down[k]))	{
					cerr << "error in back prop\n";
					for(l=0; l<nstate; l++)	{
						cerr << up[l] << '\t' << down[l] << '\t' << matrix->Stationary(l) << '\n';
					}
					exit(1);
				}
			}
			maxup = 0.0;
			for(k=0; k<nstate; k++)	{
				if (up[k] < 0.0)	{
					cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
					exit(1);
				}
				if (maxup < up[k])	{
					maxup = up[k];
				}
			}
			max = 0.0;
			for(k=0; k<nstate; k++)	{
				if (down[k] < 0.0)	{
					infprobcount++;
					down[k] = 0.0;
				}
				if (max < down[k])	{
					max = down[k];
				}
			}
			/*
			if (maxup == 0.0)	{
				cerr << "error in backward propagate: null up array\n";
				cerr << "site : " << i << '\n';
				for(l=0; l<nstate; l++)	{
					cerr << matrix->Stationary(l) << '\n';
				}
				cerr << time << '\t' << length << '\n';
				exit(1);
			}
			if (max == 0.0)	{
				cerr << "error in backward propagate: null array\n";
				for(k=0; k<nstate; k++)	{
					cerr << up[k] << '\t' << down[k] << '\n';
				}
				cerr << length << '\n';
				cerr << '\n';
				exit(1);
			}
			*/
			// this is the offset (in log)
			down[nstate] = up[nstate];
		}
	}

	delete[] aux;
}

