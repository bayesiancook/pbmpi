
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

// The matrix x vector products are done column by column: y = M . x, where mt holds M transposed
// (column l of M is mt[l*n .. l*n+n-1]). For each k, the terms are added in the order l = 0..n-1,
// starting from 0.0, with a separate multiplication and addition: this is exactly the arithmetic
// of the plain row by row loops, so the result is bit-identical (as long as the compiler is not
// allowed to contract a*b+c into a fused multiply-add: no -mfma / -march=native without -ffp-contract=off).
//
// With GCC / clang, the k are processed in blocks of 20 using explicit vector types, with the
// accumulators of a block held in registers across the whole loop over l (5 registers of 4 doubles
// with AVX, 10 registers of 2 doubles otherwise). Amino acids are one block; codons 3 blocks + 1.

#if defined(__GNUC__) && !defined(PB_NO_VECTOR_KERNEL)

#ifdef __AVX__
typedef double pbvec __attribute__((vector_size(32), aligned(8), may_alias));
#define PBVEC_SET1(x) ((pbvec){(x),(x),(x),(x)})
#define PBVEC_W 4
#else
typedef double pbvec __attribute__((vector_size(16), aligned(8), may_alias));
#define PBVEC_SET1(x) ((pbvec){(x),(x)})
#define PBVEC_W 2
#endif
#define PBVEC_LD(p) (*(const pbvec*)(p))
#define PBVEC_ST(p,v) (*(pbvec*)(p) = (v))

static inline void MatVecBlock20(const double* __restrict__ mt, const double* __restrict__ x, double* __restrict__ y, const int n)	{
	const pbvec zero = PBVEC_SET1(0.0);
#if PBVEC_W == 4
	pbvec a0 = zero, a1 = zero, a2 = zero, a3 = zero, a4 = zero;
	for(int l=0; l<n; l++)	{
		const pbvec xl = PBVEC_SET1(x[l]);
		const double* col = mt + l*n;
		a0 += PBVEC_LD(col) * xl;
		a1 += PBVEC_LD(col+4) * xl;
		a2 += PBVEC_LD(col+8) * xl;
		a3 += PBVEC_LD(col+12) * xl;
		a4 += PBVEC_LD(col+16) * xl;
	}
	PBVEC_ST(y,a0); PBVEC_ST(y+4,a1); PBVEC_ST(y+8,a2); PBVEC_ST(y+12,a3); PBVEC_ST(y+16,a4);
#else
	pbvec a0 = zero, a1 = zero, a2 = zero, a3 = zero, a4 = zero, a5 = zero, a6 = zero, a7 = zero, a8 = zero, a9 = zero;
	for(int l=0; l<n; l++)	{
		const pbvec xl = PBVEC_SET1(x[l]);
		const double* col = mt + l*n;
		a0 += PBVEC_LD(col) * xl;
		a1 += PBVEC_LD(col+2) * xl;
		a2 += PBVEC_LD(col+4) * xl;
		a3 += PBVEC_LD(col+6) * xl;
		a4 += PBVEC_LD(col+8) * xl;
		a5 += PBVEC_LD(col+10) * xl;
		a6 += PBVEC_LD(col+12) * xl;
		a7 += PBVEC_LD(col+14) * xl;
		a8 += PBVEC_LD(col+16) * xl;
		a9 += PBVEC_LD(col+18) * xl;
	}
	PBVEC_ST(y,a0); PBVEC_ST(y+2,a1); PBVEC_ST(y+4,a2); PBVEC_ST(y+6,a3); PBVEC_ST(y+8,a4);
	PBVEC_ST(y+10,a5); PBVEC_ST(y+12,a6); PBVEC_ST(y+14,a7); PBVEC_ST(y+16,a8); PBVEC_ST(y+18,a9);
#endif
}

static inline void MatVecCols(const double* __restrict__ mt, const double* __restrict__ x, double* __restrict__ y, const int n)	{
	int kb = 0;
	for(; kb+20<=n; kb+=20)	{
		MatVecBlock20(mt+kb,x,y+kb,n);
	}
	if (kb < n)	{
		for(int k=kb; k<n; k++)	{
			y[k] = 0.0;
		}
		for(int l=0; l<n; l++)	{
			const double xl = x[l];
			const double* col = mt + l*n;
			for(int k=kb; k<n; k++)	{
				y[k] += col[k] * xl;
			}
		}
	}
}

#else

static inline void MatVecCols(const double* __restrict__ mt, const double* __restrict__ x, double* __restrict__ y, const int n)	{
	for(int k=0; k<n; k++)	{
		y[k] = 0.0;
	}
	for(int l=0; l<n; l++)	{
		const double xl = x[l];
		const double* col = mt + l*n;
		for(int k=0; k<n; k++)	{
			y[k] += col[k] * xl;
		}
	}
}

#endif

void MatrixSubstitutionProcess::PropagateImpl(double*** from, double*** to, double time, bool condalloc, const int* leafstate)	{

	int i,j,k,l;
	const int nstate = GetMatrix(sitemin)->GetNstate();
	if ((int) propaux.size() < nstate)	{
		propaux.resize(nstate);
	}
	double* __restrict__ aux = &propaux[0];
	propstamp++;
	propentries.clear();
	propheads.clear();
	expbuffer.clear();
	leafindex.clear();
	leafbuffer.clear();
	leafclipped.clear();

	for(i=sitemin; i<sitemax; i++)	{
		if (ActiveSite(i))  {
			SubMatrix* matrix = GetMatrix(i);
			const double* __restrict__ eigenvect_t = matrix->GetTransposedEigenVect();
			const double* __restrict__ inveigenvect_t = matrix->GetTransposedInvEigenVect();
			double* eigenval = matrix->GetEigenVal();
			const int nrate = GetNrate(i);
			const int obs = leafstate ? leafstate[i] : -1;
			if (matrix->propstamp != propstamp)	{
				matrix->propstamp = propstamp;
				matrix->propslot = (int) propheads.size();
				propheads.push_back(-1);
			}
			const int slot = matrix->propslot;
			for(j=0; j<nrate; j++)	{
				if ((!condalloc) || (ratealloc[i] == j))	{
					const double* up = from[i][j];
					double* down = to[i][j];
					const double rate = GetRate(i,j);

					// entry of this (matrix,rate) pair, with exp(length * L)
					int e = propheads[slot];
					while ((e != -1) && (propentries[e].rate != rate))	{
						e = propentries[e].next;
					}
					if (e == -1)	{
						PropEntry entry;
						entry.rate = rate;
						entry.expoff = expbuffer.size();
						entry.next = propheads[slot];
						entry.leafbase = -1;
						const double length = time * rate;
						expbuffer.resize(expbuffer.size() + nstate);
						double* ex = &expbuffer[entry.expoff];
						for(k=0; k<nstate; k++)	{
							ex[k] = exp(length * eigenval[k]);
						}
						e = (int) propentries.size();
						propentries.push_back(entry);
						propheads[slot] = e;
					}
					const double* expdiag = &expbuffer[propentries[e].expoff];

					if (obs >= 0)	{
						// leaf with an observed state: 'up' is the unit vector of that state, so the result
						// only depends on (matrix, rate, state): it is computed once and then copied
						if (propentries[e].leafbase == -1)	{
							propentries[e].leafbase = (int) leafindex.size();
							leafindex.resize(leafindex.size() + nstate, -1);
						}
						int& idx = leafindex[propentries[e].leafbase + obs];
						if (idx == -1)	{
							// P^{-1} . up is just a column of P^{-1}
							const double* col = inveigenvect_t + obs*nstate;
							for(k=0; k<nstate; k++)	{
								// 0.0 + x, so as to get exactly what the general case gives (including when x == -0.0)
								aux[k] = (0.0 + col[k]) * expdiag[k];
							}
							MatVecCols(eigenvect_t, aux, down, nstate);
							bool isnan = false;
							for(k=0; k<nstate; k++)	{
								if (std::isnan(down[k]))	{
									isnan = true;
								}
							}
							if (isnan)	{
								cerr << "error in back prop\n";
								for(l=0; l<nstate; l++)	{
									cerr << up[l] << '\t' << down[l] << '\t' << matrix->Stationary(l) << '\n';
								}
								exit(1);
							}
							int nclipped = 0;
							for(k=0; k<nstate; k++)	{
								if (down[k] < 0.0)	{
									nclipped++;
									down[k] = 0.0;
								}
							}
							infprobcount += nclipped;
							const int newidx = (int) leafclipped.size();
							leafclipped.push_back(nclipped);
							leafbuffer.resize(leafbuffer.size() + nstate);
							double* store = &leafbuffer[(size_t) newidx * nstate];
							for(k=0; k<nstate; k++)	{
								store[k] = down[k];
							}
							// (leafindex was not resized since idx was taken)
							idx = newidx;
						}
						else	{
							const double* store = &leafbuffer[(size_t) idx * nstate];
							for(k=0; k<nstate; k++)	{
								down[k] = store[k];
							}
							infprobcount += leafclipped[idx];
						}
					}
					else	{
						bool negative = false;
						for(k=0; k<nstate; k++)	{
							negative |= (up[k] < 0.0);
						}
						if (negative)	{
							for(k=0; k<nstate; k++)	{
								if (up[k] < 0.0)	{
									cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
									exit(1);
								}
							}
						}

						// P^{-1} . up  -> aux
						MatVecCols(inveigenvect_t, up, aux, nstate);

						// exp(length * L) . aux  -> aux
						for(k=0; k<nstate; k++)	{
							aux[k] *= expdiag[k];
						}

						// P . aux -> down
						MatVecCols(eigenvect_t, aux, down, nstate);

						// rare cases: numerical errors (exit) and slightly negative values (set to 0)
						bool bad = false;
						for(k=0; k<nstate; k++)	{
							bad |= !(down[k] >= 0.0);
						}
						if (bad)	{
							bool isnan = false;
							for(k=0; k<nstate; k++)	{
								if (std::isnan(down[k]))	{
									isnan = true;
								}
							}
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

