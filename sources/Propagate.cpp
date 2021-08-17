
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

	// propchrono.Start();
	int i,j,k,l,offset;
	double length,max,maxup;
	const int nstate = GetMatrix(sitemin)->GetNstate();
	// double* bigaux = new double[(sitemax - sitemin) * GetNrate(0) * nstate];
	double* aux = new double[GetNsite() * GetNrate(0) * nstate];
	for(i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            SubMatrix* matrix = GetMatrix(i);
            double** eigenvect = matrix->GetEigenVect();
            double** inveigenvect = matrix->GetInvEigenVect();
            double* eigenval = matrix->GetEigenVal();
            for(j=0; j<GetNrate(i); j++)	{
                if ((!condalloc) || (ratealloc[i] == j))	{
                    double* up = from[i][j];
                    double* down = to[i][j];
                    //SubMatrix* matrix = GetMatrix(i);
                    length = time * GetRate(i,j);

                    //double** eigenvect = matrix->GetEigenVect();
                    //double** inveigenvect= matrix->GetInvEigenVect();
                    //double* eigenval = matrix->GetEigenVal();

                    // substitution matrix Q = P L P^{-1} where L is diagonal (eigenvalues) and P is the eigenvector matrix
                    // we need to compute 
                    // down = exp(length * Q) . up
                    // which we express as 
                    // down = P ( exp(length * L) . (P^{-1} . up) )  
            
                    // thus we successively do the following matrix.vector products

                    // P^{-1} . up  -> aux
                    // exp(length * L) . aux  -> aux 	(where exp(length*L) is diagonal, so this is linear)
                    // P . aux -> down

                    /*
                    int nstate = GetNstate();
                    double* aux = new double[nstate];
                    */
                    //double* aux = bigaux + nstate * (i*GetNrate(0)  + j);
                    offset = nstate*(i*GetNrate(0) + j);
                    // P^{-1} . up  -> aux
                    //double* tmpaux = aux;
                    for(k=0; k<nstate; k++)	{
                        //(*tmpaux++) = 0;
                        aux[offset+k] = 0.0;
                    }
                    //tmpaux -= nstate;
                    //double* tmpup = up;
                    for(k=0; k<nstate; k++)	{
                        //double* tmpinveigen = inveigenvect[i];
                        for(l=0; l<nstate; l++)	{
                            //(*tmpaux) += (*tmpinveigen++) * (*tmpup++);
                            aux[offset+k] += inveigenvect[k][l] * up[l];
                        }
                        //tmpaux++;
                        //tmpup -= nstate;
                    }
                    //tmpaux -= nstate;

                    // exp(length * L) . aux  -> aux
                    //double* tmpval = eigenval;
                    for(k=0; k<nstate; k++)	{
                        //(*tmpaux++) *= exp(length * (*tmpval++));
                        aux[offset+k] *= exp(length * eigenval[k]);
                    }
                    //tmpaux -= nstate;
                    //tmpval -= nstate;

                    // P . aux -> down
                    //double* tmpdown = down;
                    for(k=0; k<nstate; k++)	{
                        //(*tmpdown++) = 0;
                        down[k] = 0.0;
                    }
                    //tmpdown -= nstate;

                    for(k=0; k<nstate; k++)	{
                        //double* tmpeigen = eigenvect[i];
                        for(l=0; l<nstate; l++)	{
                            //(*tmpdown) += (*tmpeigen++) * (*tmpaux++);
                            down[k] += eigenvect[k][l] * aux[offset+l]; 
                        }
                        //tmpdown++;
                        //tmpaux -= nstate;
                    }
                    //tmpdown -= nstate;

                    // exit in case of numerical errors
                    for(k=0; k<nstate; k++)	{
                        if (std::isnan(down[k]))	{
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
                        cerr << GetDim() << '\n';
                        cerr << GetMinStat(i) << '\n';
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
        }
    }
    delete[] aux;
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

