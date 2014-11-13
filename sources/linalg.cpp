
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

#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

void LinAlg::QR(double** u, int dim, double** ql, double** r)	{

	double* v = new double[dim];
	double* c = new double[dim];

	double** a = r;
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			a[i][j] = u[i][j];
		}
	}

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			ql[i][j] = 0;
		}
	}
	for (int i=0; i<dim; i++)	{
		ql[i][i] = 1;
	}

	for (int k=0; k<dim-1; k++)	{

		double alpha2 = 0;
		for (int j=k; j<dim; j++)	{
			alpha2 += a[j][k] * a[j][k];
		}
		 
		double alpha = sqrt(alpha2);
		if (a[k][k] > 0)	{
			alpha = -alpha;
		}

		for (int j=0; j<k; j++)	{
			v[j] = 0;
		}

		v[k] = (a[k][k] - alpha);
		for (int j=k+1; j<dim; j++)	{
			v[j] = a[j][k];
		}
		double r = 0;
		for (int j=k; j<dim; j++)	{
			r += v[j] * v[j];
		}
		r = sqrt(r);
		for (int j=k; j<dim; j++)	{
			v[j] /= r;
		}
		
		for (int i=0; i<dim; i++)	{
			double tot = 0;
			for (int j=0; j<dim; j++)	{
				tot += a[j][i] * v[j];
			}
			c[i] = tot;
		}

		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				a[i][j] -= 2 * v[i] * c[j];
			}
		}

		for (int i=0; i<dim; i++)	{
			double tot = 0;
			for (int j=0; j<dim; j++)	{
				tot += ql[j][i] * v[j];
			}
			c[i] = tot;
		}

		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				ql[i][j] -= 2 * v[i] * c[j];
			}
		}

	}

	delete[] v;
	delete[] c;
}


void LinAlg::HouseHolder(double** u, int dim, double** a, double** ql)	{

	double* v = new double[dim];
	double* c = new double[dim];

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			a[i][j] = u[i][j];
		}
	}

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			ql[i][j] = 0;
		}
	}
	for (int i=0; i<dim; i++)	{
		ql[i][i] = 1;
	}

	for (int k=0; k<dim-2; k++)	{

		double alpha2 = 0;
		for (int j=k+1; j<dim; j++)	{
			alpha2 += a[j][k] * a[j][k];
		}
		 
		double alpha = sqrt(alpha2);
		if (a[k+1][k] > 0)	{
			alpha = -alpha;
		}

		double r = sqrt(0.5 * (alpha2 - alpha * a[k+1][k]));

		for (int j=0; j<=k; j++)	{
			v[j] = 0;
		}

		v[k+1] = (a[k+1][k] - alpha) / 2 / r;
		for (int j=k+2; j<dim; j++)	{
			v[j] = a[j][k] / 2 / r;
		}

		for (int i=0; i<dim; i++)	{
			double tot = 0;
			for (int j=k+1; j<dim; j++)	{
				tot += a[i][j] * v[j];
			}
			c[i] = tot;
		}

		double tot = 0;
		for (int j=k+1; j<dim; j++)	{
			tot += v[j] * c[j];
		}

		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				a[i][j] += - 2 * c[i] * v[j] - 2 * c[j] * v[i] + 4 * tot * v[i] * v[j];
			}
		}

		for (int i=0; i<dim; i++)	{
			double tot = 0;
			for (int j=k+1; j<dim; j++)	{
				tot += ql[j][i] * v[j];
			}
			c[i] = tot;
		}

		for (int i=k+1; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				ql[i][j] -= 2 * v[i] * c[j];
			}
		}
	}

	delete[] v;
	delete[] c;
}


int LinAlg::DiagonalizeSymmetricMatrix(double** u, int dim, int nmax, double epsilon, double* eigenval, double** eigenvect)	{

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			eigenvect[i][j] = 0;
		}
	}
	for (int i=0; i<dim; i++)	{
		eigenvect[i][i] = 1;
	}

	double premax = 0;
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			if (i!=j)	{
				double tmp = fabs(u[i][j]);
				if (premax < tmp)	{
					premax = tmp;
				}
			}
		}
	}
	if (premax < epsilon)	{
		for (int i=0; i<dim; i++)	{
			eigenval[i] = u[i][i];
		}
		return 0;
	}

	double** a = new double*[dim];
	double** q = new double*[dim];
	double** r = new double*[dim];
	for (int i=0; i<dim; i++)	{
		a[i] = new double[dim];
		q[i] = new double[dim];
		r[i] = new double[dim];
		for (int j=0; j<dim; j++)	{
			a[i][j] = u[i][j];
		}
	}

	HouseHolder(u,dim,a,r);

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			eigenvect[i][j]  = r[j][i];
		}
	}

	int n = 0;
	double max = 0;
	int s = dim-1;
	do	{
		max = 0;
		double shift = 0;
		if (s>0)	{
			shift = a[s][s];
			for (int i=0; i<s+1; i++)	{
				a[i][i] -= shift;
			}
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				q[i][j] = 0;
				r[i][j] = 0;
			}
		}
		for (int i=0; i<dim; i++)	{
			q[i][i] = 1;
		}

		QR(a,s+1,q,r);
		// QR(a,dim,q,r);
		n++;
		for (int i=0; i<s+1; i++)	{
			for (int j=0; j<s+1; j++)	{
				double total = 0;
				for (int k=0; k<s+1; k++)	{
					total += r[i][k] * q[j][k];
				}
				a[i][j] = total;
				if (i!=j)	{
					if (max < fabs(total))	{
						max = fabs(total);
					}
				}
			}
		}
		if (s>0)	{
			for (int i=0; i<s+1; i++)	{
				a[i][i] += shift;
			}
			if (fabs(a[s][s-1]) < epsilon)	{
				s--;
			}
		}
		else	{
			s--;
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				double total = 0;
				for (int k=0; k<dim; k++)	{
					total += eigenvect[i][k] * q[j][k];
				}
				r[i][j] = total;
			}
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				eigenvect[i][j] = r[i][j];
			}
		}
	}
	while ((s >= 0) && (n < nmax));
	// while ((max > epsilon) && (n < nmax));

	for (int i=0; i<dim; i++)	{
		eigenval[i] = a[i][i];
	}

	for (int i=0; i<dim; i++)	{
		delete[] a[i];
		delete[] q[i];
		delete[] r[i];
	}
	delete[] a;
	delete[] q;
	delete[] r;
	return n;
}

/*
int LinAlg::DiagonalizeSymmetricMatrix(double** u, int dim, int nmax, double epsilon, double* eigenval, double** eigenvect)	{

	double** a = new double*[dim];
	double** q = new double*[dim];
	double** r = new double*[dim];
	for (int i=0; i<dim; i++)	{
		a[i] = new double[dim];
		q[i] = new double[dim];
		r[i] = new double[dim];
		for (int j=0; j<dim; j++)	{
			a[i][j] = u[i][j];
			eigenvect[i][j] = 0;
		}
	}
	for (int i=0; i<dim; i++)	{
		eigenvect[i][i] = 1;
	}

	double premax = 0;
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			if (i!=j)	{
				double tmp = fabs(u[i][j]);
				if (premax < tmp)	{
					premax = tmp;
				}
			}
		}
	}
	if (premax < epsilon)	{
		for (int i=0; i<dim; i++)	{
			eigenval[i] = u[i][i];
		}
		// here delete !!!
		for (int i=0; i<dim; i++)	{
			delete[] a[i];
			delete[] q[i];
			delete[] r[i];
		}
		delete[] a;
		delete[] q;
		delete[] r;
		return 0;
	}

	HouseHolder(u,dim,a,r);

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			eigenvect[i][j]  = r[j][i];
		}
	}

	int n = 0;
	double max = 0;
	int s = dim-1;
	do	{
		max = 0;
		double shift = 0;
		if (s>=0)	{
			shift = a[s][s];
			for (int i=0; i<s+1; i++)	{
				a[i][i] -= shift;
			}
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				q[i][j] = 0;
				r[i][j] = 0;
			}
		}
		for (int i=0; i<dim; i++)	{
			q[i][i] = 1;
		}

		QR(a,s+1,q,r);
		// QR(a,dim,q,r);
		n++;
		for (int i=0; i<s+1; i++)	{
			for (int j=0; j<s+1; j++)	{
				double total = 0;
				for (int k=0; k<s+1; k++)	{
					total += r[i][k] * q[j][k];
				}
				a[i][j] = total;
				if (i!=j)	{
					if (max < fabs(total))	{
						max = fabs(total);
					}
				}
			}
		}
		if (s>=0)	{
			for (int i=0; i<s+1; i++)	{
				a[i][i] += shift;
			}
			// cerr << a[s][s-1] << '\n';
			if (fabs(a[s][s-1]) < epsilon)	{
				s--;
			}
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				double total = 0;
				for (int k=0; k<dim; k++)	{
					total += eigenvect[i][k] * q[j][k];
				}
				r[i][j] = total;
			}
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				eigenvect[i][j] = r[i][j];
			}
		}
	}
	while ((s >=0) && (n < nmax));
	// while ((max > epsilon) && (n < nmax));

	for (int i=0; i<dim; i++)	{
		eigenval[i] = a[i][i];
	}

	for (int i=0; i<dim; i++)	{
		delete[] a[i];
		delete[] q[i];
		delete[] r[i];
	}
	delete[] a;
	delete[] q;
	delete[] r;
	return n;
};
*/

// diagonalize a reversible rate matrix
// first transforms reversible matrix into a symmetric matrix
// then use Householder's algorithm, com
int LinAlg::DiagonalizeRateMatrix(double** u, double* pi, int dim, double* eigenval, double** eigenvect, double** inveigenvect, int nmax, double epsilon)	{

	double** a = new double*[dim];
	for (int i=0; i<dim; i++)	{
		a[i] = new double[dim];
		for (int j=0; j<dim; j++)	{
			a[i][j] = u[i][j] * sqrt(pi[i] / pi[j]);
		}
	}

	int n = DiagonalizeSymmetricMatrix(a,dim,nmax,epsilon,eigenval,eigenvect);

	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			double total = 0;
			for (int k=0; k<dim; k++)	{
				total += a[i][k] * eigenvect[k][j];
			}
		}
	}
	
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			inveigenvect[i][j] = eigenvect[j][i] * sqrt(pi[j]);
		}
	}
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			eigenvect[i][j] /= sqrt(pi[i]);
		}
	}

	for (int i=0; i<dim; i++)	{
		delete[] a[i];
	}
	delete[] a;

	return n;
}

// computes inverse of matrix given as an input (a)
// store inverse in invu
// does not corrupt matrix a
double LinAlg::Gauss(double** a, int dim, double** invu)	{

	// returns log |determinant|
	bool deleteinvu = false;
	if (! invu)	{
		deleteinvu = true;
		invu = new double*[dim];
		for (int i=0; i<dim; i++)	{
			invu[i] = new double[dim];
		}
	}
	double** u = new double*[dim];
	for (int i=0; i<dim; i++)	{
		u[i] = new double[dim];
		for (int j=0; j<dim; j++)	{
			u[i][j] = a[i][j];
		}
	}

	double logdet = 0;
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			invu[i][j] = 0;
		}
	}
	for (int i=0; i<dim; i++)	{
		invu[i][i] = 1;
	}

	for (int i=0; i<dim; i++)	{
		// find max
		double max = 0;
		int imax = 0;
		for (int j=i; j<dim; j++)	{
			double tmp = fabs(u[j][i]);
			if (max < tmp)	{
				max = tmp;
				imax = j;
			}
		}

		if (max < 1e-10)	{
			cerr << "error in Gauss: non invertible matrix\n";
			exit(1);
		}

		// swap
		if (imax != i)	{
			for (int j=0; j<dim; j++)	{
				double tmp = u[i][j];
				u[i][j] = u[imax][j];
				u[imax][j] = tmp;
				tmp = invu[i][j];
				invu[i][j] = invu[imax][j];
				invu[imax][j] = tmp;
			}
		}

		// linear combination 
		for (int j=i+1; j<dim; j++)	{
			double f = u[j][i] / u[i][i];
			for (int k=0; k<dim; k++)	{
				u[j][k] -= u[i][k] * f;
				invu[j][k] -= invu[i][k] * f;
			}
		} 

		double f = u[i][i];
		for (int k=0; k<dim; k++)	{
			u[i][k] /= f;
			invu[i][k] /= f;
		}
		logdet -= log(fabs(f));
	}

	// elimination
	for (int i=dim-1; i>=0; i--)	{

		for (int j=0; j<i; j++)	{
			double f = u[j][i];
			for (int k=0; k<dim; k++)	{
				invu[j][k] -= f * invu[i][k];
				u[j][k] -= f * u[i][k];
			}
		}
	}

	for (int i=0; i<dim; i++)	{
		delete[] u[i];
	}
	delete[] u;
	if (deleteinvu)	{
		for (int i=0; i<dim; i++)	{
			delete[] invu[i];
		}
		delete[] invu;
	}

	return logdet;
}


