
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#define MT_LEN       624

#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

// #define SAFE_EXP(x) ((x)<-200.0 ? 0.0 : exp(x))
#define SAFE_EXP(x) exp(x)


using namespace std;

const double gammacoefs[] = {0.9999999999995183,676.5203681218835,-1259.139216722289,771.3234287757674,-176.6150291498386,12.50734324009056,-0.1385710331296526,0.9934937113930748e-05,0.1659470187408462e-06};
static const double Pi = 3.1415926535897932384626;
static const double Logroot2pi =0.918938533204673;
// const double InfProb = -30;


class Random {
 
	public:

	// static const double INFPROB=250;
	
  	Random(int seed = -1);

	void InitRandom(int seed = -1);

	int GetSeed();

	double Uniform();
	double Gamma(double alpha,double beta);
  	double sNormal(void);
  	double sExpo(void);
  	double sGamma(double);
  	double sGammanew(double);

 	int Choose(int);
	int FiniteDiscrete(int n, const double* probarray);
  	void DrawFromUrn(int*, int n, int N);
	int DrawFromDiscreteDistribution(double* p, int n);
	int DrawFromLogDiscreteDistribution(double* ll, int n);

	double logGamma(double a);

	long int GetCount() {return count;}

	private:

	long int count;
	int Seed;
	int mt_index;
	unsigned long mt_buffer[MT_LEN];


  };


class rnd	{

	private:
	static Random* array;
	static int dim;
	
	public:
	/*
	rnd(int indim = 1, int inseed = -1)	{
		cerr << "in random_array::constructor\n";
		init(indim,inseed);
	}
	*/

	static void init(int indim = 1, int seed = -1)	{
		dim = indim;
		array = new Random[dim];
		for (int i=0; i<dim; i++)	{
			array[i].InitRandom(seed);
			// cerr << "seed " << i << " was : " << array[i].GetSeed() << '\n';
		}
	}

	static Random& GetRandom(int i = -1);
};

#endif // RANDOM_H
