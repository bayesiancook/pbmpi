
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionedDGamRateProcess.h"
#include "Random.h"
#include "IncompleteGamma.h"

#include <cassert>
#include "Parallel.h"
#include <string.h>

#include <cmath>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionedDGamRateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PartitionedDGamRateProcess::Create(int inncat, PartitionScheme inscheme)	{
	if (! rate)	{
		RateProcess::Create(inscheme.GetNsite());
		PartitionProcess::Create(inscheme);

		alpha = new double[GetNpart()];
		ratemult = new double[GetNpart()];

		rate = new double*[GetNpart()];
		ratesuffstatcount = new int*[GetNpart()];
		ratesuffstatbeta = new double*[GetNpart()];

		Ncat = inncat;

		for(int i = 0; i < GetNpart(); i++)
		{
			rate[i] = new double[GetNcat()];
			ratesuffstatcount[i] = new int[GetNcat()];
			ratesuffstatbeta[i] = new double[GetNcat()];
		}

		alloc = new int[GetNsite()];
		SampleRate();
	}
}


void PartitionedDGamRateProcess::Delete() 	{
	delete[] alloc;
	delete[] alpha;
	delete[] ratemult;

	for(int i = 0; i < GetNpart(); i++)
	{
		delete[] rate[i];
		delete[] ratesuffstatcount[i];
		delete[] ratesuffstatbeta[i];
	}

	delete[] rate;
	delete[] ratesuffstatcount;
	delete[] ratesuffstatbeta;

	alpha = 0;
	ratemult = 0;
	rate = 0;
	alloc = 0;
	ratesuffstatcount = 0;
	ratesuffstatbeta = 0;
	RateProcess::Delete();
}

void PartitionedDGamRateProcess::ToStream(ostream& os)	{
	for(int i = 0; i < GetNpart(); i++)
	{
		if(GetNpart() > 1)
			os << ratemult[i] << '\n';

		os << alpha[i] << '\n';
	}

	if(GetNpart() > 1)
	{
		os << alphaHyper << '\n';
		os << multHyper << '\n';
	}
}

void PartitionedDGamRateProcess::FromStream(istream& is)	{
	double tmp;
	for(int i = 0; i < GetNpart(); i++)
	{
		if(GetNpart() > 1)
			is >> ratemult[i];

		is >> tmp;
		SetAlpha(i,tmp);
	}

	if(GetNpart() > 1)
	{
		is >> alphaHyper;
		is >> multHyper;
	}
}

double PartitionedDGamRateProcess::GetMultiplierEntropy()
{
	double total = 0.0;
	for(int i = 0; i < GetNpart(); i++)
	{
		total += ratemult[i];
	}

	double ent = 0;
	for(int i = 0; i < GetNpart(); i++)
	{
		double norm = (ratemult[i] / total);
		ent -= norm * log(norm);
	}

	return ent;
}

void PartitionedDGamRateProcess::UpdateDiscreteCategories(int inpart)	{

	double* x = new double[GetNcat()];
	double* y = new double[GetNcat()];
	double lg = rnd::GetRandom().logGamma(alpha[inpart]+1.0);
	for (int i=0; i<GetNcat(); i++)	{
		x[i] = PointGamma((i+1.0)/GetNcat(), alpha[inpart], alpha[inpart]);
	}
	for (int i=0; i<GetNcat()-1; i++)	{
		y[i] = IncompleteGamma(alpha[inpart]*x[i], alpha[inpart]+1,lg);
	}
	y[GetNcat()-1] = 1.0;
	rate[inpart][0] = GetNcat() * y[0] ;
	for (int i=1; i<GetNcat(); i++)	{
		rate[inpart][i] = GetNcat() * (y[i] - y[i-1]);
	}
	delete[] x;
	delete[] y;
}

void PartitionedDGamRateProcess::SampleRate()	{
	// alpha = rnd::GetRandom().sExpo();
	for(int i = 0; i < GetNpart(); i++)
	{
		ratemult[i] = 1.0;
		alpha[i] = 1.0;
		UpdateDiscreteCategories(i);
	}

	multHyper = 1.0;
	alphaHyper = 1.0;
}

double PartitionedDGamRateProcess::Move(double tuning, int nrep)	{
	GlobalUpdateSiteRateSuffStat();

	chronorate.Start();

	GlobalUpdateRateSuffStat();

	int naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		naccepted += MoveAlphas(tuning);
	}

	if(GetNpart() > 1)
	{
		MoveMultipliers();
		MoveHyper(tuning,nrep);
	}

	chronorate.Stop();

	return ((double) naccepted) / (GetNpart() * nrep);
}

double PartitionedDGamRateProcess::NonMPIMove(double tuning, int nrep)	{
	UpdateSiteRateSuffStat();

	UpdateRateSuffStat();

	int naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		naccepted += MoveAlphas(tuning);
	}

	if(GetNpart())
		MoveMultipliers();

	return ((double) naccepted) / (nrep*GetNpart());
}

double PartitionedDGamRateProcess::MoveHyper(double tuning, int nrep)	{
	chronorate.Start();

	MoveAlphaHyper();

	int naccepted = 0;
	for (int rep=0; rep<nrep; rep++)
	{
		naccepted += MoveMultiplierHyper(tuning);
	}

	chronorate.Stop();

	return ((double) naccepted) / nrep;
}

double PartitionedDGamRateProcess::LogRatePrior()
{
	return -alphaHyper;
}

double PartitionedDGamRateProcess::LogAlphaPrior(int inpart)	{
	return -alpha[inpart];
}

double PartitionedDGamRateProcess::LogRateLikelihood(int inpart)	{
	double total = 0;
	for (int k=0; k<GetNcat(); k++)	{
		total += ratesuffstatcount[inpart][k] * log(rate[inpart][k]*ratemult[inpart]) - ratesuffstatbeta[inpart][k] * rate[inpart][k] * ratemult[inpart];
	}
	return total;
}

int PartitionedDGamRateProcess::MoveAlphas(double tuning)	{

	int naccepted = 0;
	for(int p = 0; p < GetNpart(); p++)
	{
		double bkalpha = alpha[p];
		double deltalogprob = -LogAlphaPrior(p) - LogRateLikelihood(p);
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		double newalpha = alpha[p] * e;
		SetAlpha(p, newalpha);
		deltalogprob += m + LogAlphaPrior(p) + LogRateLikelihood(p);
		bool accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted) {
			naccepted++;
		}else{
			SetAlpha(p,bkalpha);
		}
	}

	return naccepted;
}

// Gibbs alphas hyperparam move

void PartitionedDGamRateProcess::MoveAlphaHyper()
{
	double b = 0.0;
	for(int p = 0; p < GetNpart(); p++)
	{
		b += alpha[p];
	}

	alphaHyper = rnd::GetRandom().Gamma(GetNpart() + 1.0, b + 1.0);
}

// ratemult ~ iid Gamma(multHyper + 1.0, multHyper + 1.0)

double PartitionedDGamRateProcess::LogMultiplierPrior()	{

	double logProb = -multHyper;

	logProb += GetNpart()*( multHyper*log(multHyper) - rnd::GetRandom().logGamma(multHyper) );

	for(int p = 0; p < GetNpart(); p++)
	{
		logProb += (multHyper - 1.0)*log(ratemult[p]) - multHyper*ratemult[p];
	}

	return logProb;
}

// Gibbs multiplier move

void PartitionedDGamRateProcess::MoveMultipliers()
{
	for(int p = 0; p < GetNpart(); p++)
	{
		double a = 0.0;
		double b = 0.0;

		for (int k=0; k<GetNcat(); k++)	{
			a += ratesuffstatcount[p][k];
			b += ratesuffstatbeta[p][k]*rate[p][k];
		}

		ratemult[p] = rnd::GetRandom().Gamma(a + multHyper - 1.0, b + multHyper);
	}

}

// MH multiplier hyperprior move
// multHyper ~ exp(1)

int PartitionedDGamRateProcess::MoveMultiplierHyper(double tuning)	{

	int naccepted = 0;
	double bkhyper = multHyper;
	double deltalogprob = - LogMultiplierPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	multHyper *= e;
	deltalogprob += m + LogMultiplierPrior();
	bool accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
	if (accepted) {
		naccepted++;
	}else{
		multHyper = bkhyper;
	}

	return naccepted;
}


void PartitionedDGamRateProcess::GlobalUpdateRateSuffStat()	{
	assert(GetMyid() == 0);
	// MPI2
	// should ask the slaves to call their UpdateRateSuffStat
	// and then gather the statistics;
	int i,j,p,nprocs = GetNprocs(),workload = GetNcat();
	MPI_Status stat;
	MESSAGE signal = UPDATE_RATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<GetNpart(); ++i) {
		for(j=0; j<workload; j++)
		{
			ratesuffstatcount[i][j] = 0;
			ratesuffstatbeta[i][j] = 0.0;
		}
	}
#ifdef BYTE_COM
	int k,l;
	double x;
	unsigned char* bvector = new unsigned char[workload*GetNpart()*(sizeof(int)+sizeof(double))];

	for(p=1; p<nprocs; ++p) {
		MPI_Recv(bvector,workload*GetNpart()*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(i=0; i<GetNpart(); ++i) {
			for(j=0; j<workload; ++j) {
				l = 0;
				for(k=sizeof(int)-1; k>=0; --k) {
					l = (l << 8) + bvector[sizeof(int)*(i*workload + j)+k];
				}
				ratesuffstatcount[i][j] += l;
			}

			for(j=0; j<workload; ++j) {
				memcpy(&x,&bvector[sizeof(int)*workload*GetNpart()+sizeof(double)*(i*workload + j)],sizeof(double));
				ratesuffstatbeta[i][j] += x;
			}
		}
	}
	delete[] bvector;
#else
	int ivector[workload*GetNpart()];
	double dvector[workload*GetNpart()];
	for(p=1; p<nprocs; ++p) {
			MPI_Recv(ivector,workload*GetNpart(),MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			for(i=0; i<GetNpart(); ++i) {
				for(j=0; j<workload; ++j) {
					ratesuffstatcount[i][j] += ivector[i*workload + j];
				}
			}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(p=1; p<nprocs; ++p) {
			MPI_Recv(dvector,workload*GetNpart(),MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			for(i=0; i<GetNpart(); ++i) {
				for(j=0; j<workload; ++j) {
					ratesuffstatbeta[i][j] += dvector[i*workload + j];
				}
			}
	}
#endif
}

void PartitionedDGamRateProcess::UpdateRateSuffStat()	{

	for(int i=0; i<GetNpart(); ++i) {
		for (int j=0; j<GetNcat(); j++)	{
			ratesuffstatcount[i][j] = 0;
			ratesuffstatbeta[i][j] = 0.0;
		}
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		ratesuffstatcount[GetSitePart(i)][alloc[i]] += GetSiteRateSuffStatCount(i);
		ratesuffstatbeta[GetSitePart(i)][alloc[i]] += GetSiteRateSuffStatBeta(i);
	}

}	

void PartitionedDGamRateProcess::SlaveUpdateRateSuffStat()	{
	assert(GetMyid() > 0);

	UpdateRateSuffStat();

#ifdef BYTE_COM
	int n = 0;
	unsigned int j;
	unsigned char el_int[sizeof(int)],el_dbl[sizeof(double)];
	unsigned char* bvector = new unsigned char[GetNcat()*GetNpart()*(sizeof(int)+sizeof(double))];

	for(int p=0; p<GetNpart(); ++p)
		for(int i=0; i<GetNcat(); ++i) {
			convert(el_int,ratesuffstatcount[p][i]);
			for(j=0; j<sizeof(int); ++j) {
				bvector[n] = el_int[j]; n++;
			}
		}

	for(int p=0; p<GetNpart(); ++p)
		for(int i=0; i<GetNcat(); ++i) {
			convert(el_dbl,ratesuffstatbeta[p][i]);
			for(j=0; j<sizeof(double); ++j) {
				bvector[n] = el_dbl[j]; n++;
			}
		}
	MPI_Send(bvector,GetNcat()*GetNpart()*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;
#else
	int ivector[GetNcat()*GetNpart()];
	double dvector[GetNcat()*GetNpart()];

	for(int i=0; i<GetNpart(); ++i){
		memcpy(&ivector[i*GetNcat()], ratesuffstatcount[i], GetNcat()*sizeof(double));
		memcpy(&dvector[i*GetNcat()], ratesuffstatbeta[i], GetNcat()*sizeof(double));
	}

	MPI_Send(ivector,GetNcat()*GetNpart(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(dvector,GetNcat()*GetNpart(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
#endif
}
