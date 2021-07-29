
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "DGamRateProcess.h"
#include "Random.h"
#include "IncompleteGamma.h"

#include <cassert>
#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* DGamRateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void DGamRateProcess::Create(int innsite, int inncat)	{
	if (! rate)	{
		RateProcess::Create(innsite);
		Ncat = inncat;
		rate = new double[GetNcat()];
		alloc = new int[GetNsite()];
		ratesuffstatcount = new int[GetNcat()];
		ratesuffstatbeta = new double[GetNcat()];
        empalpha = 1.0;
        empbeta = 1.0;
		// SampleRate();
	}
}

void DGamRateProcess::Delete() 	{
	delete[] rate;
	delete[] alloc;
	delete[] ratesuffstatcount;
	delete[] ratesuffstatbeta;
	rate = 0;
	alloc = 0;
	ratesuffstatcount = 0;
	ratesuffstatbeta = 0;
	RateProcess::Delete();
}

void DGamRateProcess::ToStream(ostream& os)	{
	os << alpha << '\n';
}

void DGamRateProcess::FromStream(istream& is)	{
	double tmp;
	is >> tmp;
	SetAlpha(tmp);
}

void DGamRateProcess::UpdateDiscreteCategories()	{

	double* x = new double[GetNcat()];
	double* y = new double[GetNcat()];
	double lg = rnd::GetRandom().logGamma(alpha+1.0);
	for (int i=0; i<GetNcat(); i++)	{
		x[i] = PointGamma((i+1.0)/GetNcat(),alpha,alpha);
	}
	for (int i=0; i<GetNcat()-1; i++)	{
		y[i] = IncompleteGamma(alpha*x[i],alpha+1,lg);
	}
	y[GetNcat()-1] = 1.0;
	rate[0] = GetNcat() * y[0];
	for (int i=1; i<GetNcat(); i++)	{
		rate[i] = GetNcat() * (y[i] - y[i-1]);
	}
	delete[] x;
	delete[] y;
}

void DGamRateProcess::SetRateEmpiricalPrior(double inalpha, double inbeta)  {
    empalpha = inalpha;
    empbeta = inbeta;
}

void DGamRateProcess::SampleRate()	{
	// alpha = rnd::GetRandom().sExpo();
	alpha = 1;
	UpdateDiscreteCategories();
}

void DGamRateProcess::PriorSampleRate()	{
	alpha = rnd::GetRandom().Gamma(ratefrac + (1-ratefrac)*empalpha, ratefrac + (1-ratefrac)*empbeta);
	// alpha = rnd::GetRandom().sExpo();
	UpdateDiscreteCategories();
}

double DGamRateProcess::LogRatePrior()	{
    if (ratefrac == 1)  {
        return -alpha;
    }
    double a = ratefrac + (1-ratefrac)*empalpha;
    double b = ratefrac + (1-ratefrac)*empbeta;
    return a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(alpha) - b*alpha;
}

double DGamRateProcess::RateSuffStatLogProb()	{
	double total = 0;
	for (int k=0; k<GetNcat(); k++)	{
		total += ratesuffstatcount[k] * log(rate[k]) - ratesuffstatbeta[k] * rate[k];
	}
	return total;
}

double DGamRateProcess::MoveAlpha(double tuning, int nrep)	{

	GlobalUpdateRateSuffStat();
	int naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double bkalpha = alpha;
		double deltalogprob = -LogRatePrior() - RateSuffStatLogProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		double newalpha = alpha * e;
		SetAlpha(newalpha);
		deltalogprob += m + LogRatePrior() + RateSuffStatLogProb();
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			SetAlpha(bkalpha);
		}
	}
	return ((double) naccepted) / nrep;
}

double DGamRateProcess::NonMPIMoveAlpha(double tuning, int nrep)	{

	UpdateRateSuffStat();
	int naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double bkalpha = alpha;
		double deltalogprob = -LogRatePrior() - RateSuffStatLogProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		double newalpha = alpha * e;
		SetAlpha(newalpha);
		deltalogprob += m + LogRatePrior() + RateSuffStatLogProb();
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			SetAlpha(bkalpha);
		}
	}
	return ((double) naccepted) / nrep;
}


void DGamRateProcess::GlobalUpdateRateSuffStat()	{
	assert(GetMyid() == 0);
	// MPI2
	// should ask the slaves to call their UpdateRateSuffStat
	// and then gather the statistics;
	int i,j,nprocs = GetNprocs(),workload = GetNcat();
	MPI_Status stat;
	MESSAGE signal = UPDATE_RATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<workload; ++i) {
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}
	int ivector[workload];
	double dvector[workload];
        for(i=1; i<nprocs; ++i) {
                MPI_Recv(ivector,workload,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(j=0; j<workload; ++j) {
                        ratesuffstatcount[j] += ivector[j];                      
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(i=1; i<nprocs; ++i) {
                MPI_Recv(dvector,workload,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(j=0; j<workload; ++j) {
                        ratesuffstatbeta[j] += dvector[j]; 
                }
        }
}

void DGamRateProcess::UpdateRateSuffStat()	{

	for (int i=0; i<GetNcat(); i++)	{
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		ratesuffstatcount[alloc[i]] += GetSiteRateSuffStatCount(i);
		ratesuffstatbeta[alloc[i]] += GetSiteRateSuffStatBeta(i);
	}

}	

void DGamRateProcess::SlaveUpdateRateSuffStat()	{
	assert(GetMyid() > 0);

	UpdateRateSuffStat();

	MPI_Send(ratesuffstatcount,GetNcat(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(ratesuffstatbeta,GetNcat(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}	
