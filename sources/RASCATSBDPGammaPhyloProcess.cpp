
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <cassert>
#include "RASCATSBDPGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void RASCATSBDPGammaPhyloProcess::GlobalUpdateParameters()	{

	RASCATGammaPhyloProcess::GlobalUpdateParameters();
	MPI_Bcast(V,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATSBDPGammaPhyloProcess::SlaveUpdateParameters()	{

	RASCATGammaPhyloProcess::SlaveUpdateParameters();
	MPI_Bcast(V,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATSBDPGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

    /*
    case SITELOGCV:
        SlaveComputeSiteLogCVScore();
        break;
    */
    case SITELOGLCUTOFF:
        SlaveSetSiteLogLCutoff();
        break;
	case MIX_MOVE:
		SlaveMixMove();
		break;
	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	default:
		RASCATGammaPhyloProcess::SlaveExecute(signal);
	}
}

void RASCATSBDPGammaPhyloProcess::SlaveComputeCVScore()	{

    // determine cutoff on stick-breaking process
    // as a way to speed up computation
    // (most very low weight components in the trail of the mixture
    // don't make any meaningful contribution to the total log likelihood)
    // this will result in a slight underestimate of the fit of the model

    int ncomp = GetNcomponent();
    double totw = 0;
    while (ncomp && (totw < siteloglcutoff))    {
        ncomp--;
        totw += weight[ncomp];
    }

	sitemax = sitemin + testsitemax - testsitemin;
	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[ncomp];
		// sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<ncomp; k++) {
	// for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			PoissonSBDPProfileProcess::alloc[i] = k;
			UpdateZip(i);
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
		double max = 0;
        for (int k=0; k<ncomp; k++) {
		// for (int k=0; k<GetNcomponent(); k++)	{
			if ((!k) || (max < sitelogl[i][k]))	{
				max = sitelogl[i][k];
			}
		}
		double tot = 0;
		double totweight = 0;
        for (int k=0; k<ncomp; k++) {
		// for (int k=0; k<GetNcomponent(); k++)	{
			tot += weight[k] * exp(sitelogl[i][k] - max);
			totweight += weight[k];
		}
		total += log(tot) + max;
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;

	sitemax = bksitemax;

}

/*
void RASCATSBDPGammaPhyloProcess::SlaveComputeSiteLogCVScore()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

    // determine cutoff on stick-breaking process
    // as a way to speed up computation
    // (most very low weight components in the trail of the mixture
    // don't make any meaningful contribution to the total log likelihood)
    // this will result in a slight underestimate of the fit of the model

    int ncomp = GetNcomponent();
    double totw = 0;
    while (ncomp && (totw < siteloglcutoff))    {
        ncomp--;
        totw += weight[ncomp];
    }

	sitemax = sitemin + testsitemax - testsitemin;
	double** sitelogl = new double*[testnsite];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[ncomp];
		// sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<ncomp; k++)	{
	// for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			PoissonSBDPProfileProcess::alloc[i] = k;
            UpdateZip(i);
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double* meansitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meansitelogl[i] = 0;
	}

	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
		double max = 0;
		for (int k=0; k<ncomp; k++) {
		// for (int k=0; k<GetNcomponent(); k++)	{
			if ((!k) || (max < sitelogl[i][k]))	{
				max = sitelogl[i][k];
			}
		}
		double tot = 0;
		double totweight = 0;
		for (int k=0; k<ncomp; k++)	{
		// for (int k=0; k<GetNcomponent(); k++)	{
			tot += weight[k] * exp(sitelogl[i][k] - max);
			totweight += weight[k];
		}
		meansitelogl[i] = log(tot) + max;
		total += meansitelogl[i] ;
	}

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;
	delete[] meansitelogl;

	sitemax = bksitemax;
}
*/

void RASCATSBDPGammaPhyloProcess::SlaveComputeSiteLogL()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

    // determine cutoff on stick-breaking process
    // as a way to speed up computation
    // (most very low weight components in the trail of the mixture
    // don't make any meaningful contribution to the total log likelihood)
    // this will result in a slight underestimate of the fit of the model

    int ncomp = GetNcomponent();
    double totw = 0;
    while (ncomp && (totw < siteloglcutoff))    {
        ncomp--;
        totw += weight[ncomp];
    }

	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            sitelogl[i] = new double[ncomp];
        }
	}
	
	for (int k=0; k<ncomp; k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                PoissonSBDPProfileProcess::alloc[i] = k;
                UpdateZip(i);
            }
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                sitelogl[i][k] = sitelogL[i];
            }
		}
	}

	double* meansitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meansitelogl[i] = 0;
	}
	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            double max = 0;
            for (int k=0; k<ncomp; k++) {
                if ((!k) || (max < sitelogl[i][k]))	{
                    max = sitelogl[i][k];
                }
            }
            double tot = 0;
            double totweight = 0;
            for (int k=0; k<ncomp; k++)	{
                tot += weight[k] * exp(sitelogl[i][k] - max);
                totweight += weight[k];
            }
            meansitelogl[i] = log(tot) + max;
            total += meansitelogl[i] ;
        }
    }

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            delete[] sitelogl[i];
        }
	}
	delete[] sitelogl;
	delete[] meansitelogl;

}

double RASCATSBDPGammaPhyloProcess::GlobalGetSiteSteppingLogLikelihood(int site, int nrep0, int restore) {

    if (nrep0)  {
        if (! empcount) {
            cerr << "error: in IS site log likelihood, empcount not created\n";
            exit(1);
        }
    }

    int nrep_per_proc = nrep0 / (GetNprocs()-1);
    if (nrep0 % (GetNprocs() - 1))  {
        nrep_per_proc ++;
    }
    int nrep = nrep_per_proc * (GetNprocs()-1);

    MESSAGE signal = STEPPINGSITELOGL;
    MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    int param[3];
    param[0] = site;
    param[1] = nrep_per_proc;
    param[2] = restore;
    MPI_Bcast(param,3,MPI_INT,0,MPI_COMM_WORLD);
    double ret = 0;
    if (nrep0)  {
        ret = GlobalGetSiteSteppingLogLikelihoodIS(site, nrep, restore);
    }
    else    {
        ret = GlobalGetSiteSteppingLogLikelihoodNonIS(site, nrep, restore);
    }
    return ret;
}

void RASCATSBDPGammaPhyloProcess::SlaveGetSiteSteppingLogLikelihood()    {

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

    int param[3];
    MPI_Bcast(param,3,MPI_INT,0,MPI_COMM_WORLD);
    int site = param[0];
    int nrep = param[1];
    int restore = param[2];

    if (nrep)   {
        SlaveGetSiteSteppingLogLikelihoodIS(site, nrep, restore);
    }
    else    {
        SlaveGetSiteSteppingLogLikelihoodNonIS(site, nrep, restore);
    }
}


double RASCATSBDPGammaPhyloProcess::GlobalGetSiteSteppingLogLikelihoodIS(int site, int nrep, int restore) {

    int oldalloc = PoissonSBDPProfileProcess::alloc[site];

    UpdateOccupancyNumbers();

    double master_logl[2*GetNprocs()];
    double slave_logl[2];
    slave_logl[0] = 0;
    slave_logl[1] = 0;

    int master_alloc[GetNprocs()];
    int slave_alloc = -1;

    double master_profile[GetDim()*GetNprocs()];
    double slave_profile[GetDim()];
    for (int k=0; k<GetDim(); k++)  {
        slave_profile[k] = 0;
    }

    MPI_Gather(slave_logl, 2, MPI_DOUBLE, master_logl, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&slave_alloc, 1, MPI_INT, master_alloc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(slave_profile, GetDim(), MPI_DOUBLE, master_profile, GetDim(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double max1 = 0;
    for (int i=1; i<GetNprocs(); i++)   {
        if ((!max1) || (max1 < master_logl[2*i+1]))    {
            max1 = master_logl[2*i+1];
        }
    }
    double tot1 = 0;
    double post1[GetNprocs()-1];
    for (int i=1; i<GetNprocs(); i++)   {
        double tmp = exp(master_logl[2*i+1]-max1);
        post1[i-1] = tmp;
        tot1 += tmp;
    }

    double L1 = log(tot1) + max1;

    for (int i=1; i<GetNprocs(); i++)   {
        post1[i-1] /= tot1;
    }
    int procalloc1 = rnd::GetRandom().FiniteDiscrete(GetNprocs()-1, post1) + 1;

    double max0 = 0;
    for (int i=1; i<GetNprocs(); i++)   {
        if ((!max0) || (max0 < master_logl[2*i]))    {
            max0 = master_logl[2*i];
        }
    }

    double tot0 = 0;
    double post0[GetNprocs()-1];
    for (int i=1; i<GetNprocs(); i++)   {
        double tmp = exp(master_logl[2*i]-max0);
        post0[i-1] = tmp;
        tot0 += tmp;
    }

    double L0 = log(tot0/nrep) + max0;

    for (int i=1; i<GetNprocs(); i++)   {
        post0[i-1] /= tot0;
    }
    int procalloc0 = rnd::GetRandom().FiniteDiscrete(GetNprocs()-1, post0) + 1;

    double w0[GetNcomponent()];
    double totw0 = 0;
    // double totw1 = 0;
    for (int k=0; k<GetNcomponent(); k++) {
        if (occupancy[k])   {
            w0[k] = 0;
            // totw1 += weight[k];
        }
        else    {
            w0[k] = weight[k];
            totw0 += weight[k];
        }
    }
    for (int k=0; k<GetNcomponent(); k++) {
        w0[k] /= totw0;
    }
    double l0 = 0;
    double l1 = 1;
    double logl = l1;

    if (GetNcomponent() > GetNOccupiedComponent())  {
        double max = (L0 > L1) ? L0 : L1;
        l0 = totw0 * exp(L0-max);
        l1 = exp(L1-max);
        // l1 = totw1 * exp(L1-max);
        logl = log(l0 + l1) + max;
    }

    if (! restore)  {
        if ((l1 + l0)*rnd::GetRandom().Uniform() < l1)  {
            int newalloc = master_alloc[procalloc1];
            RemoveSite(site, oldalloc);
            AddSite(site, newalloc);
        }
        else    {
            int newalloc = rnd::GetRandom().FiniteDiscrete(GetNcomponent(), w0);
            RemoveSite(site, oldalloc);
            AddSite(site, newalloc);
            for (int k=0; k<GetDim(); k++)  {
                profile[newalloc][k] = master_profile[GetDim()*procalloc0+k];
            }
        }
        ResampleEmptyProfiles();
        GlobalUpdateParameters();
    }
    return logl;
}

void RASCATSBDPGammaPhyloProcess::SlaveGetSiteSteppingLogLikelihoodIS(int site, int nrep, int restore)    {

    int bkalloc = PoissonSBDPProfileProcess::alloc[site];
    double bkprofile[GetDim()];
    for (int k=0; k<GetDim(); k++)  {
        bkprofile[k] = profile[bkalloc][k];
    }

    UpdateOccupancyNumbers();
    int nocc = GetNOccupiedComponent();
	int width = nocc / (GetNprocs()-1);
    int r = nocc % (GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
    int s = 0;
	for(int i=0; i<GetNprocs()-1; i++) {
		smin[i] = s;
        if (i < r)  {
            s += width + 1;
        }
        else    {
            s += width;
        }
        smax[i] = s;
    }
    if (s != nocc)  {
        cerr << "error: nocc checksum\n";
        exit(1);
    }

    // get the range of components for this slave
    int ncomp = 0;
    int kmin = 0;
    int kmax = 0;
    for (int k=0; k<GetNcomponent(); k++)   {
        if (occupancy[k])   {
            ncomp++;
            if (ncomp == smin[myid-1])    {
                kmin = k;
            }
            if (ncomp == smax[myid-1])    {
                kmax = k;
            }
        }
    }
    int krange = kmax - kmin;

    double master_logl[2*GetNprocs()];
    double slave_logl[2];
    slave_logl[0] = 0;
    slave_logl[1] = 0;

    int master_alloc[GetNprocs()];
    int slave_alloc = -1;

    double master_profile[GetDim()*GetNprocs()];
    double slave_profile[GetDim()];

    double sitelogl1[krange];
    for (int k=0; k<krange; k++)    {
        sitelogl1[k] = 0;
    }

    for (int k=kmin; k<kmax; k++)	{
        if (occupancy[k])   {
            PoissonSBDPProfileProcess::alloc[site] = k;
            // UpdateZip(site);
            // PrepareSiteLogLikelihood(site);
            double tmp = SiteLogLikelihood(site);
            sitelogl1[k-kmin] = tmp;
        }
    }

    double max1 = 0;
    for (int k=kmin; k<kmax; k++)	{
        if (occupancy[k])   {
            if ((!max1) || (max1 < sitelogl1[k-kmin]))	{
                max1 = sitelogl1[k-kmin];
            }
        }
    }

    double post1[krange];
    double tot1= 0;
    for (int k=kmin; k<kmax; k++)   {
        if (occupancy[k])   {
            post1[k-kmin] = weight[k] * exp(sitelogl1[k-kmin] - max1);
            tot1 += post1[k-kmin];
        }
        else    {
            post1[k-kmin] = 0;
        }
    }

    slave_logl[1] = log(tot1) + max1;

    for (int k=kmin; k<kmax; k++)   {
        post1[k-kmin] /= tot1;
    }
    slave_alloc = rnd::GetRandom().FiniteDiscrete(krange, post1) + kmin;

    if (GetNcomponent() > GetNOccupiedComponent())  {

        vector<vector<double>> isprofile(nrep, vector<double>(GetDim(), 0));
        vector<double> islogl(nrep, 0);

        for (int rep=0; rep<nrep; rep++)    {
            SampleEmpiricalStat(profile[bkalloc], empcount + site*GetDim());
            for (int k=0; k<GetDim(); k++)  {
                isprofile[rep][k] = profile[bkalloc][k];
            }
            // UpdateZip(site);
            // PrepareSiteLogLikelihood(site);
            double tmp = SiteLogLikelihood(site);
            islogl[rep] = tmp;
            islogl[rep] += RASCATGammaPhyloProcess::LogStatPrior(profile[bkalloc]) - EmpiricalLogStatPrior(profile[bkalloc], empcount + site*GetDim());
        }

        double max0 = 0;
        for (int rep=0; rep<nrep; rep++)    {
            if ((!max0) || (max0 < islogl[rep]))  {
                max0 = islogl[rep];
            }
        }

        double tot0 = 0;
        double isw[nrep];
        for (int rep=0; rep<nrep; rep++)	{
            isw[rep] = exp(islogl[rep] - max0);
            tot0 += isw[rep];
        }
        slave_logl[0] = log(tot0) + max0;
        for (int rep=0; rep<nrep; rep++)	{
            isw[rep] /= tot0;
        }
        int weightindex = rnd::GetRandom().FiniteDiscrete(nrep, isw);
        for (int k=0; k<GetDim(); k++)  {
            slave_profile[k] = isprofile[weightindex][k];
        }
    }

    MPI_Gather(slave_logl, 2, MPI_DOUBLE, master_logl, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&slave_alloc, 1, MPI_INT, master_alloc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(slave_profile, GetDim(), MPI_DOUBLE, master_profile, GetDim(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (restore)    {
        PoissonSBDPProfileProcess::alloc[site] = bkalloc;
        for (int k=0; k<GetDim(); k++)  {
            profile[bkalloc][k] = bkprofile[k];
        }
        UpdateZip(site);
    }
}

double RASCATSBDPGammaPhyloProcess::GlobalGetSiteSteppingLogLikelihoodNonIS(int site, int nrep, int restore) {

    int oldalloc = PoissonSBDPProfileProcess::alloc[site];

    double master_logl[GetNprocs()];
    double slave_logl;

    int master_alloc[GetNprocs()];
    int slave_alloc = -1;

    MPI_Gather(&slave_logl, 1, MPI_DOUBLE, master_logl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&slave_alloc, 1, MPI_INT, master_alloc, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double max = 0;
    for (int i=1; i<GetNprocs(); i++)   {
        if (master_alloc[i] != -1)    {
            if ((!max) || (max < master_logl[i]))    {
                max = master_logl[i];
            }
        }
    }
    double tot = 0;
    double post[GetNprocs()-1];
    for (int i=1; i<GetNprocs(); i++)   {
        if (master_alloc[i] != -1)    {
            double tmp = exp(master_logl[i]-max);
            post[i-1] = tmp;
            tot += tmp;
        }
        else    {
            post[i-1] = 0;
        }
    }
    if (! tot)  {
        cerr << "error in stepping logl: total likelihood is 0\n";
        exit(1);
    }

    double L = log(tot) + max;

    for (int i=1; i<GetNprocs(); i++)   {
        post[i-1] /= tot;
    }
    int procalloc = rnd::GetRandom().FiniteDiscrete(GetNprocs()-1, post) + 1;

    if (! restore)  {
        int newalloc = master_alloc[procalloc];
        RemoveSite(site, oldalloc);
        AddSite(site, newalloc);
        GlobalUpdateParameters();
    }
    return L;
}

void RASCATSBDPGammaPhyloProcess::SlaveGetSiteSteppingLogLikelihoodNonIS(int site, int nrep, int restore)    {

    int bkalloc = PoissonSBDPProfileProcess::alloc[site];

	int width = GetNcomponent() / (GetNprocs()-1);
    int r = GetNcomponent() % (GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
    int s = 0;
	for(int i=0; i<GetNprocs()-1; i++) {
		smin[i] = s;
        if (i < r)  {
            s += width + 1;
        }
        else    {
            s += width;
        }
        smax[i] = s;
    }
    if (s != GetNcomponent())  {
        cerr << "error: ncomp checksum\n";
        exit(1);
    }

    int kmin = smin[myid-1];
    int kmax = smax[myid-1];
    int krange = kmax - kmin;

    double sitelogl[krange];
    for (int k=kmin; k<kmax; k++)	{
        PoissonSBDPProfileProcess::alloc[site] = k;
        double tmp = SiteLogLikelihood(site);
        sitelogl[k-kmin] = tmp;
    }

    double max = 0;
    for (int k=kmin; k<kmax; k++)	{
        if ((!max) || (max < sitelogl[k-kmin]))	{
            max = sitelogl[k-kmin];
        }
    }

    double post[krange];
    double tot= 0;
    for (int k=kmin; k<kmax; k++)   {
        post[k-kmin] = weight[k] * exp(sitelogl[k-kmin] - max);
        tot += post[k-kmin];
    }

    double slave_logl = 0;
    if (tot > 0)    {
        slave_logl = log(tot) + max;
    }

    int slave_alloc = -1;
    if (tot > 0)    {
        for (int k=kmin; k<kmax; k++)   {
            post[k-kmin] /= tot;
        }
        slave_alloc = rnd::GetRandom().FiniteDiscrete(krange, post) + kmin;
    }

    double master_logl[GetNprocs()];
    int master_alloc[GetNprocs()];

    MPI_Gather(&slave_logl, 1, MPI_DOUBLE, master_logl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&slave_alloc, 1, MPI_INT, master_alloc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (restore)    {
        PoissonSBDPProfileProcess::alloc[site] = bkalloc;
        UpdateZip(site);
    }
}

