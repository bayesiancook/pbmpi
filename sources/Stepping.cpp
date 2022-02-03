
/********************

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"

#include "Random.h"
#include "PhyloProcess.h"
#include <string>

#include <cassert>
#include "Parallel.h"

#include "TexTab.h"

bool PhyloProcess::ActiveSite(int i)  {
    return (minsitecutoff == -1) || ((steppingrank[i] >= minsitecutoff) && (steppingrank[i] < maxsitecutoff));
}


void PhyloProcess::GlobalPrepareStepping(string name, int size, int rand)   {
	assert(myid == 0);
    if (bkdata) {
        cerr << "error in SlaveBackupData: bkdata already exists\n";
        exit(1);
    }

	MESSAGE signal = PREPARESTEPPING;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    bkdata = new SequenceAlignment(GetData());

    steppingrank = new int[GetNsite()];
    if (size)   {
        ifstream is((name + ".siteranks").c_str());
        for (int i=0; i<GetNsite(); i++)    {
            is >> steppingrank[i];
        }
    }
    else    {
        if (rand)   {
            rnd::GetRandom().DrawFromUrn(steppingrank, GetNsite(), GetNsite());
        }
        else    {
            for (int i=0; i<GetNsite(); i++)    {
                steppingrank[i] = i;
            }
        }
        ofstream os((name + ".siteranks").c_str());
        for (int i=0; i<GetNsite(); i++)    {
            os << steppingrank[i] << '\t';
        }
        os << '\n';
    }
    MPI_Bcast(steppingrank,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

}

void PhyloProcess::SlavePrepareStepping()	{
    if (bkdata) {
        cerr << "error in SlaveBackupData: bkdata already exists\n";
        exit(1);
    }
    bkdata = new SequenceAlignment(GetData());
    steppingrank = new int[GetNsite()];
	MPI_Bcast(steppingrank,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
    CreateSiteConditionalLikelihoods();
}

void PhyloProcess::GlobalSetSteppingFraction(int cutoff1, int cutoff2)    {
	MESSAGE signal = SETSTEPPINGFRAC;
    int cutoff[2];
    cutoff[0] = cutoff1;
    cutoff[1] = cutoff2;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(cutoff,2,MPI_INT,0,MPI_COMM_WORLD);
    SetSteppingFraction(cutoff1, cutoff2);
}

void PhyloProcess::SlaveSetSteppingFraction()    {
    int cutoff[2];
	MPI_Bcast(cutoff,2,MPI_INT,0,MPI_COMM_WORLD);
    SetSteppingFraction(cutoff[0], cutoff[1]);
}

void PhyloProcess::SetSteppingFraction(int cutoff1, int cutoff2)  {

    for (int i=0; i<GetNsite(); i++)    {
        if ((steppingrank[i] >= cutoff1) && (steppingrank[i] < cutoff2))  {
            for (int j=0; j<GetNtaxa(); j++)    {
                GetData()->SetState(j, i, bkdata->GetState(j,i));
            }
        }
        else    {
            for (int j=0; j<GetNtaxa(); j++)    {
                GetData()->SetState(j, i, -1);
            }
        }
    }
    minsitecutoff = cutoff1;
    maxsitecutoff = cutoff2;
    // fill missing map
    FillMissingMap();
}

void PhyloProcess::GlobalSetEmpiricalFrac(double infrac)    {
	MESSAGE signal = EMPIRICALFRAC;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&infrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    SetEmpiricalFrac(infrac);
}

void PhyloProcess::SlaveSetEmpiricalFrac()  {
    double frac;
	MPI_Bcast(&frac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    SetEmpiricalFrac(frac);
}

/*
void PhyloProcess::GlobalCreateSiteDataStructures() {
	MESSAGE signal = CREATESITE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    // CreateSiteConditionalLikelihoods();
}

void PhyloProcess::SlaveCreateSiteDataStructures()  {
    CreateSiteConditionalLikelihoods();
}

void PhyloProcess::GlobalDeleteSiteDataStructures() {
	MESSAGE signal = DELETESITE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    // DeleteSiteConditionalLikelihoods();
}

void PhyloProcess::SlaveDeleteSiteDataStructures()  {
    DeleteSiteConditionalLikelihoods();
}
*/

void PhyloProcess::CreateSiteConditionalLikelihoods()	{

	if (! sitecondlmap)	{
		sitecondlmap = new double**[GetNlink()];
		for (int j=0; j<GetNlink(); j++)	{
			sitecondlmap[j] = new double*[GetMaxNrate()];
			for (int k=0; k<GetMaxNrate(); k++)	{
				sitecondlmap[j][k] = new double[GetGlobalNstate()+1];
			}
		}
	}
}

void PhyloProcess::DeleteSiteConditionalLikelihoods()	{

	if (sitecondlmap)	{
		for (int j=0; j<GetNlink(); j++)	{
			for (int k=0; k<GetMaxNrate(); k++)	{
				delete[] sitecondlmap[j][k];
			}
			delete[] sitecondlmap[j];
		}
		delete[] sitecondlmap;
		sitecondlmap = 0;
	}
}

double PhyloProcess::SiteLogLikelihood(int site)	{

	PrepareSiteLogLikelihood(site);
    if (! SumOverRateAllocations())  {
        cerr << "error in PhyloProcess::SiteLogLikelihood(site): sum over rate allocations not activated\n";
        exit(1);
    }
	// SiteActivateSumOverRateAllocation(site);
	SitePostOrderPruning(site,GetRoot());
	SiteMultiplyByStationaries(site,sitecondlmap[0]);
	double ret = SiteComputeLikelihood(site,sitecondlmap[0]);
	// SiteInactivateSumOverRateAllocation(site);
    return ret;
}

void PhyloProcess::SitePostOrderPruning(int site, const Link* from)	{

	if (from->isLeaf())	{
		SiteInitialize(site,sitecondlmap[0],GetData(from)[site]);
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			SitePostOrderPruning(site,link->Out());
			SitePropagate(site,sitecondlmap[0],sitecondlmap[GetLinkIndex(link)],GetLength(link->GetBranch()));
		}
		SiteReset(site,sitecondlmap[0]);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			SiteMultiply(site,sitecondlmap[GetLinkIndex(link)],sitecondlmap[0]);
		}
		SiteOffset(site,sitecondlmap[0]);
	}
	if (from->isRoot())	{
	}	
}

double PhyloProcess::GlobalGetSiteSteppingLogLikelihood(int site, int nrep, int restore)   {

	MESSAGE signal = STEPPINGSITELOGL;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&site,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	int maxwidth = 0;
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
		if (maxwidth < (smax[i] - smin[i]))	{
			maxwidth = smax[i] - smin[i];
		}
	}

    int slave = -1;
    for(int i=1; i<GetNprocs(); ++i) {
        if ((site >= smin[i-1]) && (site < smax[i-1]))  {
            slave = i;
        }
    }

    if (slave == -1)    {
        cerr << "error: slave proc not found\n";
        exit(1);
    }

    double ret;
	MPI_Status stat;
    MPI_Recv(&ret,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
    return ret;
}

void PhyloProcess::SlaveGetSiteSteppingLogLikelihood()  {

    int site, nrep;
	MPI_Bcast(&site,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

    double ret = 0;
    if ((site >= sitemin) && (site < sitemax))  {
        ret = SiteLogLikelihood(site);
        MPI_Send(&ret,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }
}

