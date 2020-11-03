
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
#include "Parallel.h"
#include <string>
#include "PoissonPhyloProcess.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Poisson PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PoissonPhyloProcess::Create(Tree* intree, SequenceAlignment* indata)	{
	if (! zipdata)	{
		truedata = indata;
		zipdata = new ZippedSequenceAlignment(truedata);
		PhyloProcess::Create(intree,zipdata,indata->GetNstate());
		truedata->GetEmpiricalFreq(empfreq);
		CreateZip();
		// CreateNodeStates();
		// CreateMappings();
		// condflag = false;
	}
}

void PoissonPhyloProcess::Delete()	{
	if (zipdata)	{
		// DeleteMappings();
		// DeleteNodeStates();
		PhyloProcess::Delete();
		delete zipdata;
		zipdata = 0;
	}
}

void PoissonPhyloProcess::Collapse()	{

	if (! condflag)	{
		cerr << "error in PoissonPhyloProcess::Collapse\n";
		exit(1);
	}
	DrawAllocations();
	SampleNodeStates();
	if (! dataclamped)	{
		if (rateprior)	{
			DrawAllocationsFromPrior();
		}
		if (profileprior)	{
			DrawProfileFromPrior();
		}
		SimulateForward();
	}
	DeleteCondSiteLogL();
	DeleteConditionalLikelihoods();
	InactivateSumOverRateAllocations(ratealloc);
	FillMissingMap();
	SampleSubstitutionMappings(GetRoot());
	CreateSuffStat();
	PoissonUpdateSiteProfileSuffStat();
}

void PoissonPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (siteprofilesuffstatcount)	{
		cerr << "error in PoissonPhyloProcess::CreateSuffStat\n";
		exit(1);
	}
	allocsiteprofilesuffstatcount = new int[GetNsite()*GetDim()];
	siteprofilesuffstatcount = new int*[GetNsite()];
	// for (int i=sitemin; i<sitemax; i++)	{
	for (int i=0; i<GetNsite(); i++)	{
		siteprofilesuffstatcount[i] = allocsiteprofilesuffstatcount + i*GetDim();
		// siteprofilesuffstatcount[i] = new int[GetDim()];
	}
}

void PoissonPhyloProcess::DeleteSuffStat()	{

	if (siteprofilesuffstatcount)	{
		// for (int i=sitemin; i<sitemax; i++)	{
		/*
		for (int i=0; i<GetNsite(); i++)	{
			delete[] siteprofilesuffstatcount[i];
		}
		*/
		delete[] siteprofilesuffstatcount;
		siteprofilesuffstatcount = 0;
		delete[] allocsiteprofilesuffstatcount;
		allocsiteprofilesuffstatcount = 0;
	}
	PhyloProcess::DeleteSuffStat();
}

void PoissonPhyloProcess::UpdateBranchLengthSuffStat()	{

	// double R = GetNactiveSite() * GetMeanRate();
	branchlengthsuffstatbeta[0] = 0;
	branchlengthsuffstatcount[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		// branchlengthsuffstatbeta[j] = R;
		branchlengthsuffstatbeta[j] = 0;
		branchlengthsuffstatcount[j] = 0;
		AddBranchLengthSuffStat(branchlengthsuffstatcount[j],branchlengthsuffstatbeta[j],submap[j],missingmap[j]);
	}
}

void PoissonPhyloProcess::UpdateSiteRateSuffStat()	{
	// double totallength = GetTotalLength();
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		siteratesuffstatcount[i] = 0;
		siteratesuffstatbeta[i] = 0;
		// siteratesuffstatbeta[i] = totallength;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,blarray[j],submap[j],missingmap[j]);
	}
}

int PoissonPhyloProcess::RecursiveUpdateSiteProfileSuffStat(const Link* from, int site)	{

	int state = -1;
	if (from->isLeaf())	{
		state = GetData(from)[site];
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			int tmp = RecursiveUpdateSiteProfileSuffStat(link->Out(),site);
			if (tmp != -1)	{
				if ((state != -1) && (state != tmp))	{
					cerr << "error in PoissonPhyloProcess::RecursiveUpdateSiteProfileSuffStat: state should be identical\n";
					cerr << state << '\t' << tmp << '\t' << GetZipSize(site) << '\t' << GetOrbitSize(site) << '\n';
					exit(1);
				}
				state = tmp;
			}
		}
	}
	if (state != -1)	{
		BranchSitePath* path = submap[GetBranchIndex(from->GetBranch())][site];
		if (from->isRoot() || path->GetNsub())	{
			if ((GetZipSize(site) != GetOrbitSize(site)) && (state == GetOrbitSize(site)))	{
				cerr << "error in PoissonPhyloProcess::RecursiveUpdateSiteProfileSuffStat: state should be observed at tips\n";
				exit(1);
			}
			int truestate = GetStateFromZip(site,state);
			siteprofilesuffstatcount[site][truestate]++;
			state = -1;
		}
	}
	return state;
}


void PoissonPhyloProcess::UpdateSiteProfileSuffStat()	{
}

void PoissonPhyloProcess::PoissonUpdateSiteProfileSuffStat()	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			siteprofilesuffstatcount[i][k] = 0;
		}
		RecursiveUpdateSiteProfileSuffStat(GetRoot(),i);
	}
}

void PoissonPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

	// MPI2
	// ask slaves to update siteprofilesuffstats
	// slaves should call : UpdateSiteProfileSuffStat
	// then collect all suff stats
	assert(myid == 0);
	int i,j,k,l,width,nalloc,smin[nprocs-1],smax[nprocs-1],workload[nprocs-1];
	MPI_Status stat;
	MESSAGE signal = UPDATE_SPROFILE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// suff stats are contained in 2 arrays
	// int** siteprofilesuffstatcount
	// double** siteprofilesuffstatbeta
	// [site][state]

	// each slave computes its array for sitemin <= site < sitemax
	// thus, one just needs to gather all arrays into the big master array 0 <= site < Nsite
	// (gather)
	width = GetNsite()/(nprocs-1);
	nalloc = 0;
	for(i=0; i<nprocs-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (nprocs-2)) smax[i] = GetNsite();
		workload[i] = (smax[i] - smin[i])*GetDim();
		if (workload[i] > nalloc) nalloc = workload[i];
	}
	int ivector[nalloc];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,workload[i-1],MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		l = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			for(k=0; k<GetDim(); ++k) {
				siteprofilesuffstatcount[j][k] = ivector[l]; l++;
			}
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(allocsiteprofilesuffstatcount,GetNsite()*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
}

void PoissonPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
	int i,j,workload = (sitemax - sitemin)*GetDim();
	int k = 0,ivector[workload];
	for(i=sitemin; i<sitemax; ++i) {
		for(j=0; j<GetDim(); ++j) {
			ivector[k] = siteprofilesuffstatcount[i][j]; k++;
		}
	}
	MPI_Send(ivector,workload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(allocsiteprofilesuffstatcount,GetNsite()*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
}

/*
void PoissonPhyloProcess::UpdateBranchLengthSuffStat()	{
	if (! GetMyid())	{
		cerr << "error in update branch length suff stat: master\n";
		exit(1);
	}
	double R = (GetSiteMax() - GetSiteMin()) * GetMeanRate();
	// double R = GetNsite() * GetMeanRate();
	branchlengthsuffstatbeta[0] = 0;
	branchlengthsuffstatcount[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		branchlengthsuffstatbeta[j] = R;
		branchlengthsuffstatcount[j] = 0;
		AddBranchLengthSuffStat(branchlengthsuffstatcount[j],submap[j]);
	}
}

void PoissonPhyloProcess::UpdateSiteRateSuffStat()	{
	double totallength = GetTotalLength();
	for (int i=sitemin; i<sitemax; i++)	{
		siteratesuffstatcount[i] = 0;
		siteratesuffstatbeta[i] = totallength;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		AddSiteRateSuffStat(siteratesuffstatcount,submap[j]);
	}
}

void PoissonPhyloProcess::UpdateSiteProfileSuffStat()	{

	for (int i=sitemin; i<sitemax; i++)	{
		for (int k=0; k<GetDim(); k++)	{
			siteprofilesuffstatcount[i][k] = 0;
		}
	}
	for (int j=0; j<GetNbranch(); j++)	{
		AddSiteProfileSuffStat(siteprofilesuffstatcount,submap[j],(j==0));
	}
}

*/

/*
void PoissonSubstitutionProcess::ChooseTrueStates(BranchSitePath** patharray, int* nodestateup, int* nodestatedown, bool root)	{
	for (int i=sitemin; i<sitemax; i++)	{
		int tmp = nodestateup[i];
		if (root || patharray[i]->GetNsub())	{
			tmp = GetRandomStateFromZip(i,patharray[i]->GetFinalState());
		}
		cerr  << ' '<< i  << ' '<< patharray[i]->GetNsub() << ' '<< tmp  << '\n';
		nodestatedown[i] = tmp;
	}
}

*/

void PoissonPhyloProcess::GlobalSetTestData()	{

	testnsite = testdata->GetNsite();
	int* tmp = new int[testnsite * GetNtaxa()];
	testdata->GetDataVector(tmp);
	/*
	for (int i=0; i<testnsite*GetNtaxa(); i++)	{
		tmp[i] = -1;
	}
	*/

	MESSAGE signal = SETTESTDATA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);

	delete[] tmp;

	// GlobalCollapse();
	// DeleteZip();
    zipdata->ComputeZipArrays();
    /*
	delete zipdata;
	zipdata = new ZippedSequenceAlignment(truedata);
    data = zipdata;
    */
	// CreateZip();
	// GlobalUnfold();
	
}

void PoissonPhyloProcess::SlaveSetTestData()	{

	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	int* tmp = new int[testnsite * GetNtaxa()];
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	
	SetTestSiteMinAndMax();
	truedata->SetTestData(testnsite,sitemin,testsitemin,testsitemax,tmp);

	delete[] tmp;

	Collapse();
	// DeleteZip();
    zipdata->ComputeZipArrays();
    /*
	delete zipdata;
	zipdata = new ZippedSequenceAlignment(truedata);
    data = zipdata;
    */
	// CreateZip();
	Unfold();
}

void PoissonPhyloProcess::RecursiveUnzipBranchSitePath(const Link* from){
	for (const Link* link=from->Next(); link!=from; link=link->Next()){
		RecursiveUnzipBranchSitePath(link->Out());
		UnzipBranchSitePath(submap[GetBranchIndex(link->GetBranch())], GetStates(link->GetNode()), GetStates(link->Out()->GetNode()));
	}
}


void PoissonPhyloProcess::SlaveWriteMappings(){
	SampleTrueNodeStates(GetRoot());
	RecursiveUnzipBranchSitePath(GetRoot());
	PhyloProcess::SlaveWriteMappings();
}


void PoissonPhyloProcess::SampleTrueNodeStates(const Link* from)	{
	
	if (from->isRoot())	{
		ChooseTrueStates(submap[0],GetStates(from->Out()->GetNode()),GetStates(from->GetNode()),from->isRoot());
	}
	else{
		ChooseTrueStates(submap[GetBranchIndex(from->GetBranch())],GetStates(from->Out()->GetNode()),GetStates(from->GetNode()),from->isRoot());
	}

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleTrueNodeStates(link->Out());
	}
}


/*
void PoissonPhyloProcess::GlobalSetTestData()	{

	ZippedSequenceAlignment* ziptestdata = new ZippedSequenceAlignment(testdata);
	testnsite = testdata->GetNsite();
	int* tmp = new int[testnsite * GetNtaxa()];
	ziptestdata->GetDataVector(tmp);

	MESSAGE signal = SETTESTDATA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);

	delete[] tmp;
}

void PoissonPhyloProcess::SlaveSetTestData()	{

	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	int* tmp = new int[testnsite * GetNtaxa()];
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	
	SetTestSiteMinAndMax();
	zipdata->SetTestData(testnsite,sitemin,testsitemin,testsitemax,tmp);

	delete[] tmp;
}
*/
