
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
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "Parallel.h"
#include <string.h>


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* General Path Suff Stat PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixPhyloProcess::Unfold()	{

	if (condflag)	{
		cerr << "error in PhyloProcess::Unfold\n";
		exit(1);
	}
	DeleteSuffStat();
	DeleteMappings();
	ActivateSumOverRateAllocations();

	// this will in fact create only the matrices that did not already exist
	CreateMatrices();

	// this one is important
	UpdateMatrices();

	CreateCondSiteLogL();
	CreateConditionalLikelihoods();

	UpdateConditionalLikelihoods();
}

void GeneralPathSuffStatMatrixPhyloProcess::Collapse()	{

	if (! condflag)	{
		cerr << "error in PhyloProcess::Collapse\n";
		exit(1);
	}
	DrawAllocations();
	SampleNodeStates();
	if (! dataclamped)	{
		SimulateForward();
	}
	DeleteCondSiteLogL();
	DeleteConditionalLikelihoods();
	InactivateSumOverRateAllocations(ratealloc);
	SampleSubstitutionMappings(GetRoot());
	// DeleteMatrices();
	CreateSuffStat();
}

void GeneralPathSuffStatMatrixPhyloProcess::GlobalUnfold()	{

	assert(myid == 0);
	DeleteSuffStat();
	GlobalUpdateParameters();

	CreateMatrices();

	MESSAGE signal = UNFOLD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	GlobalUpdateConditionalLikelihoods();
}


void GeneralPathSuffStatMatrixPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (sitepaircount)	{
		cerr << "error in PhyloProcess::CreateSuffStat\n";
		exit(1);
	}
	siterootstate = new int[GetNsite()];
	sitepaircount = new map<pair<int,int>, int>[GetNsite()];
	sitewaitingtime = new map<int,double>[GetNsite()];
}

void GeneralPathSuffStatMatrixPhyloProcess::DeleteSuffStat()	{

	/*
	if (!sitepaircount)	{
		cerr << "error in PhyloProcess::DeleteSuffStat\n";
		// exit(1);
	}
	*/
	delete[] siterootstate;
	delete[] sitepaircount;
	delete[] sitewaitingtime;
	siterootstate = 0;
	sitepaircount = 0;
	sitewaitingtime = 0;
	PhyloProcess::DeleteSuffStat();
}

void GeneralPathSuffStatMatrixPhyloProcess::UpdateSiteProfileSuffStat()	{

	for (int i=sitemin; i<sitemax; i++)	{
		sitepaircount[i].clear();
		sitewaitingtime[i].clear();
	}

	for (int j=0; j<GetNbranch(); j++)	{
		AddSiteProfileSuffStat(siterootstate,sitepaircount,sitewaitingtime,submap[j],blarray[j],(j == 0));
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::UpdateSiteRateSuffStat()	{

	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		siteratesuffstatcount[i] = 0;
		siteratesuffstatbeta[i] = 0;
	}

	for (int j=1; j<GetNbranch(); j++)	{
		AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,submap[j],blarray[j]);
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::UpdateBranchLengthSuffStat()	{

	branchlengthsuffstatcount[0] = 0;
	branchlengthsuffstatbeta[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		int& count = branchlengthsuffstatcount[j];
		double& beta = branchlengthsuffstatbeta[j];
		count = 0;
		beta = 0;
		AddBranchLengthSuffStat(count,beta,submap[j]);
	}
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MPI 
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

	for (int i=0; i<GetNsite(); i++)	{
		sitepaircount[i].clear();
		sitewaitingtime[i].clear();
	}

	// MPI2
	// ask slaves to update siteprofilesuffstats
	// slaves should call : UpdateSiteProfileSuffStat
	// then collect all suff stats
	assert(myid == 0);
	int width,inalloc,dnalloc,smin[nprocs-1],smax[nprocs-1],iworkload[nprocs-1],dworkload[nprocs-1];
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
	inalloc = 0;
	dnalloc = 0;
	for(int i=0; i<nprocs-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (nprocs-2)) smax[i] = GetNsite();
		// workload[i] = (smax[i] - smin[i])*(GetNstate()*GetNstate()*sizeof(int) + GetNstate()*sizeof(double));
		iworkload[i] = (smax[i] - smin[i])*(GetNstate()*GetNstate() + 1);
		if (iworkload[i] > inalloc) inalloc = iworkload[i];
		dworkload[i] = (smax[i] - smin[i])*GetNstate();
		if (dworkload[i] > dnalloc) dnalloc = dworkload[i];
	}

	int* ivector = new int[inalloc];
	double* dvector = new double[dnalloc];

	int iload = GetNsite() * (GetNstate()*GetNstate() + 1);
	int dload = GetNsite() * GetNstate();
	int* iivector = new int[iload];
	double* ddvector = new double[dload];

	int im = 0;
	for(int i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,iworkload[i-1],MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		int m = 0;
		for(int j=smin[i-1]; j<smax[i-1]; ++j) {
			siterootstate[j] = ivector[m];
			iivector[im] = ivector[m];
			im++;
			m++;
			for(int k=0; k<GetNstate(); ++k) {
				for(int l=0; l<GetNstate(); ++l) {
					if (ivector[m])	{
						sitepaircount[j][pair<int,int>(k,l)] = ivector[m];
					}
					iivector[im] = ivector[m];
					m++;
					im++;
				}
			}
		}
	}

	int dm = 0;
	for(int i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,dworkload[i-1],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		int m = 0;
		for(int j=smin[i-1]; j<smax[i-1]; ++j) {
			for(int k=0; k<GetNstate(); ++k) {
				if (dvector[m])	{
					sitewaitingtime[j][k] = dvector[m];
				}
				ddvector[dm] = dvector[m];
				m++;
				dm++;
			}
		}
	}

	// checksum
	if (im != iload)	{
		cerr << "count error for gen path suff stat counts\n";
		cerr << im << '\t' << iload << '\n';
		exit(1);
	}
	if (dm != dload)	{
		cerr << "count error for gen path suff stat waiting times\n";
		cerr << dm << '\t' << dload << '\n';
		exit(1);
	}

	MPI_Bcast(iivector,iload,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(ddvector,dload,MPI_DOUBLE,0,MPI_COMM_WORLD);

	delete[] ivector;
	delete[] dvector;
	delete[] iivector;
	delete[] ddvector;
}

void GeneralPathSuffStatMatrixPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
	int iworkload = (sitemax- sitemin)*(GetNstate()*GetNstate()+1), dworkload = (sitemax - sitemin)*GetNstate();
	int* ivector = new int[iworkload];
	int m = 0;
	for(int j=sitemin; j<sitemax; ++j) {
		ivector[m] = siterootstate[j];
		m++;
		for(int k=0; k<GetNstate(); ++k) {
			for(int l=0; l<GetNstate(); ++l) {
				ivector[m] = sitepaircount[j][pair<int,int>(k,l)];
				m++;
			}
		}
	}
	if (m != iworkload)	{
		cerr << "count error\n";
		exit(1);
	}
	MPI_Send(ivector,iworkload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	double* dvector = new double[dworkload];
	m = 0;
	for(int j=sitemin; j<sitemax; ++j) {
		for(int k=0; k<GetNstate(); ++k) {
			dvector[m] = sitewaitingtime[j][k];
			m++;
		}
	}
	if (m != dworkload)	{
		cerr << "count error\n";
		exit(1);
	}
	MPI_Send(dvector,dworkload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	for (int i=0; i<GetNsite(); i++)	{
		sitepaircount[i].clear();
		sitewaitingtime[i].clear();
	}

	int iload = GetNsite() * (GetNstate()*GetNstate() + 1);
	int dload = GetNsite() * GetNstate();
	int* iivector = new int[iload];
	double* ddvector = new double[dload];

	MPI_Bcast(iivector,iload,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(ddvector,dload,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int im = 0;
	for(int j=0; j<GetNsite(); j++)	{
		siterootstate[j] = iivector[im];
		im++;
		for(int k=0; k<GetNstate(); ++k) {
			for(int l=0; l<GetNstate(); ++l) {
				if (iivector[im])	{
					sitepaircount[j][pair<int,int>(k,l)] = iivector[im];
				}
				im++;
			}
		}
	}

	int dm = 0;
	for(int j=0; j<GetNsite(); j++)	{
		for(int k=0; k<GetNstate(); ++k) {
			if (ddvector[dm])	{
				sitewaitingtime[j][k] = ddvector[dm];
			}
			dm++;
		}
	}

	// checksum
	if (im != iload)	{
		cerr << "count error for gen path suff stat counts\n";
		cerr << im << '\t' << iload << '\n';
		exit(1);
	}
	if (dm != dload)	{
		cerr << "count error for gen path suff stat waiting times\n";
		cerr << dm << '\t' << dload << '\n';
		exit(1);
	}

	delete[] ivector;
	delete[] dvector;
	delete[] iivector;
	delete[] ddvector;
}


int GeneralPathSuffStatMatrixPhyloProcess::CountMapping(int i)	{
	cerr << "in count mapping\n";
	exit(1);
	int count = 0;
	for(int k=0; k<GetNstate(i); ++k) {
		for(int l=0; l<GetNstate(i); ++l) {
			count+=sitepaircount[i][pair<int,int>(k,l)];
		}
	}
	return count;
}

