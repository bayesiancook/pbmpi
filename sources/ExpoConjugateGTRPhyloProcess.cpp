
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
#include "ExpoConjugateGTRPhyloProcess.h"
#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugateGTRPhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void ExpoConjugateGTRPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (siteprofilesuffstatcount)	{
		cerr << "error in ExpoConjugateGTRPhyloProcess::CreateSuffStat\n";
		cerr << myid << '\n';
		exit(1);
	}
	allocsiteprofilesuffstatcount = new int[GetNsite()*GetDim()];
	allocsiteprofilesuffstatbeta = new double[GetNsite()*GetDim()];
	tmpcount = new int[GetNsite()*GetDim()];
	tmpbeta = new double[GetNsite()*GetDim()];
	siteprofilesuffstatcount = new int*[GetNsite()];
	siteprofilesuffstatbeta = new double*[GetNsite()];
	// for (int i=sitemin; i<sitemax; i++)	{
	for (int i=0; i<GetNsite(); i++)	{
		siteprofilesuffstatcount[i] = allocsiteprofilesuffstatcount + i*GetDim();
		siteprofilesuffstatbeta[i] = allocsiteprofilesuffstatbeta + i*GetDim();
	}
}

void ExpoConjugateGTRPhyloProcess::DeleteSuffStat()	{

	if (siteprofilesuffstatcount)	{
		// for (int i=sitemin; i<sitemax; i++)	{
		/*
		for (int i=0; i<GetNsite(); i++)	{
			delete[] siteprofilesuffstatcount[i];
			delete[] siteprofilesuffstatbeta[i];
		}
		*/
		delete[] siteprofilesuffstatcount;
		delete[] siteprofilesuffstatbeta;
		siteprofilesuffstatcount = 0;
		siteprofilesuffstatbeta = 0;
		delete[] allocsiteprofilesuffstatcount;
		delete[] allocsiteprofilesuffstatbeta;
		allocsiteprofilesuffstatcount = 0;
		allocsiteprofilesuffstatbeta = 0;
		delete[] tmpcount;
		delete[] tmpbeta;
		tmpcount = 0;
		tmpbeta = 0;
	}
	PhyloProcess::DeleteSuffStat();
}

void ExpoConjugateGTRPhyloProcess::UpdateRRSuffStat()	{

	for (int k=0; k<GetNrr(); k++)	{
		rrsuffstatcount[k] = 0;
		rrsuffstatbeta[k] = 0;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		// AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,submap[j],blarray[j]);
		AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateSiteRateSuffStat()	{

	// cerr << "in update site rate : " << GetTotalLength() << '\n';
	for (int i=sitemin; i<sitemax; i++)	{
		siteratesuffstatcount[i] = 0;
		siteratesuffstatbeta[i] = 0;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		// AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,submap[j],blarray[j]);
		AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateBranchLengthSuffStat()	{

	branchlengthsuffstatcount[0] = 0;
	branchlengthsuffstatbeta[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		int& count = branchlengthsuffstatcount[j];
		double& beta = branchlengthsuffstatbeta[j];
		count = 0;
		beta = 0;
		// AddBranchLengthSuffStat(count,beta,submap[j]);
		AddBranchLengthSuffStat(count,beta,submap[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateSiteProfileSuffStat()	{

	for (int i=sitemin; i<sitemax; i++)	{
		for (int k=0; k<GetDim(); k++)	{
			siteprofilesuffstatcount[i][k] = 0;
			siteprofilesuffstatbeta[i][k] = 0;
		}
	}
	for (int j=0; j<GetNbranch(); j++)	{
		// AddSiteProfileSuffStat(siteprofilesuffstatcount,siteprofilesuffstatbeta,submap[j],blarray[j], (j == 0));
		AddSiteProfileSuffStat(siteprofilesuffstatcount,siteprofilesuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

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
		workload[i] = (smax[i] - smin[i])*GetGlobalNstate();
		if (workload[i] > nalloc) nalloc = workload[i];
	}

	int ivector[nalloc];
	double dvector[nalloc];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,workload[i-1],MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		l = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			for(k=0; k<GetGlobalNstate(); ++k) {
				siteprofilesuffstatcount[j][k] = ivector[l]; l++;
			}
		}
	}
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,workload[i-1],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		l = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			for(k=0; k<GetGlobalNstate(); ++k) {
				siteprofilesuffstatbeta[j][k] = dvector[l]; l++;
			}
		}
	}

	MPI_Bcast(allocsiteprofilesuffstatcount,GetNsite()*GetGlobalNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocsiteprofilesuffstatbeta,GetNsite()*GetGlobalNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void ExpoConjugateGTRPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
	int i,j,workload = (sitemax - sitemin)*GetGlobalNstate();

	int k = 0,ivector[workload];
	for(i=sitemin; i<sitemax; ++i) {
		for(j=0; j<GetGlobalNstate(); ++j) {
			ivector[k] = siteprofilesuffstatcount[i][j]; k++;
		}
	}
	MPI_Send(ivector,workload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	double dvector[workload];
	k = 0;
	for(i=sitemin; i<sitemax; ++i) {
		for(j=0; j<GetGlobalNstate(); ++j) {
			dvector[k] = siteprofilesuffstatbeta[i][j]; k++;
		}
	}
	MPI_Send(dvector,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Bcast(allocsiteprofilesuffstatcount,GetNsite()*GetGlobalNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocsiteprofilesuffstatbeta,GetNsite()*GetGlobalNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void ExpoConjugateGTRPhyloProcess::GlobalUpdateRRSuffStat()	{

	// MPI2
	// should send message to slaves for updating their rrsuffstats
	// by calling UpdateRRSuffStat();
	// then collect all suff stats
	//
	// suff stats are contained in 2 arrays
	// int* rrsuffstatcount
	// double* rrsuffstatbeta

	// should be summed over all slaves (reduced)
	assert(myid == 0);
	int i,j,workload = Nrr;
	MPI_Status stat;
	MESSAGE signal = UPDATE_RRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<workload; ++i) {
		rrsuffstatcount[i] = 0;
		rrsuffstatbeta[i] = 0.0;
	}

	int ivector[workload];
	double dvector[workload];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,workload,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(ivector,workload,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<workload; ++j) {
			rrsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,workload,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(dvector,workload,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<workload; ++j) {
			rrsuffstatbeta[j] += dvector[j];
		}
	}

	MPI_Bcast(rrsuffstatcount,Nrr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void ExpoConjugateGTRPhyloProcess::SlaveUpdateRRSuffStat()	{

	UpdateRRSuffStat();
	int workload = Nrr;

	MPI_Send(rrsuffstatcount,workload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(rrsuffstatbeta,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Bcast(rrsuffstatcount,Nrr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

int ExpoConjugateGTRPhyloProcess::GlobalCountMapping()	{

	GlobalUpdateSiteProfileSuffStat();
	return PhyloProcess::GlobalCountMapping();
}

int ExpoConjugateGTRPhyloProcess::CountMapping()	{

	int total = 0;	
	for(int i = sitemin; i < sitemax; i++){
		total += CountMapping(i);
	}
	return total;
}

int ExpoConjugateGTRPhyloProcess::CountMapping(int i)	{

	const int* tmp = GetSiteProfileSuffStatCount(i);
	int total = 0;
	for (int k=0; k<GetNstate(i); k++)	{
		total += tmp[k];
	}
	total--;
	return total;
}

