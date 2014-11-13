
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
#include "AAMutSelSiteSpecificPhyloProcess.h"
#include "Parallel.h"
#include <string.h>


// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
// are model dependent
// should be implemented in .cpp file
void AAMutSelSiteSpecificPhyloProcess::SlaveUpdateParameters()	{

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nnucrr = GetNnucrr(),nnucstat = 4,k = 0;
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1); // check if these last terms are correct in this context...
	ni = 1 + ProfileProcess::GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[i];
	}
	for(i=0; i<nnucrr; ++i) {
		nucrr[i] = dvector[nbranch+i];
	}
	for(i=0; i<nnucstat; ++i) {
		nucstat[i] = dvector[nbranch+nnucrr+i];
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[1+nbranch+nnucrr+nnucstat+i];
	}
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[nbranch+nnucrr+nnucstat+L2+k];
			k++;
		}
		weight[i] = dvector[nbranch+nnucrr+nnucstat+L2+k];
		k++;
	}
	Ncomponent = ivector[0];
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		FiniteProfileProcess::alloc[i] = ivector[1+i];
	}
	//GetBranchLengthsFromArray();
	delete[] dvector;
	delete[] ivector;
	// this one is really important
	// in those cases where new components have appeared, or some old ones have disappeared
	// during allocation move on the master node.
	// 
	// note that CreateMatrices() in fact creates only those that are not yet allocated
	// and also deletes those that are now obsolete
	CreateMatrices();
	UpdateMatrices();
}


void AAMutSelSiteSpecificPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

	/*
	case PRINT_TREE:
		SlavePrintTree();
		break;
	*/
	case REALLOC_MOVE:
		SlaveIncrementalFiniteMove();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void AAMutSelSiteSpecificPhyloProcess::GlobalUpdateParameters() {
	// MPI2
	// should send the slaves the relevant information
	// about model parameters
	// for this model, should broadcast
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	assert(myid == 0);
	int i,j,nnucrr,nnucstat,nbranch = GetNbranch(),ni,nd,L1,L2,k = 0;
	nnucrr = GetNnucrr();
	nnucstat = 4;	
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1);  // check if these last terms are correct in this context...
	ni = 1 + ProfileProcess::GetNsite(); // 1 for the number of componenets, and the rest for allocations
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	// First we assemble the vector of doubles for distribution
	for(i=0; i<nbranch; ++i) {
		dvector[i] = blarray[i];
	}

	for(i=0; i<nnucrr; ++i) {
		dvector[nbranch+i] = nucrr[i];
	}
	for(i=0; i<nnucstat; ++i) {
		dvector[nbranch+nnucrr+i] = nucstat[i];
	}

	for (int i=0; i<GetDim(); i++)	{
		dvector[nbranch+nnucrr+nnucstat+i] = dirweight[i];
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[nbranch+nnucrr+nnucstat+L2+k] = profile[i][j];
			k++;
		}
		dvector[nbranch+nnucrr+nnucstat+L2+k] = weight[i];
		k++;
	}

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		ivector[1+i] = FiniteProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

