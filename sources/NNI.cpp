
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

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


//-------------------------------------------------RG
//-------------------------------------------------
//------------------- NNI Moves -------------------
//-------------------------------------------------
//-------------------------------------------------


double PhyloProcess::GibbsNNI(double tuning, int type){
	if(!type){tuning=0;}
	GlobalRootAtRandom();
	GlobalUpdateConditionalLikelihoods(); // It have to be call !
	GlobalComputeNodeLikelihood(GetRoot());
	int success =0;
	int moves =0;

	int anumber = rnd::GetRandom().Uniform()*6+1;
	if(anumber > 1){GlobalKnit(GetRoot());}
	if(anumber > 3){GlobalKnit(GetRoot());}
	RecursiveGibbsNNI(GetRoot()->Next()->Out(),tuning,type,success,moves);
	for (int i=anumber%2+1; i; i--){GlobalKnit(GetRoot());}
	RecursiveGibbsNNI(GetRoot()->Next()->Out(),tuning,type,success,moves);
	for (int i=anumber%2+1; i; i--){GlobalKnit(GetRoot());}
	RecursiveGibbsNNI(GetRoot()->Next()->Out(),tuning,type,success,moves);
	


	GlobalComputeNodeLikelihood(GetRoot());
	// cerr << success << ' ' << moves << '\n';
	return success/(double)moves;
}


// The Local Gibbs sampling over the 3 topologies over the branches is done recursively.
void PhyloProcess::RecursiveGibbsNNI(Link* from, double tuning, int type, int& success, int& moves){


		// Forward
	if((! from->isLeaf()) and (! from->isRoot())){
		if(rnd::GetRandom().Uniform() > 0.5){GlobalKnit(from);}
		for (int i=2; i; i--){
			GlobalKnit(from);
			if(! from->Next()->Out()->isLeaf()){
				success += GlobalNNI(from,tuning,type);
				moves++;
			}
		}
	}
	// Recursion
	if(rnd::GetRandom().Uniform() > 0.5){GlobalKnit(from);}
	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		GlobalPropagateOverABranch(link->Out());
		RecursiveGibbsNNI(link->Out(),tuning,type,success,moves);
	}

	if((! from->isLeaf()) and (! from->isRoot())){
		if(rnd::GetRandom().Uniform() > 0.5){GlobalKnit(from);}
		// Backward
		for (int i=2; i; i--){
			GlobalKnit(from);
			if(! from->Next()->Out()->isLeaf()){
				success += GlobalNNI(from,tuning,0);
				moves++;
			}
		}
		GlobalPropagateOverABranch(from->Out());
	}
}






// GLOBAL NNI
int PhyloProcess::GlobalNNI(Link* from, double tuning, int type)	{


	// MPI
	MESSAGE signal = NNI;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);


	int n =0;
	if(type){
		n = 1 + (int) 5 * rnd::GetRandom().Uniform();
	}
	int args[] = {GetLinkIndex(from),n};
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);

	Link** branches;
	double logDiffPriorAndHastings = 0.0;
	if(type){
		branches = new Link*[n];
		logDiffPriorAndHastings = SendRandomBranches(from, tuning, branches, n);
	}
	// Initialise the loglikelihood vector
	double* loglikelihood = new double[3];
	loglikelihood[0]=logL;
	loglikelihood[1]=logDiffPriorAndHastings;
	loglikelihood[2]=logDiffPriorAndHastings;

	// Receive the new loglikelihood
	double* vec = new double[2];
	MPI_Status stat;
	for(int i=1; i<nprocs; ++i) {
		MPI_Recv(vec,2,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		loglikelihood[1]+=vec[0];
		loglikelihood[2]+=vec[1];
	}
	delete[] vec;

	// Sample a configuration
	int choice = rnd::GetRandom().DrawFromLogDiscreteDistribution(loglikelihood, 3);
	MPI_Bcast(&choice,1,MPI_INT,0,MPI_COMM_WORLD);
	bool success = (choice != 0);

	// Update the logL of the model
	logL=loglikelihood[choice];
	if(success){logL-=logDiffPriorAndHastings;}
	delete[] loglikelihood;

	// Restore the branch length if the initial configuration was choosen
	if(type){
		if(! success){
			for(int i=0; i<n; ++i){
				Restore(branches[i]->GetBranch());
			}
		}
		delete[] branches;
	}

	// Go to the good the configuration
	if(success){
		for(int i = choice;i>0;i--){
			GetTree()->NNIturn(from);
		}
	}
	return success;
}


// This function fill the table of link with n link pointing on the branches that will move
// It also send the movements to each of the slaves.
double PhyloProcess::SendRandomBranches(Link* from, double tuning, Link**& branches, int n){

	// Here we draw the n branches
	int* br = new int[n];
	rnd::GetRandom().DrawFromUrn (br , n, 5);

	double* m = new double[n];
	double logDiffPriorAndHastings = 0;
	for(int i = 0; i<n; i++){
		switch(br[i]){
			case 0:	branches[i] = from->Next()->Out()->Next()->Next();
				break;
			case 1: branches[i] = from->Next()->Out()->Next();
				break;
			case 2: branches[i] = from->Next()->Next();
				break;
			case 3: branches[i] = from;
				break;
			case 4: branches[i] = from->Next();
				break;
			default:cerr << "Error in PhyloProcess::SendRandomBranches " << br[i] << ' ';exit(1);
		}
		m[i] = tuning * (rnd::GetRandom().Uniform() - 0.5);
		logDiffPriorAndHastings += m[i];
		Branch* b = branches[i]->GetBranch();
		logDiffPriorAndHastings += LogBranchLengthPrior(b);
		MoveBranch(b,m[i]);
		logDiffPriorAndHastings -= LogBranchLengthPrior(b);
		br[i]=GetLinkIndex(branches[i]);
	}
	MPI_Bcast(br,n,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(m,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

	delete[] m;
	delete[] br;
	return logDiffPriorAndHastings;
}




void PhyloProcess::SlaveNNI(Link* from, int n){
	
	Link* up = from->Next();
	int* br;
	double* m;
	if(n){
		br = new int[n];
		m = new double[n];
		MPI_Bcast(br,n,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(m,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
		for(int i=0; i<n; ++i){
			Link* link = GetLinkForGibbs(br[i]);
			MoveBranch(link->GetBranch(),m[i]);
			if(link!=up){PropagateOverABranch(link);}
		}
		PropagateOverABranch(up);
	}
	


	int choice = SlaveSendNNILikelihood(from);

	

	if(n){
		if(choice == 0){
			for(int i=0; i<n; ++i){
				Link* link = GetLinkForGibbs(br[i]);
				Restore(link->GetBranch());
				if(link!=up){PropagateOverABranch(link);}
			}
		}
		delete[] br;
		delete[] m;
	}


	if(choice == 1){
		GetTree()->NNIturn(from);	
	}
	if(choice != 2){
		GetTree()->NNIturn(from);
		PropagateOverABranch(up);
	}
}




int PhyloProcess::SlaveSendNNILikelihood(Link* from){
	
	double* loglikelihood = new double[2];

	GetTree()->NNIturn(from);
	PropagateOverABranch(from->Next());
	loglikelihood[0]= ComputeNodeLikelihood(from);

	GetTree()->NNIturn(from);
	PropagateOverABranch(from->Next());
	loglikelihood[1]= ComputeNodeLikelihood(from);

	MPI_Send(loglikelihood,2,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	delete[] loglikelihood;

	int choice;
	MPI_Bcast(&choice,1,MPI_INT,0,MPI_COMM_WORLD);
	return choice;
}







void PhyloProcess::GlobalPropagateOverABranch(Link* from)	{
	assert(myid == 0);
	assert(!from->isRoot());
	if(! from->isLeaf() ){
		MESSAGE signal = BRANCHPROPAGATE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		int args[] = {GetLinkIndex(from)};
		MPI_Bcast(args,1,MPI_INT,0,MPI_COMM_WORLD);
	}
}



// Update the conditional Likelihoods associated to from
// Assume that everything else is updated
void PhyloProcess::PropagateOverABranch(const Link* from){


	if( from->Out()->isLeaf() ){
		Initialize(condlmap[0],GetData(from->Out()));
	}
	else{
		Reset(condlmap[0]);
		for (const Link* link=from->Out()->Next(); link!=from->Out(); link=link->Next())	{
			if(!link->isRoot()){
				Multiply(GetConditionalLikelihoodVector(link),condlmap[0]);
			}
		}
		Offset(condlmap[0]);
	}
	Propagate(condlmap[0],GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()));
}



void PhyloProcess::GlobalKnit(Link* from)	{

	assert(myid == 0);
	MESSAGE signal = KNIT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int args[] = {GetLinkIndex(from)};
	MPI_Bcast(args,1,MPI_INT,0,MPI_COMM_WORLD);
	from->Knit();

}


// --------------------------------------------
// --------------------------------------------
// --------------------------------------------
// --------------------------------------------
// --------------------------------------------

/*



void PhyloProcess::SlaveCheckLinkLikelihood(Link* from, int id){
	double d = ComputeNodeLikelihood(from);
	UpdateConditionalLikelihoods();
	double err = abs(d - ComputeNodeLikelihood(from));
	if(err>1e-10)
	{
		cout.flush();cout<< myid << ' ' << id << ' '  << std::setprecision(40) << err <<'\n';
		exit(1);
	}
}	



double*** PhyloProcess::SaveCondl(const Link* from){
	double*** condl = new double**[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		condl[i] = new double*[GetNrate(i)];
		for (int j=0; j<GetNrate(i); j++)	{
			condl[i][j] = new double[GetNstate(i) + 1];
			for (int k=0; k<GetNstate(i); k++)	{
				condl[i][j][k] = GetConditionalLikelihoodVector(from)[i][j][k];
			}
		}
	}
	return condl;
}


bool PhyloProcess::TestAndDeleteCondl(const Link* from, double*** condl){
	bool usefull = false;
	for (int i=sitemin; i<sitemax; i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			for (int k=0; k<GetNstate(i); k++)	{
				if(! usefull){
					usefull = (condl[i][j][k] != GetConditionalLikelihoodVector(from)[i][j][k]);
				}
			}
			delete[] condl[i][j];
		}
		delete[] condl[i];
	}
	delete[] condl;
	return usefull;
}*/
