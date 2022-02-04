
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

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PhyloProcess::Unfold()	{
	if (condflag)	{
		cerr << "error in PhyloProcess::Unfold\n";
		exit(1);
	}
	DeleteSuffStat();
	DeleteMappings();
	ActivateSumOverRateAllocations();
	CreateCondSiteLogL();
	CreateConditionalLikelihoods();
	UpdateConditionalLikelihoods();
}

void PhyloProcess::Collapse()	{

	if (! condflag)	{
		cerr << "error in PhyloProcess::Collapse\n";
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
}

void PhyloProcess::CreateMappings()	{

	for (int j=0; j<GetNbranch(); j++)	{
		submap[j] = 0;
	}
}

void PhyloProcess::DeleteMappings()	{

	for (int j=0; j<GetNbranch(); j++)	{
		if (submap[j])	{
			for (int i=sitemin; i<sitemax; i++)	{
				delete submap[j][i];
			}
			delete[] submap[j];
			submap[j] = 0;
		}
	}
}

void PhyloProcess::CreateSuffStat()	{

	if (! siteratesuffstatcount)	{
		if (GetMyid())	{
			siteratesuffstatcount = new int[GetNsite()];
			siteratesuffstatbeta = new double[GetNsite()];
		}
	}
	if (! branchlengthsuffstatcount)	{
		branchlengthsuffstatcount = new int[GetNbranch()];
		branchlengthsuffstatbeta = new double[GetNbranch()];
	}
	activesuffstat = true;
}

void PhyloProcess::DeleteSuffStat()	{

	if (GetMyid())	{
		delete[] siteratesuffstatcount;
		delete[] siteratesuffstatbeta;
		siteratesuffstatcount = 0;
		siteratesuffstatbeta = 0;
	}
	delete[] branchlengthsuffstatcount;
	delete[] branchlengthsuffstatbeta;
	branchlengthsuffstatcount = 0;
	branchlengthsuffstatbeta = 0;
	activesuffstat = false;
}

void PhyloProcess::CreateNodeStates()	{

	for (int j=0; j<GetNnode(); j++)	{
		nodestate[j] = new int[GetNsite()];
	}
}

void PhyloProcess::DeleteNodeStates()	{

	for (int j=0; j<GetNnode(); j++)	{
		delete[] nodestate[j];
	}
}
void PhyloProcess::CreateConditionalLikelihoods()	{

	// do not create for leaves
	if (! condflag)	{
		for (int j=0; j<GetNlink(); j++)	{
			condlmap[j] =  CreateConditionalLikelihoodVector();
		}
	}
	condflag = true;
}

void PhyloProcess::GlobalResetAllConditionalLikelihoods()  {
	assert(myid == 0);
	MESSAGE signal = RESETALL;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveResetAllConditionalLikelihoods()	{
    for (int j=0; j<GetNlink(); j++)	{
        Reset(condlmap[j], false, true);
    }
}

void PhyloProcess::DeleteConditionalLikelihoods()	{

	if (condflag)	{
		for (int j=0; j<GetNlink(); j++)	{
			DeleteConditionalLikelihoodVector(condlmap[j]);
		}
	}
	condflag = false;
}

void PhyloProcess::UpdateConditionalLikelihoods()	{
	PostOrderPruning(GetRoot(),condlmap[0]);

	// not necessary
	MultiplyByStationaries(condlmap[0]);
	ComputeLikelihood(condlmap[0]);

	PreOrderPruning(GetRoot(),condlmap[0]);
}

void PhyloProcess::GlobalCheckLikelihood()	{

	vector<double> logl;
	GlobalRecursiveComputeLikelihood(GetRoot(),0,logl);
	double max = 0;
	for (unsigned int i=0; i<logl.size(); i++)	{
		double tmp = fabs(logl[i] - logl[0]);
		if (max < tmp)	{
			max = tmp;
		}
	}
	if (max > 1e-10)	{
		cerr << "error in check likelihoods\n";
		cerr << max << '\n';
		cerr.precision(25);
		for (unsigned int i=0; i<logl.size(); i++)	{
			cerr << logl[i] << '\n';
		}
		exit(1);
	}
}

void PhyloProcess::CheckLikelihood()	{

	vector<double> logl;
	RecursiveComputeLikelihood(GetRoot(),0,logl);
	double max = 0;
	for (unsigned int i=0; i<logl.size(); i++)	{
		double tmp = fabs(logl[i] - logl[0]);
		if (max < tmp)	{
			max = tmp;
		}
	}
	if (max > 1e-10)	{
		cerr << "error in check likelihoods\n";
		cerr << max << '\n';
		cerr.precision(25);
		for (unsigned int i=0; i<logl.size(); i++)	{
			cerr << logl[i] << '\n';
		}
		exit(1);
	}
}

double PhyloProcess::ComputeNodeLikelihood(const Link* from, int auxindex)	{

	if (! myid)	{
		cerr << "error : master doing slave's work\n";
		exit(1);
	}
	double*** aux = 0;
	bool localaux = false;
	if (auxindex != -1)	{
		aux = condlmap[auxindex];
	}
	else	{
		localaux = true;
		aux = CreateConditionalLikelihoodVector();
	}

	if (from->isLeaf())	{
		Initialize(aux,GetData(from));
	}
	else	{
		Reset(aux);
	}
	if (! from->isRoot())	{
		Multiply(GetConditionalLikelihoodVector(from),aux);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! link->isRoot())	{
			Multiply(GetConditionalLikelihoodVector(link),aux);
		}
	}
	MultiplyByStationaries(aux);
	double lnL = ComputeLikelihood(aux);
	if (localaux)	{
		DeleteConditionalLikelihoodVector(aux);
	}
	if (std::isnan(lnL))	{
		cerr << "in PhyloProcess::ComputeNodeLikelihood: nan\n";
		exit(1);
	}
	return lnL;
}

void PhyloProcess::PostOrderPruning(const Link* from, double*** aux)	{

	if (from->isLeaf())	{
        Initialize(aux,GetData(from));
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			PostOrderPruning(link->Out(),aux);
			Propagate(aux,GetConditionalLikelihoodVector(link),GetLength(link->GetBranch()));
		}
		Reset(aux);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			Multiply(GetConditionalLikelihoodVector(link),aux);
		}
		Offset(aux);
	}
	if (from->isRoot())	{
		// copy aux into GetConditionalLikelihoodVector(root) ?
		// or aux IS the conditional likelihood vector of the root ?
	}	
}

void PhyloProcess::PreOrderPruning(const Link* from, double*** aux)	{

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Reset(aux);
		for (const Link* link2=link->Next(); link2!=link; link2=link2->Next())	{
			if (! link2->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link2),aux);
			}
		}
		// Here, in principle
		// should be done even if link->Out()->isLeaf()
		// in order for all the conditional likelihood vectors, including those at the leaves, to be updated
		// but in practice, the leaf likelihood vectors are not used anyway (and they represent half of the whole set of likelihood vectors)
		// so not computing them saves 50% CPU time
		if (! link->Out()->isLeaf())	{
			Propagate(aux,GetConditionalLikelihoodVector(link->Out()),GetLength(link->GetBranch()));
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! link->Out()->isLeaf())	{
			PreOrderPruning(link->Out(),aux);
		}
	}
}

void PhyloProcess::GlobalRecursiveComputeLikelihood(const Link* from, int auxindex, vector<double>& logl)	{

	double lnL = GlobalComputeNodeLikelihood(from,auxindex);
	logl.push_back(lnL);
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// WARNING: preorder pruning does not update leaf condtional likelihood vectors (not necessary in the present case)
		// so the following will issue an error message if tried on leaf
		if (! link->Out()->isLeaf())	{
			GlobalRecursiveComputeLikelihood(link->Out(),auxindex,logl);
		}
	}
}

void PhyloProcess::RecursiveComputeLikelihood(const Link* from, int auxindex, vector<double>& logl)	{

	double lnL = ComputeNodeLikelihood(from,auxindex);
	// double lnL = GlobalComputeNodeLikelihood(from,auxindex);
	logl.push_back(lnL);
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// WARNING: preorder pruning does not update leaf condtional likelihood vectors (not necessary in the present case)
		// so the following will issue an error message if tried on leaf
		if (! link->Out()->isLeaf())	{
			RecursiveComputeLikelihood(link->Out(),auxindex,logl);
		}
	}
}

void PhyloProcess::SampleNodeStates()	{
	SampleNodeStates(GetRoot(),condlmap[0]);
}


void PhyloProcess::SimulateForward()	{
	if (rootprior)	{
		ChooseStatesAtEquilibrium(GetStates(GetRoot()->GetNode()));
	}
	RecursiveSimulateForward(GetRoot());
}

void PhyloProcess::RecursiveSimulateForward(const Link* from)	{
	
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SimuPropagate(GetStates(from->GetNode()),GetStates(link->Out()->GetNode()),GetLength(link->GetBranch()));
		RecursiveSimulateForward(link->Out());
	}
}


void PhyloProcess::SampleNodeStates(const Link* from, double*** aux)	{
	
	if (from->isLeaf())	{
		Initialize(aux,GetData(from));
	}
	else	{
		Reset(aux,true);
	}
	// make product of conditional likelihoods around node
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Multiply(GetConditionalLikelihoodVector(link),aux,true);
	}
	if (!from->isRoot())	{
		Multiply(GetConditionalLikelihoodVector(from),aux,true);
	}
	MultiplyByStationaries(aux,true);
	// let substitution process choose states based on this vector
	// this should collapse the vector into 1s and 0s
	ChooseStates(aux,GetStates(from->GetNode()));

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// propagate forward
		Propagate(aux,GetConditionalLikelihoodVector(link->Out()),GetLength(link->GetBranch()),true);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleNodeStates(link->Out(),aux);
	}
}

void PhyloProcess::SampleSubstitutionMappings(const Link* from)	{

	if (from->isRoot())	{
		submap[0] = SampleRootPaths(GetStates(from->GetNode()));
	}
	else	{
		submap[GetBranchIndex(from->GetBranch())] = SamplePaths(GetStates(from->Out()->GetNode()), GetStates(from->GetNode()), GetLength(from->GetBranch()));
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleSubstitutionMappings(link->Out());
	}
}

// MPI master functions
double PhyloProcess::BranchLengthMove(double tuning)	{

	// uses condlmap[0] as auxiliary variable
	int n = 0;
	double total = RecursiveBranchLengthMove(GetRoot(),tuning,n);
	return total / n;
}

// assumes aux contains the product of incoming likelihoods
double PhyloProcess::RecursiveBranchLengthMove(const Link* from, double tuning, int& n)	{

	// uses condlmap[0] as auxiliary variable
	double total = 0;

	if (! from->isRoot())	{
		total += LocalBranchLengthMove(from,tuning);
		n++;
	}
	
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		GlobalReset(0);
		// Reset(condlmap[0]);
		for (const Link* link2=link->Next(); link2!=link; link2=link2->Next())	{
			if (! link2->isRoot())	{
				GlobalMultiply(link2,0);
				// Multiply(GetConditionalLikelihoodVector(link2),condlmap[0]);
			}
		}
		total += RecursiveBranchLengthMove(link->Out(),tuning,n);
	}

	if (from->isLeaf())	{
		GlobalInitialize(0,from);
		// Initialize(condlmap[0],GetData(from));
	}
	else	{
		GlobalReset(0);
		// Reset(condlmap[0]);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (! link->isRoot())	{
				GlobalMultiply(link,0);
				// Multiply(GetConditionalLikelihoodVector(link),condlmap[0]);
			}
		}
	}
	
	if (! from->isRoot())	{
		total += LocalBranchLengthMove(from->Out(),tuning); // Pourquoi 'Out' ??? 
		n++;
	}

	return total;
};

double PhyloProcess::LocalBranchLengthMove(const Link* from, double tuning)	{

	// uses condlmap[0] as auxiliary variable

	double currentloglikelihood = logL;
	double currentlogprior = LogBranchLengthPrior(from->GetBranch());
	double loghastings = GlobalProposeMove(from->GetBranch(),tuning);
	// double loghastings = ProposeMove(from->GetBranch(),tuning);

	GlobalPropagate(0,from,GetLength(from->GetBranch()));
	// Propagate(aux,GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()));

	double newloglikelihood = GlobalComputeNodeLikelihood(from);
	// double newloglikelihood = ComputeNodeLikelihood(from);
	double newlogprior = LogBranchLengthPrior(from->GetBranch());
	double delta = newlogprior + newloglikelihood - currentlogprior - currentloglikelihood + loghastings;
	
	int accepted = (log(rnd::GetRandom().Uniform()) < delta);
	if (!accepted)	{
		GlobalRestore(from->GetBranch());
		// Restore(from->GetBranch());
		GlobalPropagate(0,from,GetLength(from->GetBranch()));
		// Propagate(aux,GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()));
		GlobalComputeNodeLikelihood(from);
		// ComputeNodeLikelihood(from);
		// not useful: done by ComputeNodeLikelihood(from) just above
		// logL = currentloglikelihood;
	}
	return (double) accepted;
}

// MPI master functions
double PhyloProcess::NonMPIBranchLengthMove(double tuning)	{

	UpdateConditionalLikelihoods();
	// uses condlmap[0] as auxiliary variable
	int n = 0;
	double total = RecursiveNonMPIBranchLengthMove(GetRoot(),tuning,n);
	return total / n;
}

// assumes aux contains the product of incoming likelihoods
double PhyloProcess::RecursiveNonMPIBranchLengthMove(const Link* from, double tuning, int& n)	{

	// uses condlmap[0] as auxiliary variable
	double total = 0;

	if (! from->isRoot())	{
		total += LocalNonMPIBranchLengthMove(from,tuning);
		n++;
	}
	//This update should be in the previous "if" ?
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Reset(condlmap[0]);
		for (const Link* link2=link->Next(); link2!=link; link2=link2->Next())	{
			if (! link2->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link2),condlmap[0]);
			}
		}
		total += RecursiveNonMPIBranchLengthMove(link->Out(),tuning,n);
	}

	if (from->isLeaf())	{
		Initialize(condlmap[0],GetData(from));
	}
	else	{
		Reset(condlmap[0]);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (! link->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link),condlmap[0]);
			}
		}
	}
	
	if (! from->isRoot())	{
		total += LocalNonMPIBranchLengthMove(from->Out(),tuning);
		n++;
	}

	return total;
};

double PhyloProcess::LocalNonMPIBranchLengthMove(const Link* from, double tuning)	{

	// uses condlmap[0] as auxiliary variable

	double currentloglikelihood = logL;
	double currentlogprior = LogBranchLengthPrior(from->GetBranch());
	double loghastings = ProposeMove(from->GetBranch(),tuning);

	Propagate(condlmap[0],GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()));

	double newloglikelihood = ComputeNodeLikelihood(from);
	double newlogprior = LogBranchLengthPrior(from->GetBranch());
	double delta = newlogprior + newloglikelihood - currentlogprior - currentloglikelihood + loghastings;
	
	int accepted = (log(rnd::GetRandom().Uniform()) < delta);
	if (!accepted)	{
		Restore(from->GetBranch());
		Propagate(condlmap[0],GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()));
		ComputeNodeLikelihood(from);
		// not useful: done by ComputeNodeLikelihood(from) just above
		// logL = currentloglikelihood;
	}
	return (double) accepted;
}

double PhyloProcess::MoveTopo(int spr, int nni){
	double success = 0;

	if (size >= topoburnin)	{
		success += GibbsSPR(spr);
		for(int i=0; i<nni; i++){
			success += GibbsNNI(0.1,1);
		}
	}
	return success;
}

double PhyloProcess::GibbsSPR(int nrep)	{
	// useless, assuming that preceding move maintains conditinal likelihoods correctly updated
	GlobalUpdateConditionalLikelihoods();
	/*
	GlobalComputeNodeLikelihood(GetRoot()->Next());
	if (logL >  -2000)	{
		cerr << "error before gibbs\n";
		exit(1);
	}
	*/

	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		naccepted += GibbsSPR();
	}
	// necessary, assuming that following move assumes likelihoods are updated
	/*
	GlobalComputeNodeLikelihood(GetRoot()->Next());
	if (logL >  -2000)	{
		cerr << "error after gibbs\n";
		exit(1);
	}
	*/

	/*
	MPI_Status stat;
	MESSAGE signal = BCAST_TREE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	GlobalBroadcastTree();
	*/
	
	GlobalUpdateConditionalLikelihoods();



	return naccepted / nrep;
}

int PhyloProcess::GibbsSPR()	{

	Link* up = 0;
	Link* down = 0;
	GlobalRootAtRandom();
	GetTree()->DrawSubTree(down,up);

	if (down->isRoot())	{
		cerr << "down is root\n";
		return 0;
	}
	if (up->isRoot())	{
		cerr << "up is root\n";
		return 0;
	}
	int sizebefore = GetTree()->GetSize();
	int subtreesize = GetTree()->GetSize(down);

	Link* fromdown = GlobalDetach(down,up);
	
	int sizeafter = GetTree()->GetSize();

	if (sizebefore-subtreesize-sizeafter)	{
		cerr << "error in gibbs spr: non matching subtree sizes\n";
		exit(1);
	}

	GlobalUpdateConditionalLikelihoods();
	
	GlobalGibbsSPRScan(down,up,loglarray);
	map<pair<Link*,Link*>, double> loglmap;
	int n = 0;
	RecursiveGibbsFillMap(GetRoot(),GetRoot(),loglmap,loglarray,n);

	double max = 0;
	int j = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (std::isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (std::isinf(i->second))	{
			cerr << "inf log prob in gibbs\n";
			exit(1);
		}
		if ((i==loglmap.begin()) || (max < i->second))	{
			max = i->second;
		}
		j++;
	}
	if (j != 2*sizeafter-3)	{
		cerr << "error in gibbs: number of cases not matching 2p-3\n";
		cerr << j << '\t' << sizeafter << '\n';
		exit(1);
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		total += exp(i->second - max);
		if (std::isinf(total))	{
			cerr << "error in gibbs: inf\n";
		}
		if (std::isnan(total))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = exp(i->second-max);
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		cumul += exp(i->second-max);
	}
	
	int accepted = (fromdown != i->first.first);
	if (i->second - max < -40)	{
		cerr << "suspicious choice in gibbs : " << i->second - max << '\t' << exp(i->second - max) << '\n';
		exit(1);
	}
	GlobalAttach(down,up,i->first.first,i->first.second);
	GlobalUpdateConditionalLikelihoods();
	return accepted;
}

void PhyloProcess::RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, double* loglarray, int& n)	{

	if (! from->isRoot())	{
		GetTree()->Attach(down,up,from,fromup);
		double*** aux = condlmap[0];
		Reset(aux);
		for (const Link* link=up->Next(); link!=up; link=link->Next())	{
			if (link->isRoot())	{
				cerr << "ROOT\n";
				exit(1);
			}
			Multiply(GetConditionalLikelihoodVector(link),aux);
		}
		Propagate(aux,GetConditionalLikelihoodVector(up->Out()),GetLength(up->GetBranch()));
		double logl = ComputeNodeLikelihood(up->Out(),0);
		if (n >= GetNbranch())	{
			cerr << "branch overflow\n";
			exit(1);
		}
		loglarray[n] = logl;
		n++;
		GetTree()->Detach(down,up);
		// Link* tmp1 = GetTree()->Detach(down,up);
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveGibbsSPRScan(link->Out(),trailer,down,up,loglarray,n);
		trailer = trailer->Next();
	}
}

void PhyloProcess::RecursiveGibbsFillMap(Link* from, Link* fromup, map<pair<Link*,Link*>,double>& loglmap, double* loglarray, int& n)	{

	if (! from->isRoot())	{
		if (n >= GetNbranch())	{
			cerr << "branch overflow\n";
			exit(1);
		}
		loglmap[pair<Link*,Link*>(from,fromup)] = loglarray[n];
		n++;
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveGibbsFillMap(link->Out(),trailer,loglmap,loglarray,n);
		trailer = trailer->Next();
	}
}

double PhyloProcess::NonMPIGibbsSPR()	{

	GetTree()->RootAtRandom();

	Link* up = 0;
	Link* down = 0;

	GetTree()->DrawSubTree(down,up);
	
	int sizebefore = GetTree()->GetSize();
	int subtreesize = GetTree()->GetSize(down);
	Link* fromdown = GetTree()->Detach(down,up);
	
	int sizeafter = GetTree()->GetSize();

	if (sizebefore-subtreesize-sizeafter)	{
		cerr << "error in gibbs spr: non matching subtree sizes\n";
		exit(1);
	}

	UpdateConditionalLikelihoods();
	
	map<pair<Link*,Link*>, double> loglmap;
	RecursiveNonMPIGibbsSPRScan(GetRoot(),GetRoot(),down,up,loglmap);

	double max=0;
	int j = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if ((i==loglmap.begin()) || (max < i->second))	{
			max = i->second;
		}
		j++;
	}
	if (j != 2*sizeafter-3)	{
		cerr << "error in gibbs: number of cases not matching 2p-3\n";
		cerr << j << '\t' << sizeafter << '\n';
		exit(1);
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		total += exp(i->second - max);
		if (std::isinf(total))	{
			cerr << "error in gibbs: inf\n";
		}
		if (std::isnan(total))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = exp(i->second-max);
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		cumul += exp(i->second-max);
	}
	
	int accepted = (fromdown != i->first.first);
	if (i->second - max < -40)	{
		cerr << "suspicious choice in gibbs : " << i->second - max << '\t' << exp(i->second - max) << '\n';
		exit(1);
	}
	GetTree()->Attach(down,up,i->first.first,i->first.second);
	UpdateConditionalLikelihoods();
	return accepted;
}

void PhyloProcess::RecursiveNonMPIGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap)	{

	if (! from->isRoot())	{
		GetTree()->Attach(down,up,from,fromup);
		double*** aux = condlmap[0];
		Reset(aux);
		for (const Link* link=up->Next(); link!=up; link=link->Next())	{
			if (link->isRoot())	{
				cerr << "ROOT\n";
				exit(1);
			}
			Multiply(GetConditionalLikelihoodVector(link),aux);
		}
		Propagate(aux,GetConditionalLikelihoodVector(up->Out()),GetLength(up->GetBranch()));
		double logl = ComputeNodeLikelihood(up->Out(),0);
		loglmap[pair<Link*,Link*>(from,fromup)] = logl;
		GetTree()->Detach(down,up);
		// Link* tmp1 = GetTree()->Detach(down,up);
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveNonMPIGibbsSPRScan(link->Out(),trailer,down,up,loglmap);
		trailer = trailer->Next();
	}
}

void PhyloProcess::Create(Tree* intree, SequenceAlignment* indata,int indim)	{

	if (! data)	{
		data = indata;
		// MPI : master and slaves
		RateProcess::Create(data->GetNsite());
		ProfileProcess::Create(data->GetNsite(),indim);
		BranchProcess::Create(intree);

		empfreq = new double[data->GetNstate()];
		data->GetEmpiricalFreq(empfreq);

		loglarray = new double[GetNbranch()];
		// MPI : slaves only
		// for each slave, should specify the range of sites (sitemin <= i < sitemax)
		// SubstitutionProcess::Create(data->GetNsite(),indim, sitemin, sitemax);
		if (myid > 0) {
			int sitemin,sitemax,width = data->GetNsite()/(nprocs-1);
			sitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				sitemax = data->GetNsite();
			}
			else {
				sitemax = myid*width;
			} 
			SubstitutionProcess::Create(data->GetNsite(),indim,sitemin,sitemax);

			submap = new BranchSitePath**[GetNbranch()];
			for (int j=0; j<GetNbranch(); j++)	{
				submap[j] = 0;
			}
			nodestate = new int*[GetNnode()];
			condlmap = new double***[GetNlink()];
			CreateNodeStates();
			CreateMissingMap();
			CreateMappings();
			condflag = false;
		}
		else {
			sitemin = -1;
			sitemax = -1;
			nodestate = new int*[GetNnode()];
			CreateNodeStates();
		}
	}

	if (mintotweight == -1)	{
		mintotweight = ((double) GetDim()) / 4;
	}
}

void PhyloProcess::CreateMissingMap()	{

	missingmap = new int*[GetNbranch()];
	for (int j=0; j<GetNnode(); j++)	{
		missingmap[j] = new int[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			missingmap[j][i] = -1;
		}
	}
}

void PhyloProcess::DeleteMissingMap()	{

	for (int j=0; j<GetNbranch(); j++)	{
		delete[] missingmap[j];
	}
	delete[] missingmap;
}

void PhyloProcess::FillMissingMap()	{
	BackwardFillMissingMap(GetRoot());
	ForwardFillMissingMap(GetRoot(),GetRoot());
}

void PhyloProcess::BackwardFillMissingMap(const Link* from)	{

	int index = GetBranchIndex(from->GetBranch());
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		missingmap[index][i] = 0;
	}
	if (from->isLeaf())	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			int state = GetData(from)[i];
			if (state != -1)	{
				missingmap[index][i] = 1;
			}
		}
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			BackwardFillMissingMap(link->Out());
			int j = GetBranchIndex(link->Out()->GetBranch());
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				if (missingmap[j][i])	{
					missingmap[index][i] ++;
				}
			}
		}
	}
}

void PhyloProcess::ForwardFillMissingMap(const Link* from, const Link* up)	{

	int index = GetBranchIndex(from->GetBranch());
	int upindex = GetBranchIndex(up->GetBranch());
	if (from->isRoot())	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (missingmap[index][i] <= 1)	{
				missingmap[index][i] = 0;
			}
			else	{
				missingmap[index][i] = 2;
			}
		}
	}
	else	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (missingmap[index][i] > 0)	{
				if (missingmap[upindex][i])	{
					missingmap[index][i] = 1;
				}
				else	{
					if (from->isLeaf() || (missingmap[index][i] > 1))	{
					// if (missingmap[index][i] > 1)	{
						missingmap[index][i] = 2;
					}
					else	{
						missingmap[index][i] = 0;
					}
				}
			}
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		ForwardFillMissingMap(link->Out(),from);
	}
}

void PhyloProcess::SetTestSiteMinAndMax()	{

	bksitemax = sitemax;
	if (myid > 0) {
		int width = testnsite/(nprocs-1);
		testsitemin = (myid-1)*width;
		testsitemax = 0;
		if (myid == (nprocs-1)) {
			testsitemax = testnsite;
		}
		else {
			testsitemax = myid*width;
		} 
		// sitemax = sitemin + (testsitemax - testsitemin);
	}
}


void PhyloProcess::Delete() {

	if (data)	{
		// MPI slaves only
		if (myid > 0) {
			DeleteConditionalLikelihoods();
			DeleteNodeStates();
			DeleteMappings();
			delete[] submap;
			delete[] nodestate;
			delete[] condlmap;
			SubstitutionProcess::Delete();
		}
		// MPI master and slaves
		BranchProcess::Delete();
		ProfileProcess::Delete();
		RateProcess::Delete();

		delete[] empfreq;
	}
}

void PhyloProcess::GlobalUnclamp()	{

	assert(myid == 0);
	MESSAGE signal = UNCLAMP;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	dataclamped = 0;
}

void PhyloProcess::GlobalRestoreData()	{

	assert(myid == 0);
	MESSAGE signal = RESTOREDATA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	dataclamped = 1;
}

void PhyloProcess::GlobalSetDataFromLeaves()	{

	assert(myid == 0);
	MESSAGE signal = SETDATA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

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

	MPI_Status stat;
	int* tmp = new int[maxwidth * GetNtaxa()];

	for(int i=0; i<GetNprocs()-1; ++i) {
		MPI_Recv(tmp,(smax[i] - smin[i]) * GetNtaxa(),MPI_INT,i+1,TAG1,MPI_COMM_WORLD,&stat);
		int k = 0;
		for (int l=0; l<GetNtaxa(); l++)	{
			for (int j=smin[i]; j<smax[i]; j++)	{
				data->SetState(l,j,tmp[k]);
				k++;
			}
		}
	}
	delete[] tmp;
}

void PhyloProcess::GlobalSetNodeStates()	{

	assert(myid == 0);
	MESSAGE signal = SETNODESTATES;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

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

	MPI_Status stat;
	int* tmp = new int[maxwidth * GetNnode()];

	for(int i=0; i<GetNprocs()-1; ++i) {
		MPI_Recv(tmp,(smax[i] - smin[i]) * GetNnode(),MPI_INT,i+1,TAG1,MPI_COMM_WORLD,&stat);
		int k = 0;
		for (int l=0; l<GetNnode(); l++)	{
			for (int j=smin[i]; j<smax[i]; j++)	{
				nodestate[l][j] = tmp[k];
				k++;
			}
		}
	}
	delete[] tmp;
}


/*
double PhyloProcess::GlobalGetMeanDiversity()	{

	assert(myid == 0);
	MESSAGE signal = GETDIV;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;

	double total = 0;
	for(int i=1; i<nprocs; ++i) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		total += tmp;
	}
	return total / GetNsite();
}

double PhyloProcess::GlobalGetMeanSquaredFreq()	{

	assert(myid == 0);
	MESSAGE signal = GETSQUAREDFREQ;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;

	double total = 0;
	for(int i=1; i<nprocs; ++i) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		total += tmp;
	}
	return total / GetNsite();
}

double PhyloProcess::GlobalGetMeanFreqVariance()	{

	assert(myid == 0);
	MESSAGE signal = GETFREQVAR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;

	double total1[GetDim()];
	double total2[GetDim()];
    for (int k=0; k<GetDim(); k++)  {
        total1[k] = total2[k] = 0;
    }

    double* tmp = new double[2*GetDim()];
	for(int i=1; i<nprocs; ++i) {
		MPI_Recv(tmp,2*GetDim(),MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
        for (int k=0; k<GetDim(); k++)  {
            total1[k] += tmp[k];
            total2[k] += tmp[k+GetDim()];
        }
	}
    delete[] tmp;

    double meanvar = 0;
    for (int k=0; k<GetDim(); k++)  {
        total1[k] /= GetNsite();
        total2[k] /= GetNsite();
        total2[k] -= total1[k]*total1[k];
        meanvar += total2[k];
    }
    meanvar /= GetDim();
	return meanvar;
}
*/

void PhyloProcess::GlobalSetRatePrior(int inrateprior)	{

	assert(myid == 0);
	rateprior = inrateprior;
	MESSAGE signal = SETRATEPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&rateprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveSetRatePrior()	{

	MPI_Bcast(&rateprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalSetProfilePrior(int inprofileprior)	{

	assert(myid == 0);
	profileprior = inprofileprior;
	MESSAGE signal = SETPROFILEPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&profileprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveSetProfilePrior()	{

	MPI_Bcast(&profileprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalSetRootPrior(int inrootprior)	{

	assert(myid == 0);
	rootprior = inrootprior;
	MESSAGE signal = SETROOTPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&rootprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveSetRootPrior()	{

	MPI_Bcast(&rootprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveRestoreData()	{
	GetData()->Restore();
	dataclamped = 1;
}

void PhyloProcess::SlaveUnclamp()	{
	GetData()->Unclamp();
	dataclamped = 0;
}

void PhyloProcess::SlaveSetDataFromLeaves()	{
	SetDataFromLeaves();

	// mpi send the array
	int* tmp = new int[(sitemax - sitemin) * GetNtaxa()];
	int k = 0;
	for (int i=0; i<GetNtaxa(); i++)	{
		for (int j=sitemin; j<sitemax; j++)	{
			tmp[k] = data->GetState(i,j);
			k++;
		}
	}
	MPI_Send(tmp,(sitemax-sitemin)*GetNtaxa(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

	delete[] tmp;

}

void PhyloProcess::SlaveSetNodeStates()	{

	// mpi send the array
	int* tmp = new int[(sitemax - sitemin) * GetNnode()];
	int k = 0;
	for (int i=0; i<GetNnode(); i++)	{
		for (int j=sitemin; j<sitemax; j++)	{
			tmp[k] = nodestate[i][j];
			k++;
		}
	}
	MPI_Send(tmp,(sitemax-sitemin)*GetNnode(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

	delete[] tmp;

}

/*
void PhyloProcess::SlaveGetMeanDiversity()	{

	double div = GetData()->GetTotalDiversity(sitemin,sitemax);
	MPI_Send(&div,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
}

void PhyloProcess::SlaveGetMeanSquaredFreq()	{

	double div = GetData()->GetTotalSquaredFreq(sitemin,sitemax);
	MPI_Send(&div,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
}

void PhyloProcess::SlaveGetMeanFreqVariance()	{

    double* m = new double[2*GetDim()];
    for (int k=0; k<2*GetDim(); k++)    {
        m[k] = 0;
    }
	GetData()->GetTotalFreqMoments(m,sitemin,sitemax);
	MPI_Send(m,2*GetDim(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] m;
}
*/

void PhyloProcess::GlobalSimulateForward()	{

	// MPI
	assert(myid == 0);
	MESSAGE signal = SIMULATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}


void PhyloProcess::GlobalUnfold()	{

	assert(myid == 0);
	DeleteSuffStat();
	GlobalUpdateParameters();

	MESSAGE signal = UNFOLD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	GlobalUpdateConditionalLikelihoods();
}

void PhyloProcess::GlobalCollapse()	{

	// MPI
	// call Collapse() on slaves only
	// as for the master: should take care of one or two flags
	// conflag = false;
	assert(myid == 0);
	MESSAGE signal = COLLAPSE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	CreateSuffStat();
}

double PhyloProcess::GlobalComputeNodeLikelihood(const Link* from, int auxindex)	{ 
	// MPI
	// send messages to slaves : message "compute likelihood", with 2 arguments: GetLinkIndex(from) and auxindex
	// slaves: upon receiving message with two arguments fromindex and auxindex
	// call ComputeNodeLikelihood(GetLink(fromindex),auxindex)
	// return the value
	assert(myid == 0);
	MESSAGE signal = LIKELIHOOD;
	MPI_Status stat;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int i,args[] = {GetLinkIndex(from),auxindex};
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
	// master : sums up all values sent by slaves
	// store this sum into member variable logL
	// and return it

	logL = 0.0;
	double sum;
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(&sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		logL += sum;
	}
	return logL;
}

void PhyloProcess::GlobalReset(const Link* link, bool condalloc)	{

	// MPI
	// send a Reset message with GetLinkIndex(link) as argument
	// slaves: upon receiving message
	// call the Reset function with link corresponding to index received as argument of the message
	assert(myid == 0);
	MESSAGE signal = RESET;
	int args[2];
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	args[0] = GetLinkIndex(link);
	args[1] = (condalloc) ? 1 : 0;
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
}


void PhyloProcess::GlobalMultiply(const Link* from, const Link* to, bool condalloc)	{

	// MPI
	// send a Multiply message with GetLinkIndex(from) and GetLinkIndex(to) as argument
	// slaves: upon receiving message
	// call the Multiply function with links corresponding to the two indices received as argument
	assert(myid == 0);
	MESSAGE signal = MULTIPLY;
	int args[3];
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	args[0] = GetLinkIndex(from);
	args[1] = GetLinkIndex(to);
	args[2] = (condalloc) ? 1 : 0;
	MPI_Bcast(args,3,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalMultiplyByStationaries(const Link* from, bool condalloc)	{

	// MPI
	assert(myid == 0);
	MESSAGE signal = SMULTIPLY;
	int args[2];
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	args[0] = GetLinkIndex(from);
	args[1] = (condalloc) ? 1 : 0;
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalInitialize(const Link* from, const Link* link, bool condalloc)	{

	// MPI
	assert(myid == 0);
	MESSAGE signal = INITIALIZE;
	int args[3];
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	args[0] = GetLinkIndex(from);
	args[1] = GetLinkIndex(link);
	args[2] = (condalloc) ? 1 : 0;
	MPI_Bcast(args,3,MPI_INT,0,MPI_COMM_WORLD);
}


void PhyloProcess::GlobalPropagate(const Link* from, const Link* to, double time, bool condalloc)	{

	// MPI
	assert(myid == 0);
	MESSAGE signal = PROPAGATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    int args[3];
	args[0] = GetLinkIndex(from);
	args[1] = GetLinkIndex(to);
	args[2] = (condalloc) ? 1 : 0;
	MPI_Bcast(args,3,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

double PhyloProcess::GlobalProposeMove(const Branch* branch, double tuning)	{

	// MPI
	// master and all slaves should all call MoveBranch(branch,m)
	// should send a message with arguments: GetBranchIndex(branch), m
	// slaves should interpret the message, and apply on branch with index received as message argument
	assert(myid == 0);
	MESSAGE signal = PROPOSE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	int index = branch->GetIndex();
	MPI_Bcast(&index,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&m,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MoveBranch(branch,m);
	return m;
}

void PhyloProcess::GlobalRestore(const Branch* branch)	{

	// MPI
	// master and all slaves should all call RestoreBranch(branch)
	assert(myid == 0);
	MESSAGE signal = RESTORE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int n = branch->GetIndex();
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	Restore(branch);
}

void PhyloProcess::GlobalUpdateConditionalLikelihoods()	{

	// MPI
	// just send Updateconlikelihood message to all slaves
	assert(myid == 0);
	MESSAGE signal = UPDATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	GlobalComputeNodeLikelihood(GetRoot(),0);
	// GlobalCheckLikelihood();
}

Link* PhyloProcess::GlobalDetach(Link* down, Link* up)	{

	// MPI
	// master and all slaved should call 
	// GetTree()->Detach(down,up,fromdown,fromup);
	// but message passing will again  use link to index, then index to link, translations.
	assert(myid == 0);
	MESSAGE signal = DETACH;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int args[] = {GetLinkIndex(down),GetLinkIndex(up)};
	// int args[] = {down->GetIndex(),up->GetIndex(),fromdown->GetIndex(),fromup->GetIndex()};
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
	return GetTree()->Detach(down,up);
}

void PhyloProcess::GlobalAttach(Link* down, Link* up, Link* fromdown, Link* fromup)	{

	// MPI
	// same thing as for detach
	assert(myid == 0);
	MESSAGE signal = ATTACH;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int args[] = {GetLinkIndex(down),GetLinkIndex(up),GetLinkIndex(fromdown),GetLinkIndex(fromup)};
	// int args[] = {down->GetIndex(),up->GetIndex(),fromdown->GetIndex(),fromup->GetIndex()};
	MPI_Bcast(args,4,MPI_INT,0,MPI_COMM_WORLD);
	GetTree()->Attach(down,up,fromdown,fromup);
}

void PhyloProcess::GlobalRootAtRandom()	{

	//Link* newroot = GetTree()->ChooseLinkAtRandom();
	// I have to do it this way because GetLink returns a const Link pointer, 
	// whereas RootAt requires a non-const Link pointer...
	assert(myid == 0);
	// except root
	int n = GetTree()->CountInternalNodes(GetRoot());
	int choose = (int) (n * rnd::GetRandom().Uniform());
	// because root is first in choose internal node
	// choose ++;
	// MPI
	// call slaves, send a reroot message with argument newroot
	MESSAGE signal = ROOT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&choose,1,MPI_INT,0,MPI_COMM_WORLD);

	Link* tmp = 0;
	Link* newroot = GetTree()->ChooseInternalNode(GetRoot(),tmp,choose);
	if (newroot->isLeaf())	{
		cerr << "error : root at leaf\n";
		exit(1);
	}
	GetTree()->RootAt(newroot);
	GlobalUpdateConditionalLikelihoods();	
	
}


void PhyloProcess::GlobalGibbsSPRScan(Link* down, Link* up, double* loglarray)  {
	assert(myid == 0);
	int i,j,args[2],nbranch = GetNbranch();
	MPI_Status stat;
	MESSAGE signal = SCAN;
	args[0] = GetLinkIndex(down);
	args[1] = GetLinkIndex(up);

	// MPI3 : send message : GibbsSPRScan(idown,iup);
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);

	//
	// gather all slaves'arrays
	// into loglarray
	// of size GetNbranch();
	// (actually shorter than that, but should be ok)
	double dvector[nbranch];
	for(i=0; i<nbranch; ++i) {
		loglarray[i] = 0.0;
	}
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,nbranch,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			loglarray[j] += dvector[j];
		}
	}
}


// MPI
// for slaves
// I guess there should be some general function here
// for parsing and executing all messages received from master
// this should more or less implement a switch
// scanning all the possibilities mentioned above


// MPI: slave execute fonction
// waits for messages
// and call SlaveExecute();
void PhyloProcess::WaitLoop()	{
	MESSAGE signal;
	do {
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		if (signal == KILL) break;
		SlaveExecute(signal);
	} while(true);
}

// MPI: slave execute fonction
// waits for messages
// and call SlaveExecute();

void PhyloProcess::SlaveExecute(MESSAGE signal)	{
	int n,arg[4];
    int branchindex;
    double time;
	bool tvalue;

	switch(signal) {
	case SETRATEPRIOR:
		SlaveSetRatePrior();
		break;
	case SETPROFILEPRIOR:
		SlaveSetProfilePrior();
		break;
	case SETROOTPRIOR:
		SlaveSetRootPrior();
		break;
	case ROOT:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveRoot(n);
		break;
	case LIKELIHOOD:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveLikelihood(arg[0],arg[1]);
		break;
	case SCAN:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveGibbsSPRScan(arg[0],arg[1]);
		break;
	case PROPOSE:
		MPI_Bcast(&branchindex,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		SlavePropose(branchindex,time);
		break;
	case RESTORE:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveRestore(n);
		break;
    case PREPARESTEPPING:
		SlavePrepareStepping();
		break;
    case SETSTEPPINGFRAC:
        SlaveSetSteppingFraction();
        break;
	case RESET:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[1] == 1) ? true : false;
		SlaveReset(arg[0],tvalue);
		break;
    case RESETALL:
        SlaveResetAllConditionalLikelihoods();
        break;
	case MULTIPLY:
		MPI_Bcast(arg,3,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[2] == 1) ? true : false;
		SlaveMultiply(arg[0],arg[1],tvalue);
		break;
	case SMULTIPLY:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[1] == 1) ? true : false;
		SlaveSMultiply(arg[0],tvalue);
		break;
	case INITIALIZE:
		MPI_Bcast(arg,3,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[2] == 1) ? true : false;
		SlaveInitialize(arg[0],arg[1],tvalue);
		break;
	case PROPAGATE:
		MPI_Bcast(arg,3,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		tvalue = (arg[2] == 1) ? true : false;
		SlavePropagate(arg[0], arg[1], tvalue, time);
		break;
	case ATTACH:
		MPI_Bcast(arg,4,MPI_INT,0,MPI_COMM_WORLD);
		SlaveAttach(arg[0],arg[1],arg[2],arg[3]);
		break;
	case DETACH:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveDetach(arg[0],arg[1]);
		break;
	case NNI:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveNNI(GetLinkForGibbs(arg[0]),arg[1]);
		break;
	case KNIT:
		MPI_Bcast(arg,1,MPI_INT,0,MPI_COMM_WORLD);
		GetLinkForGibbs(arg[0])->Knit();
		break;
	case BRANCHPROPAGATE:
		MPI_Bcast(arg,1,MPI_INT,0,MPI_COMM_WORLD);
		PropagateOverABranch(GetLinkForGibbs(arg[0]));
		break;
	case UNFOLD:
		Unfold();
		break;
	case COLLAPSE:
		Collapse();
		break;
	case 	UPDATE:
		UpdateConditionalLikelihoods();
		break;
	case UPDATE_SRATE:
		SlaveUpdateSiteRateSuffStat();
		break;
	case UPDATE_SPROFILE:
		SlaveUpdateSiteProfileSuffStat();
		break;
	case UPDATE_BLENGTH:
		SlaveUpdateBranchLengthSuffStat();
		break;
	case PARAMETER_DIFFUSION:
		SlaveUpdateParameters();
		break;
	case BCAST_TREE:
		SlaveBroadcastTree();
		break;
	case UNCLAMP:
		SlaveUnclamp();
		break;
	case RESTOREDATA:
		SlaveRestoreData();
		break;
	case SETDATA:
		SlaveSetDataFromLeaves();
		break;
	case SETNODESTATES:
		SlaveSetNodeStates();
		break;
    /*
	case GETDIV:
		SlaveGetMeanDiversity();
		break;
	case GETSQUAREDFREQ:
		SlaveGetMeanSquaredFreq();
		break;
	case GETFREQVAR:
		SlaveGetMeanFreqVariance();
		break;
    */
	case CVSCORE:
		SlaveComputeCVScore();
		break;
	case SITELOGL:
		SlaveComputeSiteLogL();
		break;
	case STEPPINGSITELOGL:
		SlaveGetSiteSteppingLogLikelihood();
		break;
	case STATEPOSTPROBS:
		SlaveComputeStatePostProbs();
		break;
	case SITERATE:
		SlaveSendMeanSiteRate();
		break;
	case SETTESTDATA:
		SlaveSetTestData();
		break;
	case WRITE_MAPPING:
		SlaveWriteMappings();
		break;
	case COUNTMAPPING:
		SlaveCountMapping();
		break;
	case SIMULATE:
		SimulateForward();
		break;
    case EMPIRICALFRAC:
        SlaveSetEmpiricalFrac();
        break;
    case EMPIRICALPRIOR:
        SlaveSetEmpiricalPrior();
        break;
    /*
    case CREATESITE:
        SlaveCreateSiteDataStructures();
        break;
    case DELETESITE:
        SlaveDeleteSiteDataStructures();
        break;
    */
	
	default:
		// or : SubstitutionProcess::SlaveExecute?
		cerr << "slave could not process signal : " << signal << '\n';
		exit(1);
	}
}

void PhyloProcess::SlaveRoot(int n) {
	assert(myid > 0);
	Link* tmp = 0;
	Link* newroot = GetTree()->ChooseInternalNode(GetRoot(),tmp,n);
	GetTree()->RootAt(newroot);
}

void PhyloProcess::SlaveLikelihood(int fromindex,int auxindex) {
	assert(myid > 0);
	double lvalue = ComputeNodeLikelihood(GetLinkForGibbs(fromindex),auxindex);
	MPI_Send(&lvalue,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveGibbsSPRScan(int idown, int iup)	{
	assert(myid > 0);

	int n = 0;
	Link* down = GetLink(idown);
	Link* up = GetLink(iup);
	RecursiveGibbsSPRScan(GetRoot(),GetRoot(),down,up,loglarray,n);

	// MPI3 : send loglarray
	MPI_Send(loglarray,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::SlavePropose(int n,double x) {
	assert(myid > 0);
	const Branch* br = GetBranch(n);
	MoveBranch(br,x);	
}

void PhyloProcess::SlaveRestore(int n) {
	assert(myid > 0);
	const Branch* br = GetBranch(n);
	Restore(br);
}

void PhyloProcess::SlaveReset(int n,bool v) {
	assert(myid > 0);
	// const Link* link = GetLink(n);
	Reset(condlmap[n],v);	
}

void PhyloProcess::SlaveMultiply(int n,int m,bool v) {
	assert(myid > 0);
	// const Link* from = GetLink(n);
	// const Link* to = GetLink(m);
	// Multiply(condlmap[from],condlmap[to],v);
	Multiply(condlmap[n],condlmap[m],v);
}

void PhyloProcess::SlaveSMultiply(int n,bool v) {
	assert(myid > 0);
	// const Link* from = GetLink(n);
	// MultiplyByStationaries(condlmap[from],v);
	MultiplyByStationaries(condlmap[n],v);
}

void PhyloProcess::SlaveInitialize(int n,int m,bool v) {
	assert(myid > 0);
	// const Link* from = GetLink(n);
	const Link* link = GetLink(m);
	// Initialize(condlmap[from],GetData(link),v);
	Initialize(condlmap[n],GetData(link),v);
	// Initialize(condlmap[n],GetData(m),v);
}

void PhyloProcess::SlavePropagate(int n,int m,bool v,double t) {
	assert(myid > 0);
	// const Link* from = GetLink(n);
	// const Link* to = GetLink(m);
	// Propagate(condlmap[from],condlmap[to],t,v);
	Propagate(condlmap[n],condlmap[m],t,v);
	Offset(condlmap[m]);
}

void PhyloProcess::SlaveDetach(int n,int m) {
	assert(myid > 0);
	Link* down = GetLinkForGibbs(n);
	Link* up = GetLinkForGibbs(m);
	GetTree()->Detach(down,up);
}

void PhyloProcess::SlaveAttach(int n,int m,int p,int q) {
	assert(myid > 0);
	Link* down = GetLinkForGibbs(n);
	Link* up = GetLinkForGibbs(m);
	Link* fromdown = GetLinkForGibbs(p);
	Link* fromup = GetLinkForGibbs(q);
	GetTree()->Attach(down,up,fromdown,fromup);
}


//-----------------------------PhyloProcess--------------------------------------------
//-------------------------------------------------------------------------
//	* MPI SuffStat 
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PhyloProcess::GlobalUpdateBranchLengthSuffStat()	{

	assert(myid == 0);
	int i,j,nbranch = GetNbranch();
	MPI_Status stat;
	MESSAGE signal = UPDATE_BLENGTH;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<nbranch; ++i) {
		branchlengthsuffstatcount[i] = 0;
		branchlengthsuffstatbeta[i] = 0.0;
	}

	int ivector[nbranch];
	double dvector[nbranch];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,nbranch,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchlengthsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchlengthsuffstatbeta[j] += dvector[j];
		}
	}

	if (branchlengthsuffstatcount[0])	{
		cerr << "error at root\n";
		cerr << branchlengthsuffstatcount[0] << '\n';
	}
	if (branchlengthsuffstatbeta[0])	{
		cerr << "error at root\n";
		cerr << branchlengthsuffstatbeta[0] << '\n';
	}
	// check for nan
	for(int j=0; j<GetNbranch(); ++j) {
		if (std::isnan(branchlengthsuffstatbeta[j]))	{
			cerr << "in PhyloProcess::GlobalUpdateBranchLengthSuffStat: nan\n";
			exit(1);
		}
	}
}

void PhyloProcess::SlaveUpdateBranchLengthSuffStat()	{

	UpdateBranchLengthSuffStat();
	if (branchlengthsuffstatcount[0])	{
		cerr << "error at root in slave " << GetMyid() << "\n";
		cerr << branchlengthsuffstatcount[0] << '\n';
	}
	if (branchlengthsuffstatbeta[0])	{
		cerr << "error at root in slave " << GetMyid() << "\n";
		cerr << branchlengthsuffstatbeta[0] << '\n';
	}
	MPI_Send(branchlengthsuffstatcount,GetNbranch(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(branchlengthsuffstatbeta,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalUpdateSiteRateSuffStat()	{

	MESSAGE signal = UPDATE_SRATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveUpdateSiteRateSuffStat()	{

	UpdateSiteRateSuffStat();
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (std::isnan(siteratesuffstatbeta[i]))	{
			cerr << "in PhyloProcess::GlobalUpdateSiteRateSuffStat: nan ratesuffstatbeta\n";
			exit(1);
		}
	}
}

void PhyloProcess::GlobalGetMeanSiteRate()	{

	if (! meansiterate)	{
		meansiterate = new double[GetNsite()];
	}

	assert(myid == 0);
	int i,width,smin[nprocs-1],smax[nprocs-1];
	MPI_Status stat;
	MESSAGE signal = SITERATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	width = GetNsite()/(nprocs-1);
	for(i=0; i<nprocs-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (nprocs-2)) smax[i] = GetNsite();
	}
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(meansiterate+smin[i-1],smax[i-1]-smin[i-1],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
	}
}

void PhyloProcess::SlaveSendMeanSiteRate()	{
	assert(myid > 0);
	MPI_Send(meansiterate+sitemin,sitemax-sitemin,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalBroadcastTree()	{

	// tree->RegisterWith(tree->GetTaxonSet());
	// SetNamesFromLengths();	
	ostringstream os;
	tree->ToStream(os);
	string s = os.str();
	unsigned int len = s.length();
	unsigned char* bvector = new unsigned char[len];
	for (unsigned int i=0; i<len; i++)	{
		bvector[i] = s[i];
	}
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	delete[] bvector;

}


void PhyloProcess::SlaveBroadcastTree()	{

	int len;
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	unsigned char* bvector = new unsigned char[len];
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	ostringstream os;
	for (int i=0; i<len; i++)	{
		os << bvector[i];
	}
	istringstream is(os.str());
	tree->ReadFromStream(is);
	// SetLengthsFromNames();
	delete[] bvector;
}


void PhyloProcess::ReadPB(int argc, char* argv[])	{

	cerr << "in PhyloProcess::ReadPB\n";
	exit(1);
}

void PhyloProcess::Read(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << '\n';
	cerr << "burnin : " << burnin << "\n";
	cerr << "every  : " << every << '\n'; 
	cerr << "until  : " << until << '\n';
	cerr << '\n';

	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	list<double> lengthlist;
	list<double> alphalist;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;
		double alpha = GetAlpha();
		alphalist.push_back(alpha);
		double length = GetRenormTotalLength();
		lengthlist.push_back(length);
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	cerr << '\n';
	cerr << "		   post. mean (95 % CI)\n";
	cerr << "tree length      : ";
	printCI(lengthlist, cerr);
	cerr << '\n';
	cerr << "alpha paarameter : ";
	printCI(alphalist, cerr);
	cerr << '\n';
	cerr << '\n';
}

void PhyloProcess::ReadSiteRates(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	double* meanrate = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meanrate[i] = 0;
	}

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		QuickUpdate();

		GlobalGetMeanSiteRate();

		// double length = GetRenormTotalLength();
		for (int i=0; i<GetNsite(); i++)	{
			// meansiterate[i] *= length;
			meanrate[i] += meansiterate[i];
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	ofstream os((name + ".meansiterates").c_str());
	for (int i=0; i<GetNsite(); i++)	{
		meanrate[i] /= samplesize;
		os << i << '\t' << meanrate[i] << '\n';
	}
	cerr << "posterior mean site rates in " << name << ".meansiterates\n";

	delete[] meanrate;

}

void PhyloProcess::AllPostPred(string name, int burnin, int every, int until, int inrateprior, int inprofileprior, int inrootprior)	{

	GlobalSetRatePrior(inrateprior);
	GlobalSetProfilePrior(inprofileprior);
	GlobalSetRootPrior(inrootprior);

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	double* obstaxstat = new double[GetNtaxa()];
	int nstat = 5;
	double obsarray[nstat];
	obsarray[0] = data->GetMeanDiversity();
	obsarray[1] = data->GetMeanSquaredFreq();
	obsarray[2] = data->GetMeanFreqVariance();
	obsarray[3] = GetObservedCompositionalHeterogeneity(obstaxstat,obsarray[4]);

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	// cerr << "number of points : " << (until - burnin)/every << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;
	double meanstatarray[nstat];
	double varstatarray[nstat];
	double ppstatarray[nstat];
	for (int k=0; k<nstat; k++)	{
		meanstatarray[k] = 0;
		varstatarray[k] = 0;
		ppstatarray[k] = 0;
	}

	double* meantaxstat = new double[GetNtaxa()];
	double* vartaxstat = new double[GetNtaxa()];
	double* pptaxstat = new double[GetNtaxa()];
	double* taxstat = new double[GetNtaxa()];
	for (int j=0; j<GetNtaxa(); j++)	{
		meantaxstat[j] = 0;
		vartaxstat[j] = 0;
		pptaxstat[j] = 0;
	}
	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;

		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalUnclamp();
		GlobalCollapse();
		GlobalSetDataFromLeaves();

		double statarray[nstat];
		statarray[0] = data->GetMeanDiversity();
		statarray[1] = data->GetMeanSquaredFreq();
		statarray[2] = data->GetMeanFreqVariance();
		statarray[3] = GetCompositionalHeterogeneity(taxstat,statarray[4]);
		for (int j=0; j<GetNtaxa(); j++)	{
			meantaxstat[j] += taxstat[j];
			vartaxstat[j] += taxstat[j] * taxstat[j];
			if (taxstat[j] > obstaxstat[j])	{
				pptaxstat[j] ++;
			}
		}

		for (int k=0; k<nstat; k++)	{
			meanstatarray[k] += statarray[k];
			varstatarray[k] += statarray[k] * statarray[k];
			if (statarray[k] > obsarray[k])	{
				ppstatarray[k]++;
			}
		}

		GlobalRestoreData();
		GlobalUnfold();

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	cerr << '\n';
	for (int k=0; k<nstat; k++)	{
		meanstatarray[k] /= samplesize;
		varstatarray[k] /= samplesize;
		varstatarray[k] -= meanstatarray[k] * meanstatarray[k];
		ppstatarray[k] /= samplesize;
	}

	ofstream os((name + ".ppred").c_str());

	os << "diversity test\n";
	os << "obs div : " << obsarray[0] << '\n';
	os << "mean div: " << meanstatarray[0] << " +/- " << sqrt(varstatarray[0]) << '\n';
	os << "z-score : " << (meanstatarray[0] - obsarray[0]) / sqrt(varstatarray[0]) << '\n';
	os << "pp      : " << 1-ppstatarray[0] << '\n';

	os << '\n';

	os << "empirical convergence probability test\n";
	os << "obs     : " << obsarray[1] << '\n';
	os << "mean    : " << meanstatarray[1] << " +/- " << sqrt(varstatarray[1]) << '\n';
	os << "z-score : " << (obsarray[1] - meanstatarray[1]) / sqrt(varstatarray[1]) << '\n';
	os << "pp      : " << ppstatarray[1] << '\n';

	os << '\n';

	os << "across-site compositional heterogeneity test\n";
	os << "obs     : " << obsarray[2] << '\n';
	os << "mean    : " << meanstatarray[2] << " +/- " << sqrt(varstatarray[2]) << '\n';
	os << "z-score : " << (obsarray[2] - meanstatarray[2]) / sqrt(varstatarray[2]) << '\n';
	os << "pp      : " << ppstatarray[2] << '\n';
	cerr << "result of test in " << name << ".sitecomp\n";

	os << '\n';

	os << "compositional homogeneity test\n";
	os << '\n';
	os << "max heterogeneity across taxa\n";
	os << "obs comp : " << obsarray[3] << '\n';
	os << "mean comp: " << meanstatarray[3] << " +/- " << sqrt(varstatarray[3]) << '\n';
	os << "z-score : " << (obsarray[3] - meanstatarray[3]) / sqrt(varstatarray[3]) << '\n';
	os << "pp      : " << ppstatarray[3] << '\n';

	os << '\n';

	os << "mean squared heterogeneity across taxa\n";
	os << "obs comp : " << obsarray[4] << '\n';
	os << "mean comp: " << meanstatarray[4] << " +/- " << sqrt(varstatarray[4]) << '\n';
	os << "z-score : " << (obsarray[4] - meanstatarray[4]) / sqrt(varstatarray[4]) << '\n';
	os << "pp      : " << ppstatarray[4] << '\n';

	os << '\n';

	os << "taxonname\tobs\tmean pred\tz-score\tpp\n";
	for (int j=0; j<GetNtaxa(); j++)	{
		meantaxstat[j] /= samplesize;
		vartaxstat[j] /= samplesize;
		pptaxstat[j] /= samplesize;
		vartaxstat[j] -= meantaxstat[j] * meantaxstat[j];
		os << GetTaxonSet()->GetTaxon(j) << '\t' << obstaxstat[j] << '\t' << meantaxstat[j] << '\t' << (obstaxstat[j] - meantaxstat[j])/sqrt(vartaxstat[j]) << '\t' << pptaxstat[j] << '\n';
	}

	cerr << "results of all posterior predictive tests in " << name << ".ppred\n";
	cerr << '\n';
}

void PhyloProcess::PostPred(int ppredtype, string name, int burnin, int every, int until, int inrateprior, int inprofileprior, int inrootprior, int savetrees)	{

	GlobalSetRatePrior(inrateprior);
	GlobalSetProfilePrior(inprofileprior);
	GlobalSetRootPrior(inrootprior);

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	double* obstaxstat = new double[GetNtaxa()];
	double obs = 0;
	double obs2 = 0;
	if (ppredtype == 2)	{
		obs = data->GetMeanDiversity();
		// obs = GlobalGetMeanDiversity();
	}
	else if (ppredtype == 3)	{
		obs = GetObservedCompositionalHeterogeneity(obstaxstat,obs2);
	}
    else if (ppredtype == 4)    {
        obs = data->GetMeanSquaredFreq();
    }
    else if (ppredtype == 5)    {
        obs = data->GetMeanFreqVariance();
    }
	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	// cerr << "number of points : " << (until - burnin)/every << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;
	double meanstat = 0;
	double varstat = 0;
	double ppstat = 0;
	double meanstat2 = 0;
	double varstat2 = 0;
	double ppstat2 = 0;
	double* meantaxstat = new double[GetNtaxa()];
	double* vartaxstat = new double[GetNtaxa()];
	double* pptaxstat = new double[GetNtaxa()];
	double* taxstat = new double[GetNtaxa()];
	for (int j=0; j<GetNtaxa(); j++)	{
		meantaxstat[j] = 0;
		vartaxstat[j] = 0;
		pptaxstat[j] = 0;
	}
	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;

		if (savetrees)	{
			// output tree
			ostringstream s;
			s << name << "_ppred" << samplesize << ".tree";
			ofstream os(s.str().c_str());
			SetNamesFromLengths();
			RenormalizeBranchLengths();
			GetTree()->ToStream(os);
			DenormalizeBranchLengths();
			os.close();
		}

		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalUnclamp();
		GlobalCollapse();
		GlobalSetDataFromLeaves();

		if (ppredtype > 1)	{
			double stat = 0;
			double stat2 = 0;
			if (ppredtype == 2)	{
				stat = data->GetMeanDiversity();
			}
			else if (ppredtype == 3)	{
				stat = GetCompositionalHeterogeneity(taxstat,stat2);
				for (int j=0; j<GetNtaxa(); j++)	{
					meantaxstat[j] += taxstat[j];
					vartaxstat[j] += taxstat[j] * taxstat[j];
					if (taxstat[j] > obstaxstat[j])	{
						pptaxstat[j] ++;
					}
				}
			}
		    else if (ppredtype == 4)    {
			stat = data->GetMeanSquaredFreq();
		    }
		    else if (ppredtype == 5)    {
			stat = data->GetMeanFreqVariance();
		    }
			meanstat += stat;
			varstat += stat * stat;
			if (stat < obs)	{
				ppstat++;
			}
			meanstat2 += stat2;
			varstat2 += stat2 * stat2;
			if (stat2 < obs2)	{
				ppstat2++;
			}
		}
		else	{
			// write datafile
			ostringstream s;
			s << name << "_ppred" << samplesize << ".ali";
			ofstream os(s.str().c_str());
			data->ToStream(os);
			os.close();
		}

		GlobalRestoreData();
		GlobalUnfold();

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	cerr << '\n';
	if (ppredtype > 1)	{
		meanstat /= samplesize;
		varstat /= samplesize;
		varstat -= meanstat * meanstat;
		ppstat /= samplesize;

		meanstat2 /= samplesize;
		varstat2 /= samplesize;
		varstat2 -= meanstat2 * meanstat2;
		ppstat2 /= samplesize;
	}

	if (ppredtype == 1)	{
		cerr << "datasets in " << name << "_ppred<rep>.ali\n";
	}
	if (ppredtype == 2)	{
		ofstream os((name + ".div").c_str());
		os << "diversity test\n";
		os << "obs div : " << obs << '\n';
		os << "mean div: " << meanstat << " +/- " << sqrt(varstat) << '\n';
		os << "z-score : " << (meanstat - obs) / sqrt(varstat) << '\n';
		os << "pp      : " << ppstat << '\n';
		cerr << "result of diversity test in " << name << ".div\n";
	}
	else if (ppredtype == 3)	{
		ofstream os((name + ".comp").c_str());
		os << "compositional homogeneity test\n";

		os << '\n';
		os << "max heterogeneity across taxa\n";
		os << "obs comp : " << obs << '\n';
		os << "mean comp: " << meanstat << " +/- " << sqrt(varstat) << '\n';
		os << "z-score : " << (obs - meanstat) / sqrt(varstat) << '\n';
		os << "pp      : " << (1 - ppstat) << '\n';

		os << '\n';
		os << "mean squared heterogeneity across taxa\n";
		os << "obs comp : " << obs2 << '\n';
		os << "mean comp: " << meanstat2 << " +/- " << sqrt(varstat2) << '\n';
		os << "z-score : " << (obs2 - meanstat2) / sqrt(varstat2) << '\n';
		os << "pp      : " << (1 - ppstat2) << '\n';

		os << '\n';
		os << "taxonname\tobs\tmean pred\tz-score\tpp\n";
		for (int j=0; j<GetNtaxa(); j++)	{
			meantaxstat[j] /= samplesize;
			vartaxstat[j] /= samplesize;
			pptaxstat[j] /= samplesize;
			vartaxstat[j] -= meantaxstat[j] * meantaxstat[j];
			os << GetTaxonSet()->GetTaxon(j) << '\t' << obstaxstat[j] << '\t' << meantaxstat[j] << '\t' << (obstaxstat[j] - meantaxstat[j])/sqrt(vartaxstat[j]) << '\t' << pptaxstat[j] << '\n';
		}
		cerr << "result of compositional homogeneity test in " << name << ".comp\n";
	}
    else if (ppredtype == 4)	{
		ofstream os((name + ".siteconvprob").c_str());
		os << "empirical convergence probability test\n";
		os << "obs     : " << obs << '\n';
		os << "mean    : " << meanstat << " +/- " << sqrt(varstat) << '\n';
		os << "z-score : " << (obs - meanstat) / sqrt(varstat) << '\n';
		os << "pp      : " << 1-ppstat << '\n';
		cerr << "result of test in " << name << ".siteconvprob\n";
	}
    else if (ppredtype == 5)	{
		ofstream os((name + ".sitecomp").c_str());
		os << "across-site compositional heterogeneity test\n";
		os << "obs     : " << obs << '\n';
		os << "mean    : " << meanstat << " +/- " << sqrt(varstat) << '\n';
		os << "z-score : " << (obs - meanstat) / sqrt(varstat) << '\n';
		os << "pp      : " << 1-ppstat << '\n';
		cerr << "result of test in " << name << ".sitecomp\n";
	}
	cerr << '\n';
}

void PhyloProcess::GlobalSetTestData()	{
	testnsite = testdata->GetNsite();
    if (testnsite > GetNsite()) {
        cerr << "error in CV: validation dataset should have less columns than training dataset\n";
        exit(1);
    }
	int* tmp = new int[testnsite * GetNtaxa()];
	testdata->GetDataVector(tmp);

	MESSAGE signal = SETTESTDATA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);

	delete[] tmp;
}

void PhyloProcess::SlaveSetTestData()	{

	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	int* tmp = new int[testnsite * GetNtaxa()];
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	
	SetTestSiteMinAndMax();
	data->SetTestData(testnsite,sitemin,testsitemin,testsitemax,tmp);

	delete[] tmp;
}

void PhyloProcess::ReadCV(string testdatafile, string name, int burnin, int every, int until, int iscodon, GeneticCodeType codetype)	{
	
	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	if (iscodon)	{
		SequenceAlignment* tempdata = new FileSequenceAlignment(testdatafile,0,myid);
		testdata = new CodonSequenceAlignment(tempdata,true,codetype);
	}
	else	{
		testdata = new FileSequenceAlignment(testdatafile,0,myid);
	}
	GlobalSetTestData();

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;
	vector<double> scorelist;

	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		MPI_Status stat;
		MESSAGE signal = CVSCORE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		double tmp = 0;
		double score = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			score += tmp;
		}
		scorelist.push_back(score);
		// cerr << score << '\n';
		
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}

	cerr << '\n';

	double max = 0;
	for (int j=0; j<samplesize; j++)	{
		if ((!j) || (max < scorelist[j]))	{
			max = scorelist[j];
		}
	}

	double tot = 0;
	for (int j=0; j<samplesize; j++)	{
		tot += exp(scorelist[j] - max);
	}
	tot /= samplesize;
	
	double meanscore = log(tot) + max;
	
	ofstream os((name + ".cv").c_str());
	os << meanscore << '\n';
	cerr << meanscore << '\n';
}

void PhyloProcess::ReadSiteCV(string testdatafile, string name, int burnin, int every, int until, int iscodon, GeneticCodeType codetype)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	if (iscodon)	{
		SequenceAlignment* tempdata = new FileSequenceAlignment(testdatafile,0,myid);
		testdata = new CodonSequenceAlignment(tempdata,true,codetype);
	}
	else	{
		testdata = new FileSequenceAlignment(testdatafile,0,myid);
	}
	GlobalSetTestData();

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

    int testwidth = testnsite/(nprocs-1);
	int testsmin[GetNprocs()-1];
	int testsmax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		testsmin[i] = testwidth*i;
		testsmax[i] = testwidth*(i+1);
		if (i == (GetNprocs()-2)) testsmax[i] = testnsite;
	}

    int width = GetNsite()/(nprocs-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	int maxwidth = 0;
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = smin[i] + testsmax[i] - testsmin[i];
		if (maxwidth < (smax[i] - smin[i]))	{
			maxwidth = smax[i] - smin[i];
		}
	}

	double* tmp = new double[GetNsite()];
	vector<double>* logl = new vector<double>[testnsite];

	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		MPI_Status stat;
		MESSAGE signal = SITELOGL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        int count = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			for (int j=smin[i-1]; j<smax[i-1]; j++)	{
				if (std::isnan(tmp[j]))	{
					cerr << "error: nan logl received by master\n";
					cerr << "site : " << j << '\n';
					cerr << "proc : " << i << '\n';
					exit(1);
				}
				logl[count].push_back(tmp[j]);
                count++;
			}
		}
        if (count != testnsite) {
            cerr << "error in read site cv: non matching number of sites (testnsite)\n";
            exit(1);
        }
		
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
    cerr << '\n';

    /*
    ofstream sos((name + ".cvsitelogl").c_str());
    for (int j=0; j<samplesize; j++)    {
        for (int i=0; i<testnsite; i++) {
            sos << logl[i][j] << '\t';
        }
        sos << '\n';
    }

    cerr << "site log likelihoods over the mcmc in " << name << ".cvsitelogl\n";
    */

    ofstream los((name + ".cvlog").c_str());
    los << samplesize << '\t' << testnsite << '\n';

    double* cvscore = new double[testnsite];
    double meancvscore = 0;

	double* siteess = new double[testnsite];
    double meaness = 0;
    double miness = 10;
    double nminess = 0;

    for (int i=0; i<testnsite; i++) {
        double max = 0;
        for (int j=0; j<samplesize; j++)	{
            if ((!j) || (max < logl[i][j]))	{
                max = logl[i][j];
            }
        }

        double tot = 0;
        for (int j=0; j<samplesize; j++)	{
            tot += exp(logl[i][j] - max);
        }

        double invess = 0;
        for (int j=0; j<samplesize; j++)	{
            double w = exp(logl[i][j] - max) / tot;
            invess += w*w;
        }
        siteess[i] = 1.0 / invess;
        if (siteess[i] < miness)   {
            nminess++;
        }
        meaness += siteess[i];

        tot /= samplesize;
        
        cvscore[i] = log(tot) + max;
        meancvscore += cvscore[i];
    }
    meancvscore /= testnsite;
    meaness /= testnsite;
    nminess /= testnsite;
        
	ofstream scos((name + ".sitecv").c_str());
    scos << "site\tcv\tess\n";
	for (int i=0; i<testnsite; i++) {
		scos << i+1 << '\t' << cvscore[i] << '\t' << siteess[i] << '\n';
	}

    // joint cv
	double max = 0;
    double* scorelist = new double[samplesize];

	for (int j=0; j<samplesize; j++)	{
        double tot = 0;
        for (int i=0; i<testnsite; i++) {
            tot += logl[i][j];
        }
        scorelist[j] = tot;
		if ((!j) || (max < scorelist[j]))	{
			max = scorelist[j];
		}
	}

    double meanlog = 0;
    double varlog = 0;
    for (int j=0; j<samplesize; j++)    {
        meanlog += scorelist[j];
        varlog += scorelist[j]*scorelist[j];
    }
    meanlog /= samplesize;
    varlog /= samplesize;
    varlog -= meanlog*meanlog;

	double tot = 0;
	for (int j=0; j<samplesize; j++)	{
		tot += exp(scorelist[j] - max);
	}

    double invess = 0;
	for (int j=0; j<samplesize; j++)	{
        double w = exp(scorelist[j] - max) / tot;
        invess += w*w;
	}
	tot /= samplesize;
	
	ofstream cos((name + ".cv").c_str());
	cos << "site cv (posterior average taken separately for each site of the validation set)\n";
    cos << "mcmc estimate  " << testnsite * meancvscore << '\n';
    cos << "per site       " << meancvscore << '\n';
    cos << "ESS            " << meaness << '\n';
    cos << "%(ESS<10)      " << 100*nminess << '\n';
    cerr << "cv results in " << name << ".cv\n";
}

void PhyloProcess::ReadSiteLogL(string name, int burnin, int every, int until, int verbose)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

    // posterior mean and variance (across the chain) of site-specific logls
    vector<double> site_postmeanlogl(GetNsite(),0);
    vector<double> site_postvarlogl(GetNsite(),0);

	double* tmp = new double[GetNsite()];
	vector<double>* logl = new vector<double>[GetNsite()];

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

	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		MPI_Status stat;
		MESSAGE signal = SITELOGL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		double total = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			for (int j=smin[i-1]; j<smax[i-1]; j++)	{
				if (std::isnan(tmp[j]))	{
					cerr << "error: nan logl received by master\n";
					cerr << "site : " << j << '\n';
					cerr << "proc : " << i << '\n';
					exit(1);
				}
				logl[j].push_back(tmp[j]);
				site_postmeanlogl[j] += tmp[j];
                site_postvarlogl[j] += tmp[j]*tmp[j];
				total += tmp[j];
			}
		}
		
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
    cerr << '\n';

    if (verbose)    {
        ofstream sos((name + ".mcmcsitelogls").c_str());
        for (int j=0; j<samplesize; j++)    {
            for (int i=0; i<GetNsite(); i++) {
                sos << logl[i][j] << '\t';
            }
            sos << '\n';
        }

        cerr << "site log likelihoods over the mcmc in " << name << ".mcmcsitelogls\n";
    }

    double mean_postmeanlogl = 0;
    double var_postmeanlogl = 0;
    double mean_postvarlogl = 0;
    double var_postvarlogl = 0;

    for (int i=0; i<GetNsite(); i++)    {
        site_postmeanlogl[i] /= samplesize;
        site_postvarlogl[i] /= samplesize;
        site_postvarlogl[i] -= site_postmeanlogl[i] * site_postmeanlogl[i];
        mean_postmeanlogl += site_postmeanlogl[i];
        var_postmeanlogl += site_postmeanlogl[i] * site_postmeanlogl[i];
        mean_postvarlogl += site_postvarlogl[i];
        var_postvarlogl += site_postvarlogl[i] * site_postvarlogl[i];
    }

    mean_postmeanlogl /= GetNsite();
    var_postmeanlogl /= GetNsite();
    var_postmeanlogl -= mean_postmeanlogl * mean_postmeanlogl;

    mean_postvarlogl /= GetNsite();
    var_postvarlogl /= GetNsite();
    var_postvarlogl -= mean_postvarlogl * mean_postvarlogl;

    // site cpo score (in log)
	double* site_logcpo = new double[GetNsite()];
    // mean and var across sites
	double mean_logcpo = 0;
	double var_logcpo = 0;
    // corresponding ESS
	double* cpo_siteess = new double[GetNsite()];
    // mean ESS
    double cpo_meaness = 0;
    // frac ESS<10
    double cpo_ness10 = 0;

    // log of site posterior mean likelihoods (for waic1)
    double* site_logpostmeanl = new double[GetNsite()];
    // mean and var across sites
    double mean_logpostmeanl= 0;
    double var_logpostmeanl = 0;
    // corresponding ESS
    double* postmeanl_siteess = new double[GetNsite()];
    // mean ESS
    double postmeanl_meaness = 0;
    // frac ESS < 10
    double postmeanl_ness10 = 0;

	for (int i=0; i<GetNsite(); i++)	{
		double min = 0;
        double max = 0;
		int count = 0;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			if ((!count) || (min > (*j)))	{
				min = *j;
			}
			if ((!count) || (max < (*j)))	{
				max = *j;
			}
			count++;
		}

        // cpo
        double site_cpo = 0;
        vector<double> cpo_weight;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			site_cpo += exp(min - (*j));
            cpo_weight.push_back(exp(min-(*j)));
		}

        double cpo_invess = 0;
        for (int j=0; j<count; j++) {
            cpo_weight[j] /= site_cpo;
            cpo_invess += cpo_weight[j] * cpo_weight[j];
        }
		site_cpo /= count;
		site_logcpo[i] = min - log(site_cpo);
		mean_logcpo += site_logcpo[i];
		var_logcpo += site_logcpo[i] * site_logcpo[i];
        cpo_siteess[i] = 1.0 / cpo_invess;
        cpo_meaness += cpo_siteess[i];
        if (cpo_siteess[i] < 10.0)  {
            cpo_ness10 ++;
        }

        // postmeanl
        double site_postmeanl = 0;
        vector<double> postmeanl_weight;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			site_postmeanl += exp((*j) - max);
            postmeanl_weight.push_back(exp((*j) - max));
		}

        double postmeanl_invess = 0;
        for (int j=0; j<count; j++) {
            postmeanl_weight[j] /= site_postmeanl;
            postmeanl_invess += postmeanl_weight[j] * postmeanl_weight[j];
        }
		site_postmeanl /= count;
		site_logpostmeanl[i] = log(site_postmeanl) + max;
        mean_logpostmeanl += site_logpostmeanl[i];
        var_logpostmeanl += site_logpostmeanl[i] * site_logpostmeanl[i];
        
        postmeanl_siteess[i] = 1.0 / postmeanl_invess;
        postmeanl_meaness += postmeanl_siteess[i];
        if (postmeanl_siteess[i] < 10.0)  {
            postmeanl_ness10 ++;
        }
	}

	mean_logcpo /= GetNsite();
	var_logcpo /= GetNsite();
	var_logcpo -= mean_logcpo * mean_logcpo;
    cpo_meaness /= GetNsite();
    cpo_ness10 /= GetNsite();

	mean_logpostmeanl /= GetNsite();
	var_logpostmeanl /= GetNsite();
	var_logpostmeanl -= mean_logpostmeanl * mean_logpostmeanl;
    postmeanl_meaness /= GetNsite();
    postmeanl_ness10 /= GetNsite();

	ofstream os((name + ".sitelogl").c_str());
    os << "site\tlogl\tvar\tlogcpo\tess\tlogpostmeanl\tess\n";
	for (int i=0; i<GetNsite(); i++)	{
		os << i+1 << '\t' << site_postmeanlogl[i] << '\t' << site_postvarlogl[i] << '\t' << site_logcpo[i] << '\t' << cpo_siteess[i] << '\t' << site_logpostmeanl[i] << '\t' << postmeanl_siteess[i] << '\n';
	}

	ofstream cos((name + ".cpo").c_str());

    cos << "wAIC          " << mean_logpostmeanl - mean_postvarlogl << '\n';
    cos << "mean(ESS)      " << postmeanl_meaness << '\n';
    cos << "%(ESS<10)      " << 100*postmeanl_ness10 << '\n';
    cos << '\n';
	cos << "LOO-CV        " << mean_logcpo << '\n';
    cos << "mean(ESS)      " << cpo_meaness << '\n';
    cos << "%(ESS<10)      " << 100*cpo_ness10 << '\n';
    cos << '\n';
    cos << "note: these are raw point estimates\n";
    cos << "to get an estimate of the bias / stdev, at least 2 independent runs are needed\n";
    cos << "these can then be post-analysed using the python script pbmpi/scripts/read_loocv_waic.py\n";
}

void PhyloProcess::ReadAncestral(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	double* allocmeanstatepostprob = new double[GetNsite()*GetNnode()*GetGlobalNstate()];
	double* allocstatepostprob = new double[GetNsite()*GetNnode()*GetGlobalNstate()];

	double*** meanstatepostprob = new double**[GetNsite()];
	double*** statepostprob = new double**[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meanstatepostprob[i] = new double*[GetNnode()];
		statepostprob[i] = new double*[GetNnode()];
		for (int j=0; j<GetNnode(); j++)	{
			meanstatepostprob[i][j] = allocmeanstatepostprob + (i*GetNnode() + j)*GetGlobalNstate();
			statepostprob[i][j] = allocstatepostprob + (i*GetNnode() + j)*GetGlobalNstate();
			for (int k=0; k<GetGlobalNstate(); k++)	{
				meanstatepostprob[i][j][k] = 0;
			}
		}
	}

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

	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		MPI_Status stat;
		MESSAGE signal = STATEPOSTPROBS;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		for(int proc=1; proc<GetNprocs(); proc++) {
			MPI_Recv(allocstatepostprob+smin[proc-1]*GetNnode()*GetGlobalNstate(),(smax[proc-1]-smin[proc-1])*GetNnode()*GetGlobalNstate(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
			for (int i=smin[proc-1]; i<smax[proc-1]; i++)	{
				for (int j=0; j<GetNnode(); j++)	{
					double tot = 0;
					for (int k=0; k<GetGlobalNstate(); k++)	{
						meanstatepostprob[i][j][k] += statepostprob[i][j][k];
						tot += statepostprob[i][j][k];
					}
					/*
					if (fabs(tot-1) > 1e-6)	{
						cerr << "error in receive post probs\n";
						cerr << tot << '\n';
						exit(1);
					}
					*/
				}
			}
		}
		
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}

	for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNnode(); j++)	{
			for (int k=0; k<GetGlobalNstate(); k++)	{
				meanstatepostprob[i][j][k] /= samplesize;
			}
		}
	}

	WriteStatePostProbs(meanstatepostprob,name,GetRoot());
	cerr << '\n';
	cerr << "ancestral state posterior probabilities in " << name << "_nodelabel_taxon1_taxon2_.ancstatepostprob\n";
	cerr << "for MRCA of taxon1 and taxon2\n";
	cerr << '\n';

	for (int i=0; i<GetNsite(); i++)	{
		delete[] meanstatepostprob[i];
		delete[] statepostprob[i];
	}
	delete[] meanstatepostprob;
	delete[] statepostprob;
	delete[] allocmeanstatepostprob;
	delete[] allocstatepostprob;
}

void PhyloProcess::SlaveComputeStatePostProbs()	{

	UpdateConditionalLikelihoods();

	// allocate arrays
	double* allocstatepostprob = new double[(GetSiteMax() - GetSiteMin())*GetNnode()*GetGlobalNstate()];
	double*** statepostprob = new double**[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		statepostprob[i] = new double*[GetNnode()];
		for (int j=0; j<GetNnode(); j++)	{
			statepostprob[i][j] = allocstatepostprob + ((i-GetSiteMin())*GetNnode() + j)*GetGlobalNstate();
		}
	}

	RecursiveComputeStatePostProbs(statepostprob,GetRoot(),0);
	/*
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		for (int j=0; j<GetNnode(); j++)	{
			double tot = 0;
			for (int k=0; k<GetGlobalNstate(); k++)	{
				tot += statepostprob[i][j][k];
			}
			if (fabs(tot-1) > 1e-6)	{
				cerr << "error in send post probs\n";
				cerr << tot << '\n';
				exit(1);
			}
		}
	}
	*/

	// send
	MPI_Send(allocstatepostprob,(GetSiteMax()-GetSiteMin())*GetNnode()*GetGlobalNstate(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	// delete
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		delete[] statepostprob[i];
	}
	delete[] statepostprob;
	delete[] allocstatepostprob;
}

void PhyloProcess::ComputeStatePostProbs(double*** statepostprob, const Link* from, int auxindex)	{

	if (! myid)	{
		cerr << "error : master doing slave's work\n";
		exit(1);
	}
	double*** aux = 0;
	bool localaux = false;
	if (auxindex != -1)	{
		aux = condlmap[auxindex];
	}
	else	{
		localaux = true;
		aux = CreateConditionalLikelihoodVector();
	}

	if (from->isLeaf())	{
		Initialize(aux,GetData(from));
	}
	else	{
		Reset(aux);
	}
	if (! from->isRoot())	{
		Multiply(GetConditionalLikelihoodVector(from),aux);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! link->isRoot())	{
			Multiply(GetConditionalLikelihoodVector(link),aux);
		}
	}
	MultiplyByStationaries(aux);
	int j = GetNodeIndex(from->GetNode());
	ConditionalLikelihoodsToStatePostProbs(aux,statepostprob,j);
	if (localaux)	{
		DeleteConditionalLikelihoodVector(aux);
	}
}

void PhyloProcess::RecursiveComputeStatePostProbs(double*** statepostprob, const Link* from, int auxindex)	{

	ComputeStatePostProbs(statepostprob,from,auxindex);
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// WARNING: preorder pruning does not update leaf condtional likelihood vectors (not necessary in the present case)
		// so the following will issue an error message if tried on leaf
		if (! link->Out()->isLeaf())	{
			RecursiveComputeStatePostProbs(statepostprob,link->Out(),auxindex);
		}
	}
}

void PhyloProcess::WriteStatePostProbs(double*** statepostprob, string name, const Link* from)	{

	ostringstream s;
	int nodelabel = GetNodeIndex(from->GetNode());
	s << name << "_" << nodelabel << "_" << GetLeftMost(from) << "_" << GetRightMost(from) << ".ancstatepostprob";
	ofstream os(s.str().c_str());

	os << GetNsite() << '\t' << GetGlobalNstate();
	for (int k=0; k<GetGlobalNstate(); k++)	{
		os << '\t' << GetStateSpace()->GetState(k);
	}
	os << '\n';
	for (int i=0; i<GetNsite(); i++)	{
		os << i+1;
		for (int k=0; k<GetGlobalNstate(); k++)	{
			os << '\t' << statepostprob[i][nodelabel][k];
		}
		os << '\n';
	}

	for (const Link* link=from->Next(); link!=from; link=link->Next()){
		if (! link->Out()->isLeaf())	{
			WriteStatePostProbs(statepostprob,name,link->Out());
		}
	}
}

double PhyloProcess::GlobalGetFullLogLikelihood()  {

    MPI_Status stat;
    MESSAGE signal = SITELOGL;
    MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

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

	double* tmp = new double[GetNsite()];
    double total = 0;
    for(int i=1; i<GetNprocs(); ++i) {
        MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        for (int j=smin[i-1]; j<smax[i-1]; j++)	{
            if (ActiveSite(j))  {
                if (std::isnan(tmp[j]))	{
                    cerr << "error: nan logl received by master\n";
                    cerr << "site : " << j << '\n';
                    cerr << "proc : " << i << '\n';
                    exit(1);
                }
                total += tmp[j];
            }
        }
    }
    delete[] tmp;
    return total;
}


void PhyloProcess::ReadMap(string name, int burnin, int every, int until){
  	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}
	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;
	double meandiff = 0;
	double vardiff = 0;
	double meanobs = 0;
	for(int i = 0; i < GetNsite(); i++){
		stringstream osfmap;
		osfmap << name << '_' << i << ".map";
		ofstream osmap((osfmap.str()).c_str());
		osmap.close();
	}
	while (i < until)	{
		cerr << ".";
		// cerr << i << '\t' << rnd::GetRandom().Uniform() << '\n';

		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		// prepare file for ancestral node states
		ostringstream s;
		s << name << "_" << samplesize << ".nodestates";
		ofstream sos(s.str().c_str());

		// quick update and mapping on the fly
		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalCollapse();

		// write posterior mappings
		GlobalWriteMappings(name);

		// write posterior ancestral node states
		GlobalSetNodeStates();
		WriteNodeStates(sos,GetRoot());
		sos << '\n';

		double obs = GlobalCountMapping();

		//Posterior Predictive Mappings
		GlobalUnfold();
		GlobalUnclamp();
		GlobalCollapse();

		GlobalSetDataFromLeaves();

		// write posterior predictive mappings
		GlobalWriteMappings(name);

		// write posterior predictive ancestral node states
		GlobalSetNodeStates();
		WriteNodeStates(sos,GetRoot());

		double pred = GlobalCountMapping();

		obs /= GetNsite();
		pred /= GetNsite();

		meandiff += obs - pred;
		vardiff += (obs-pred)*(obs-pred);
		meanobs += obs;

		GlobalRestoreData();
		GlobalUnfold();

		for(int i = 0; i < GetNsite(); i++){
			stringstream osfmap;
			osfmap << name << '_' << i << ".map";
			ofstream osmap((osfmap.str()).c_str(), ios_base::app);
			osmap << '\n';
			osmap.close();
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	meandiff /= samplesize;
	vardiff /= samplesize;
	vardiff -= meandiff*meandiff;
	meanobs /= samplesize;
	cerr << "mean obs : " << meanobs << '\n';
	cerr << meandiff << '\t' << sqrt(vardiff) << '\n';
}

void PhyloProcess::GlobalWriteMappings(string name){
	MESSAGE signal = WRITE_MAPPING;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	 //send the chain name
	ostringstream os;
	os << name;
	string s = os.str();
	unsigned int len = s.length();
	unsigned char* bvector = new unsigned char[len];
	for (unsigned int i=0; i<len; i++)	{
		bvector[i] = s[i];
	}
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	delete[] bvector;

}

void PhyloProcess::SlaveWriteMappings(){

	int len;
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	unsigned char* bvector = new unsigned char[len];
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	ostringstream os;
	for (int i=0; i<len; i++)	{
		os << bvector[i];
	}
	string name = os.str();
	delete[] bvector;

	for(int i = sitemin; i < sitemax; i++){
		stringstream osfmap;
		osfmap << name << '_' << i << ".map";
		ofstream osmap((osfmap.str()).c_str(), ios_base::app);
		WriteTreeMapping(osmap, GetRoot(), i);
		osmap.close();
	}
}


void PhyloProcess::WriteTreeMapping(ostream& os, const Link* from, int i){
	if(from->isLeaf()){
		os << from->GetNode()->GetName();
	}
	else{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next()){
			WriteTreeMapping(os, link->Out(), i);
			if (link->Next() != from)       {
				os << ',';
			}
		}
		os << ')';
	}
	if(from->isRoot()){
		BranchSitePath* mybsp = submap[GetBranchIndex(from->Next()->GetBranch())][i];
		os << '_' << GetStateSpace()->GetState(mybsp->Init()->GetState()) << ";\n";     
		/*
		BranchSitePath* mybsp = submap[0][i];
		os << '_' << GetStateSpace()->GetState(mybsp->Last()->GetState()) << ";\n";		
		*/
	}
	else{
		BranchSitePath* mybsp = submap[GetBranchIndex(from->GetBranch())][i];
		double l = GetLength(from->GetBranch());
		os << '_' << GetStateSpace()->GetState(mybsp->Last()->GetState());
		for(Plink* plink = mybsp->Last(); plink ; plink = plink->Prev()){
			os << ':' << plink->GetRelativeTime() * l << ':' << GetStateSpace()->GetState(plink->GetState());
		}
	}
}

void PhyloProcess::WriteNodeStates(ostream& os, const Link* from)	{
	os << GetLeftMost(from) << '\t' << GetRightMost(from) << '\t';
	int nodelabel = GetNodeIndex(from->GetNode());
	for (int i=0; i<GetNsite(); i++)	{
		os << GetStateSpace()->GetState(nodestate[nodelabel][i]);
	}
	os << '\n';
		
	for (const Link* link=from->Next(); link!=from; link=link->Next()){
		WriteNodeStates(os,link->Out());
	}
}

int PhyloProcess::CountMapping()	{

	int total = 0;	
	for(int i = sitemin; i < sitemax; i++){
		total += CountMapping(i);
	}
	return total;
}

int PhyloProcess::CountMapping(int i)	{
	return 0;
}

int PhyloProcess::GlobalCountMapping()	{

	assert(myid==0);
	MESSAGE signal = COUNTMAPPING;
	MPI_Status stat;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	int i, count, totalcount=0;
	for (i=1; i<nprocs; ++i)	{
		MPI_Recv(&count,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD, &stat);
		totalcount += count;
	}
	return totalcount;

}

void PhyloProcess::SlaveCountMapping()	{

	int count = CountMapping();
	MPI_Send(&count,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);

}


