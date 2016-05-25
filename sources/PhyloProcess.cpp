
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
extern MPI_Datatype Propagate_arg;

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
		siteratesuffstatcount = new int[GetNsite()];
		siteratesuffstatbeta = new double[GetNsite()];
	}
	if (! branchlengthsuffstatcount)	{
		branchlengthsuffstatcount = new int[GetNbranch()];
		branchlengthsuffstatbeta = new double[GetNbranch()];
	}
	activesuffstat = true;
}

void PhyloProcess::DeleteSuffStat()	{

	delete[] siteratesuffstatcount;
	delete[] siteratesuffstatbeta;
	siteratesuffstatcount = 0;
	siteratesuffstatbeta = 0;
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

	// not necessary, but
	// should be updated when entering this function
	// UpdateConditionalLikelihoods();
	// 
	// this is important
	// for the conditional likelihoods of the cut-then-regrafted subtree to be OK

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
	// UpdateConditionalLikelihoods();
	
	// double* loglarray = new double[GetNbranch()];
	GlobalGibbsSPRScan(down,up,loglarray);
	map<pair<Link*,Link*>, double> loglmap;
	int n = 0;
	RecursiveGibbsFillMap(GetRoot(),GetRoot(),loglmap,loglarray,n);

	double max = 0;
	int j = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (isinf(i->second))	{
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
		if (isinf(total))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total))	{
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
	// GetTree()->Attach(down,up,i->first.first,i->first.second);
	GlobalUpdateConditionalLikelihoods();
	// UpdateConditionalLikelihoods();
	// delete[] loglarray;
	return accepted;
}

void PhyloProcess::RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, double* loglarray, int& n)	{

	if (! from->isRoot())	{
		// ostringstream s1;
		// NewickTree::ToStream(s1);
		// GlobalAttach(down,up,from,fromup);
		GetTree()->Attach(down,up,from,fromup);
		// UpdateConditionalLikelihoods();
		// GlobalReset(0);
		double*** aux = condlmap[0];
		Reset(aux);
		for (const Link* link=up->Next(); link!=up; link=link->Next())	{
			if (link->isRoot())	{
				cerr << "ROOT\n";
				exit(1);
			}
			// GlobalMultiply(link,0);
			Multiply(GetConditionalLikelihoodVector(link),aux);
		}
		// GlobalPropagate(0,up->Out(),GetLength(up->GetBranch()));
		Propagate(aux,GetConditionalLikelihoodVector(up->Out()),GetLength(up->GetBranch()));
		double logl = ComputeNodeLikelihood(up->Out(),0);
		// double logl = GlobalComputeNodeLikelihood(up->Out(),0);
		// UpdateConditionalLikelihoods();
		// double logl2 = ComputeNodeLikelihood(up->Out(),aux);
		// cerr << logl << '\t' << logl2 << '\n';
		if (n >= GetNbranch())	{
			cerr << "branch overflow\n";
			exit(1);
		}
		loglarray[n] = logl;
		n++;
		// loglmap[pair<Link*,Link*>(from,fromup)] = logl;
		Link* tmp1 = GetTree()->Detach(down,up);
		// GlobalDetach(down,up,tmp1,tmp2);
		
		// UpdateConditionalLikelihoods();
		// ostringstream s2;
		// NewickTree::ToStream(s2);
		/*
		if (s1.str() != s2.str())	{
			cerr << "error: two different trees\n";
			cerr << s1.str() << '\n';
			cerr << s2.str() << '\n';
			exit(1);
		}
		if (from != tmp1)	{
			cerr << "error in gibbs: " << from << '\t' << tmp1 << '\n';
			exit(1);
		}
		if (fromup != tmp2)	{
			cerr << "error in gibbs: " << fromup << '\t' << tmp2 << '\n';
			exit(1);
		}
		*/
		// UpdateConditionalLikelihoods();
		// cerr << from << '\t' << tmp1 << '\t' << fromup << '\t' << tmp2 << '\n';
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// RecursiveGibbsSPRScan(link->Out(),from,down,up,loglmap);
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

	// not necessary, but
	// should be updated when entering this function
	// UpdateConditionalLikelihoods();
	// 
	// this is important
	// for the conditional likelihoods of the cut-then-regrafted subtree to be OK

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
		if (isinf(total))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total))	{
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
		// ostringstream s1;
		// NewickTree::ToStream(s1);
		GetTree()->Attach(down,up,from,fromup);
		// UpdateConditionalLikelihoods();
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
		// UpdateConditionalLikelihoods();
		// double logl2 = ComputeNodeLikelihood(up->Out(),aux);
		// cerr << logl << '\t' << logl2 << '\n';
		loglmap[pair<Link*,Link*>(from,fromup)] = logl;
		Link* tmp1 = GetTree()->Detach(down,up);
		// UpdateConditionalLikelihoods();
		// ostringstream s2;
		// NewickTree::ToStream(s2);
		/*
		if (s1.str() != s2.str())	{
			cerr << "error: two different trees\n";
			cerr << s1.str() << '\n';
			cerr << s2.str() << '\n';
			exit(1);
		}
		if (from != tmp1)	{
			cerr << "error in gibbs: " << from << '\t' << tmp1 << '\n';
			exit(1);
		}
		if (fromup != tmp2)	{
			cerr << "error in gibbs: " << fromup << '\t' << tmp2 << '\n';
			exit(1);
		}
		*/
		// UpdateConditionalLikelihoods();
		// cerr << from << '\t' << tmp1 << '\t' << fromup << '\t' << tmp2 << '\n';
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// RecursiveGibbsSPRScan(link->Out(),from,down,up,loglmap);
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

void PhyloProcess::SlaveGetMeanDiversity()	{

	double div = GetData()->GetTotalDiversity(sitemin,sitemax);
	MPI_Send(&div,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
}

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
	prop_arg args;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	args.from = GetLinkIndex(from);
	args.to = GetLinkIndex(to);
	args.condalloc = (condalloc) ? 1 : 0;
	args.time = time;
	MPI_Bcast(&args,1,Propagate_arg,0,MPI_COMM_WORLD);
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
	prop_arg args;
	args.time = m;
	args.condalloc = branch->GetIndex();
	MPI_Bcast(&args,1,Propagate_arg,0,MPI_COMM_WORLD);
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
	prop_arg alpha;
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
		MPI_Bcast(&alpha,1,Propagate_arg,0,MPI_COMM_WORLD);
		SlavePropose(alpha.condalloc,alpha.time);
		break;
	case RESTORE:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveRestore(n);
		break;
	case RESET:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[1] == 1) ? true : false;
		SlaveReset(arg[0],tvalue);
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
		MPI_Bcast(&alpha,1,Propagate_arg,0,MPI_COMM_WORLD);
		tvalue = (alpha.condalloc == 1) ? true : false;
		SlavePropagate(alpha.from,alpha.to,tvalue,alpha.time);
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
	case GETDIV:
		SlaveGetMeanDiversity();
		break;
	case CVSCORE:
		SlaveComputeCVScore();
		break;
	case SITELOGL:
		SlaveComputeSiteLogL();
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

	// MPI2
	// should send message to slaves for updating their siteprofilesuffstats
	// by calling UpdateSiteProfileSuffStat()
	// then collect all suff stats
	//
	// suff stats are contained in 2 arrays
	// int* branchlengthsuffstatcount
	// double* branchlengthsuffstatbeta
	assert(myid == 0);
	int i,j,nbranch = GetNbranch();
	MPI_Status stat;
	MESSAGE signal = UPDATE_BLENGTH;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<nbranch; ++i) {
		branchlengthsuffstatcount[i] = 0;
		branchlengthsuffstatbeta[i] = 0.0;
	}
	#ifdef BYTE_COM
	// should be summed over all slaves (reduced)
	int k,l;
	double x;
	unsigned char* bvector = new unsigned char[nbranch*(sizeof(int)+sizeof(double))];

	for(i=1; i<nprocs; ++i) {
		MPI_Recv(bvector,nbranch*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			l = 0;
			for(k=sizeof(int)-1; k>=0; --k) {
				l = (l << 8) + bvector[sizeof(int)*j+k]; 
			}
			branchlengthsuffstatcount[j] += l;
		}
		for(j=0; j<nbranch; ++j) {
			memcpy(&x,&bvector[sizeof(int)*nbranch+sizeof(double)*j],sizeof(double));
			branchlengthsuffstatbeta[j] += x;
		}
	}
	delete[] bvector;
	#else
	int ivector[nbranch];
	double dvector[nbranch];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,nbranch,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(ivector,nbranch,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchlengthsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(dvector,nbranch,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
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
	// finally, sync all processes on same suffstat values 
	/*
	MPI_Bcast(branchlengthsuffstatcount,GetNbranch(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(branchlengthsuffstatbeta,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/
	#endif
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
	int workload = GetNbranch();
	#ifdef BYTE_COM
	int i,n = 0;
	unsigned int j;
	unsigned char el_int[sizeof(int)],el_dbl[sizeof(double)];
	unsigned char* bvector = new unsigned char[workload*(sizeof(int)+sizeof(double))];
	for(i=0; i<workload; ++i) {
		convert(el_int,branchlengthsuffstatcount[i]);
		for(j=0; j<sizeof(int); ++j) {
			bvector[n] = el_int[j]; n++;
		}
	}
	for(i=0; i<workload; ++i) {
		convert(el_dbl,branchlengthsuffstatbeta[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_dbl[j]; n++;
		}
	}
	MPI_Send(bvector,workload*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;
	#else
	MPI_Send(branchlengthsuffstatcount,workload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(branchlengthsuffstatbeta,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	// finally, sync all processes on same suffstat values 
	/*
	MPI_Bcast(branchlengthsuffstatcount,GetNbranch(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(branchlengthsuffstatbeta,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/
	#endif
}

void PhyloProcess::GlobalUpdateSiteRateSuffStat()	{

	// MPI2
	// ask slaves to update siteratesuffstat
	// slaves should call UpdateSiteRateSuffStat()
	// then collect all suff stats
	// suff stats are contained in 2 arrays
	// int* siteratesuffstatcount
	// double* siteratesuffstatbeta
	// [site]
	assert(myid == 0);
	//cerr << "global update site rate\n";
	// each slave computes its array for sitemin <= site < sitemax
	// thus, one just needs to gather all arrays into the big master array 0 <= site < Nsite
	// (gather)
	int i,j,k,width,nalloc,smin[nprocs-1],smax[nprocs-1],workload[nprocs-1];
	MPI_Status stat;
	MESSAGE signal = UPDATE_SRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	width = GetNsite()/(nprocs-1);
	nalloc = 0;
	for(i=0; i<nprocs-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (nprocs-2)) smax[i] = GetNsite();
		workload[i] = smax[i] - smin[i];
		if (workload[i] > nalloc) nalloc = workload[i];
	}
	#ifdef BYTE_COM
	unsigned char* bvector = new unsigned char[nalloc*(sizeof(int)+sizeof(double))];
	int l,n;
	double x;
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(bvector,workload[i-1]*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
		n = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			l = 0;
			for(k=sizeof(int)-1; k>=0; --k) {
				l = (l << 8) + bvector[sizeof(int)*n+k]; 
			}
			siteratesuffstatcount[j] = l; n++;			
		}
		n = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			memcpy(&x,&bvector[workload[i-1]*sizeof(int)+n*sizeof(double)],sizeof(double));
			siteratesuffstatbeta[j] = x; n++;			
		}
	}
	delete[] bvector;
	#else
	int ivector[nalloc];
	double dvector[nalloc];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,workload[i-1],MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		k = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			siteratesuffstatcount[j] = ivector[k]; k++;
			
		}
	}
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,workload[i-1],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		k = 0;
		for(j=smin[i-1]; j<smax[i-1]; ++j) {
			siteratesuffstatbeta[j] = dvector[k]; k++;
		}
	}


	// finally, sync all processes on same suffstat values 
	// but this seems to corrupt something somewhere
	// I have tried to put this barrier, but to no avail
	// MPI_Barrier(MPI_COMM_WORLD);

	// activate these two lines, and the problem will appear
	/*
	MPI_Bcast(siteratesuffstatcount,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(siteratesuffstatbeta,GetNsite(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/
	#endif
}

void PhyloProcess::SlaveUpdateSiteRateSuffStat()	{

	UpdateSiteRateSuffStat();
	int i,workload = sitemax - sitemin;
	#ifdef BYTE_COM
	unsigned char* bvector = new unsigned char[workload*(sizeof(int)+sizeof(double))];
	unsigned char el_int[sizeof(int)],el_dbl[sizeof(double)];
	unsigned int j,n = 0;

	for(i=sitemin; i<sitemax; ++i) {
		convert(el_int,siteratesuffstatcount[i]);
		for(j=0; j<sizeof(int); ++j) {
			bvector[n] = el_int[j]; n++;
		}
	}
	for(i=sitemin; i<sitemax; ++i) {
		convert(el_dbl,siteratesuffstatbeta[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_dbl[j]; n++;
		}
	}		
	MPI_Send(bvector,workload*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;		
	#else
	int j = 0,ivector[workload];
	for(i=sitemin; i<sitemax; ++i) {
		ivector[j] = siteratesuffstatcount[i]; j++;
	}
	double dvector[workload];
	MPI_Send(siteratesuffstatcount,workload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	j = 0;
	for(i=sitemin; i<sitemax; ++i) {
		dvector[j] = siteratesuffstatbeta[i]; j++;
	}
	MPI_Send(siteratesuffstatbeta,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);


	// finally, sync all processes on same suffstat values 
	// but this seems to corrupt something somewhere
	// I have tried to put this barrier, but to no avail
	// MPI_Barrier(MPI_COMM_WORLD);

	// activate these two lines, and the problem will appear
	/*
	MPI_Bcast(siteratesuffstatcount,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(siteratesuffstatbeta,GetNsite(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/
	#endif
}

void PhyloProcess::GlobalGetMeanSiteRate()	{

	if (! meansiterate)	{
		meansiterate = new double[GetNsite()];
	}

	assert(myid == 0);
	int i,width,smin[nprocs-1],smax[nprocs-1],workload[nprocs-1];
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

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	int cv = 0;
	int sitelogl = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;

	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-div")	{
				ppred = 2;
			}
			else if (s == "-comp")	{
				ppred = 3;
			}
			else if (s == "-nsub")	{
				ppred = 4;
			}
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-ppredrate")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rateprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rateprior = 0;
				}
				else	{
					cerr << "error after ppredrate: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredprofile")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					profileprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					profileprior = 0;
				}
				else	{
					cerr << "error after ppredprofile: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredroot")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rootprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rootprior = 0;
				}
				else	{
					cerr << "error after ppredroot: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				string tmp = argv[i];
				if (IsFloat(tmp))	{
					until = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		MESSAGE signal = KILL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}
	if (burnin == -1)	{
		burnin = GetSize() / 5;
	}

	if ((GetNprocs() == 1) && (ppred || cv || sitelogl))	{
		cerr << "error : should run readpb_mpi in mpi mode, with at least 2 processes\n";
		MPI_Finalize();
		exit(1);
	}

	if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
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

void PhyloProcess::PostPred(int ppredtype, string name, int burnin, int every, int until, int inrateprior, int inprofileprior, int inrootprior)	{

	GlobalSetRatePrior(inrateprior);
	GlobalSetProfilePrior(inprofileprior);
	GlobalSetRootPrior(inrootprior);

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	double* obstaxstat = new double[GetNtaxa()];
	SequenceAlignment* datacopy  = new SequenceAlignment(GetData());
	double obs = 0;
	double obs2 = 0;
	if (ppredtype == 2)	{
		obs = data->GetMeanDiversity();
		// obs = GlobalGetMeanDiversity();
	}
	else if (ppredtype == 3)	{
		obs = GetObservedCompositionalHeterogeneity(obstaxstat,obs2);
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

		// output tree
		ostringstream s;
		s << name << "_ppred" << samplesize << ".tree";
		ofstream os(s.str().c_str());
		SetNamesFromLengths();
		RenormalizeBranchLengths();
		GetTree()->ToStream(os);
		DenormalizeBranchLengths();
		os.close();

		MPI_Status stat;
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
	cerr << '\n';
}

void PhyloProcess::GlobalSetTestData()	{
	testnsite = testdata->GetNsite();
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
		cout << "before FromStream...\n";
		cout.flush();
		FromStream(is);
		cout << "after FromStream...\n";
		cout.flush();
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
		// Trace(cerr);
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

void PhyloProcess::ReadSiteLogL(string name, int burnin, int every, int until)	{

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

	double* tmp = new double[GetNsite()];
	double* mean = new double[GetNsite()];
	vector<double>* logl = new vector<double>[GetNsite()];

	for (int i=0; i<GetNsite(); i++)	{
		mean[i] = 0;
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
		// Trace(cerr);
		MPI_Status stat;
		MESSAGE signal = SITELOGL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		double total = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			// for (int i=0; i<GetNsite(); i++)	{
			for (int j=smin[i-1]; j<smax[i-1]; j++)	{
			// for (int i=GetSiteMin(i); i<GetSiteMax(i); i++)	{
				logl[j].push_back(tmp[j]);
				mean[j] += tmp[j];
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

	double meancpo = 0;
	double varcpo = 0;
	double* cpo = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		double min = 0;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			if (min > (*j))	{
				min = *j;
			}
		}
		double hmean = 0;
		int count = 0;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			hmean += exp(min - (*j));
			count++;
		}
		hmean /= count;
		cpo[i] = min - log(hmean);
		meancpo += cpo[i];
		varcpo += cpo[i] * cpo[i];
	}
	meancpo /= GetNsite();
	varcpo /= GetNsite();
	varcpo -= meancpo * meancpo;

	ofstream os((name + ".sitelogl").c_str());
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		mean[i] /= samplesize;
		total += mean[i];
		os << i+1 << '\t' << mean[i] << '\t' << cpo[i] << '\n';
	}

	ofstream cos((name + ".cpo").c_str());
	cos << "posterior mean ln L : " << total << '\n';
	cos << "CPO : " << GetNsite() * meancpo << '\t' << meancpo << '\t' << sqrt(varcpo) << '\n';
	
	cerr << '\n';
	cerr << "posterior mean ln L : " << total << '\n';
	cerr << "site-specific posterior mean ln L in " << name << ".sitelogl\n";
	cerr << "CPO: " << GetNsite() * meancpo << '\t' << meancpo << '\t' << sqrt(varcpo) << '\n';
	cerr << '\n';

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
		MPI_Status stat;
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
	MPI_Status stat;
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
