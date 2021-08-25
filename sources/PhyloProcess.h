
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "ZippedSequenceAlignment.h"

#include "SubstitutionProcess.h"
#include "BranchProcess.h"

#include "Parallel.h"

#include <map>
#include <vector>

class PhyloProcess : public virtual SubstitutionProcess, public virtual BranchProcess {

	public:

	virtual void WaitLoop();

	virtual void SlaveExecute(MESSAGE);

        virtual void SlaveRoot(int);
	virtual void SlaveGibbsSPRScan(int,int);
	virtual void SlaveLikelihood(int,int);
	virtual void SlavePropose(int,double);
	virtual void SlaveRestore(int);
	virtual void SlaveReset(int,bool);
	virtual void SlaveSMultiply(int,bool);
	virtual void SlaveMultiply(int,int,bool);
	virtual void SlaveInitialize(int,int,bool);
	virtual void SlavePropagate(int,int,bool,double);
	virtual void SlaveDetach(int,int);
	virtual void SlaveAttach(int,int,int,int);

	// virtual void SlaveUpdate();

	// default constructor: pointers set to nil
	PhyloProcess() : missingmap(0), sitecondlmap(0), condlmap(0), siteratesuffstatcount(0), siteratesuffstatbeta(0), branchlengthsuffstatcount(0), branchlengthsuffstatbeta(0), condflag(false), data(0), bkdata(0), steppingrank(0), minsitecutoff(-1), maxsitecutoff(-1), myid(-1), nprocs(0), size(0), version("1.9"), totaltime(0), dataclamped(1), rateprior(0), profileprior(0), rootprior(1), topoburnin(0) {
		fixbl = 0;
		sitesuffstat = 1;
	}
	virtual ~PhyloProcess() {}

	string GetVersion() {return version;}
	/*
	int GetSiteMin(int proc);
	int GetSiteMax(int proc);
	*/

	// performs one full cycle of MCMC
	// returns average success rate
	virtual double Move(double tuning = 1.0) = 0;

    virtual void GlobalSetEmpiricalPrior(istream& is)  {
        cerr << "in PhyloProcess::GlobalSetEmpiricalPrior\n";
        exit(1);
    }

    virtual void SlaveSetEmpiricalPrior()   {
        cerr << "in PhyloProcess::SlaveSetEmpiricalPrior\n";
        exit(1);
    }

    void SetEmpiricalFrac(double infrac)    {
        SetRateFrac(infrac);
        SetProfileFrac(infrac);
        SetLengthFrac(infrac);
    }

    void GlobalSetEmpiricalFrac(double infrac);
    void SlaveSetEmpiricalFrac();

	void SetFixBL(int in)	{
		fixbl = in;
	}

	// sample from prior
	virtual void Sample()	{
		SampleRate();
		SampleLength();
		SampleProfile();
	}

    virtual void PriorSample()  {
        PriorSampleRate();
        PriorSampleLength();
        PriorSampleProfile();
    }

	// print out the first line (header) of the trace file
	virtual void TraceHeader(ostream& os) = 0;

	// print out one line of trace (summary statistics such as logprob, treelength, totaltime, etc)
	virtual void Trace(ostream& os) = 0;

	virtual void Monitor(ostream& os)  {
		os << "matrix uni" << '\t' << SubMatrix::GetUniSubCount() << '\n';
		os << "inf prob  " << '\t' << GetInfProbCount() << '\n';
		os << "stat inf  " << '\t' << GetStatInfCount() << '\n';
	}

	virtual void ToStreamHeader(ostream& os)	{
		os << version << '\n';
		propchrono.ToStream(os);
		chronototal.ToStream(os);
	}

	virtual void FromStreamHeader(istream& is)	{
		is >> version;
		if (atof(version.substr(0,3).c_str()) < 1.2)	{
			cerr << "error: chain run under older version not anymore supported: " << version << '\n';
			exit(1);
		}
		propchrono.FromStream(is);
		chronototal.FromStream(is);
	}

	virtual void ToStream(ostream& os)	{
		cerr << "error: in phyloprocess::ToStream\n";
		exit(1);
	}

	virtual void FromStream(istream& is)	{
		cerr << "error: in phyloprocess::FromStream\n";
		exit(1);
	}

	// translation tables : from pointers of type Link* Branch* and Node* to their index and vice versa
	// this translation is built when the Tree::RegisterWithTaxonSet method is called (in the model, in PB.cpp)
	Link* GetLink(int linkindex)	{
		if (! linkindex)	{
			//return GetRoot();
			return 0;
		}
		return GetTree()->GetLink(linkindex);
	}

	Link* GetLinkForGibbs(int linkindex)	{
		if (! linkindex)	{
			return GetRoot();
		}
		return GetTree()->GetLink(linkindex);
	}

	const Branch* GetBranch(int branchindex)	{
		return GetTree()->GetBranch(branchindex);
	}

	const Node* GetNode(int nodeindex)	{
		return GetTree()->GetNode(nodeindex);
	}


	int GetLinkIndex(const Link* link)	{
		return link ? link->GetIndex() : 0;
	}

	int GetBranchIndex(const Branch* branch)	{
		if (! branch)	{
			return 0;
			/*
			cerr << "error in get branch index\n";
			exit(1);
			*/
		}
		return branch->GetIndex();
	}

	int GetNodeIndex(const Node* node)	{
		return node->GetIndex();
	}

	/*
	// not useful. link->GetIndex() is used instead (same for branches and nodes)
	int GetBranchIndex(const Branch* branch);
	int GetNodeIndex(const Node* node);
	*/

	const TaxonSet* GetTaxonSet() const {return data->GetTaxonSet();}

	void GlobalUpdateConditionalLikelihoods();
	double GlobalComputeNodeLikelihood(const Link* from, int auxindex = -1);

    /*
    void GlobalCreateSiteDataStructures();
    void SlaveCreateSiteDataStructures();

    void GlobalDeleteSiteDataStructures();
    void SlaveDeleteSiteDataStructures();
    */

	void CreateSiteConditionalLikelihoods();
	void DeleteSiteConditionalLikelihoods();

	virtual void PrepareSiteLogLikelihood(int site) {
		cerr << "in default PrepareSiteLogLikelihood\n";
		exit(1);
	}

    // for simple models: sums over rate allocations
    // for mixture models: sums over rate and profile allocations
    // + importance sampling for sbdp models
    // this can be parallelized:
    // each slave can deal with some of the components
    //

    virtual double GlobalGetSiteSteppingLogLikelihood(int site, int nrep, int restore);
    virtual void SlaveGetSiteSteppingLogLikelihood();

	virtual double SiteLogLikelihood(int site);
	void SitePostOrderPruning(int site, const Link* from);

    void SetSteppingFraction(int cutoff1, int cutoff2);
    void GlobalSetSteppingFraction(int cutoff1, int cutoff2);
    void SlaveSetSteppingFraction();

    void GlobalPrepareStepping(string name, int size, int rand);
    void SlavePrepareStepping();

    bool ActiveSite(int i);

	// returns total number of taxa in the analysis
	int GetNtaxa()	{
		return GetData()->GetNtaxa();
	}

	protected:

	SequenceAlignment* GetData() {return data;}
	StateSpace* GetStateSpace() {return data->GetStateSpace();}
	virtual int GetGlobalNstate() {return GetStateSpace()->GetNstate();}

	double* GetEmpiricalFreq()	{
		return empfreq;
	}

	int isObserved(int site, int k)	{
		int obs = 0;
		for (int j=0; j<GetNtaxa(); j++)	{
			if (GetData()->GetState(j,site) == k)	{
				obs = 1;
			}
		}
		return obs;
	}

	// the following methods are particularly important for MPI
	// Create / Delete / Unfold and Collapse should probably be specialized
	// according to whether this is a slave or the master processus
	virtual void Create(Tree* intree, SequenceAlignment* indata, int indim);
	virtual void Delete();

	public :

	void GlobalSetRatePrior(int inrateprior);
	void SlaveSetRatePrior();
	
	void GlobalSetProfilePrior(int inprofileprior);
	void SlaveSetProfilePrior();
	
	void GlobalSetRootPrior(int inrootprior);
	void SlaveSetRootPrior();
	

	virtual void Unfold();
	virtual void Collapse();

	void GlobalBroadcastTree();
	void SlaveBroadcastTree();

	// MPI Global dispatcher functions
	// all methods with a "Global" prefix are the ones through which MPI parallelization is supposed to go
	// these methods should dispatch the computation over slaves, and possibly, collect the result of the computation
	// in the non MPI version, they just call their non MPI counterpart (or nearly so)

	public: 

	virtual void QuickUpdate()	{

		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		
		GlobalUpdateConditionalLikelihoods();
		GlobalCollapse();
		GlobalUnfold();
	}

	void GlobalSimulateForward();
	void SimulateForward();
	virtual void RecursiveSimulateForward(const Link* from);

	virtual void SetDataFromLeaves()	{
		for (int i=sitemin; i<sitemax; i++)	{
			RecursiveSetDataFromLeaves(i,GetRoot());
		}
	}

	void RecursiveSetDataFromLeaves(int site, const Link* from)	{
		if (from->isLeaf())	{
			if (data->GetBKState(GetNodeIndex(from->GetNode()),site) != -1)	{
				data->SetState(GetNodeIndex(from->GetNode()),site,nodestate[GetNodeIndex(from->GetNode())][site]);
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetDataFromLeaves(site, link->Out());
		}
	}

	void GlobalUnclamp();
	void GlobalRestoreData();
	void GlobalSetDataFromLeaves();
    /*
	double GlobalGetMeanDiversity();
	double GlobalGetMeanSquaredFreq();
	double GlobalGetMeanFreqVariance();
    */
	void SlaveUnclamp();
	void SlaveRestoreData();
	void SlaveSetDataFromLeaves();
    /*
	void SlaveGetMeanDiversity();
	void SlaveGetMeanSquaredFreq();
	void SlaveGetMeanFreqVariance();
    */

	void GlobalSetNodeStates();
	void SlaveSetNodeStates();
	void WriteNodeStates(ostream& os, const Link* from);

	void ReadAncestral(string name, int burnin, int every, int until);
	void SlaveComputeStatePostProbs();
	void ComputeStatePostProbs(double*** statepostprobs, const Link* from, int auxindex);
	void RecursiveComputeStatePostProbs(double*** statepostprob, const Link* from, int auxindex);
	void WriteStatePostProbs(double*** statepostprob, string name, const Link* from);

	virtual void ReadPB(int argc, char* argv[]);
	virtual void Read(string name, int burnin, int every, int until);
	virtual void ReadSiteLogL(string name, int burnin, int every, int until, int verbose);
	virtual void ReadCV(string testdatafile, string name, int burnin, int every, int until, int iscodon = 0, GeneticCodeType codetype = Universal);
	virtual void ReadSiteCV(string testdatafile, string name, int burnin, int every, int until, int iscodon = 0, GeneticCodeType codetype = Universal);
	virtual void AllPostPred(string name, int burnin, int every, int until, int rateprior, int profileprior, int rootprior);
	virtual void PostPred(int ppredtype, string name, int burnin, int every, int until, int rateprior, int profileprior, int rootprior, int savetrees);

	void ReadSiteRates(string name, int burnin, int every, int until);

	// The following methids are here to write the mappings.
	void ReadMap(string name, int burnin, int every, int until);
	void ReadPostPredMap(string name, int burnin, int every, int until);
	void GlobalWriteMappings(string name);
	virtual void SlaveWriteMappings();
	void WriteTreeMapping(ostream& os, const Link* from, int i);



	virtual void GlobalSetTestData();
	virtual void SlaveSetTestData();
	void SetTestSiteMinAndMax();
	virtual void SlaveComputeCVScore() {
		cerr << "slave compute cv score\n";
		exit(1);
	}
	
	virtual void SlaveComputeSiteLogL() {
		cerr << "slave compute site logL\n";
		exit(1);
	}
	
	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
        return LogBranchPrior() + LogRatePrior() + LogProfilePrior();
	}

	double GetLogLikelihood()	{
		return logL;
	}

    double GlobalGetFullLogLikelihood();

    virtual double GlobalGetSteppingLogLikelihood(int nrep, int restore) {
        double tot = 0;
        for (int i=0; i<GetNsite(); i++)    {
            if (ActiveSite(i))  {
                double tmp = GlobalGetSiteSteppingLogLikelihood(i, nrep, restore);
                tot += tmp;
            }
        }
        return tot;
        // return GlobalGetFullLogLikelihood();
    }

	virtual void GlobalUnfold();
	virtual void GlobalCollapse();

    void GlobalResetAllConditionalLikelihoods();
    void SlaveResetAllConditionalLikelihoods();

	void GlobalReset(const Link* from, bool condalloc = false);
	void GlobalMultiply(const Link* from, const Link* to, bool condalloc = false);
	void GlobalMultiplyByStationaries(const Link* from, bool condalloc = false);
	void GlobalInitialize(const Link* from, const Link* link, bool condalloc = false);

	void GlobalPropagate(const Link* from, const Link* to, double time, bool condalloc = false);
	double GlobalProposeMove(const Branch* branch, double tuning);
	void GlobalRestore(const Branch* branch);

	void GlobalRootAtRandom();
	Link* GlobalDetach(Link* down, Link* up);
	// void GlobalDetach(Link* down, Link* up);
	void GlobalAttach(Link* down, Link* up, Link* fromdown, Link* fromup);

	//NNI functions ( in NNI.cpp )
	void RecursiveGibbsNNI(Link* from, double tuning, int type, int& success, int& moves);
	double GibbsNNI(double tuning, int);
	int  GlobalNNI(Link*,double,int);
	void GlobalKnit(Link*);
	void GlobalPropagateOverABranch(Link*);
	void SlaveNNI(Link*,int);
	void PropagateOverABranch(const Link*);
	int SlaveSendNNILikelihood(Link*);
	double SendRandomBranches(Link*,double,Link**&, int);
	double MoveTopo(int spr, int nni);

	// MCMC on branch lengths
	double BranchLengthMove(double tuning);
	double NonMPIBranchLengthMove(double tuning);

	// MCMC on topology
	double GibbsSPR(int nrep);
	double OldMPIGibbsSPR(int nrep);

	void CheckLikelihood();
	void GlobalCheckLikelihood();

	virtual void CreateSuffStat();
	virtual void DeleteSuffStat();

	// if  i == -1 then create auxiliary array, otherwise use condlmap[auxindex] as the auxiliary array
	double ComputeNodeLikelihood(const Link* from, int auxindex = -1);

	// assumes that rate allocations have already been defined
	// and that conditional likelihoods are updated
	// those conditional likelihoods will be corrupted
	void SampleNodeStates();
	void SampleNodeStates(const Link* from, double*** aux);

	// assumes that states at nodes have been sampled (using ResampleState())
	void SampleSubstitutionMappings(const Link* from);

	// conditional likelihood propagations
	void PostOrderPruning(const Link* from, double*** aux);
	void PreOrderPruning(const Link* from, double*** aux);
	void RecursiveComputeLikelihood(const Link* from, int auxindex, vector<double>& logl);
	void GlobalRecursiveComputeLikelihood(const Link* from, int auxindex, vector<double>& logl);

	double RecursiveBranchLengthMove(const Link* from, double tuning, int& n);
	double RecursiveNonMPIBranchLengthMove(const Link* from, double tuning, int& n);
	double LocalBranchLengthMove(const Link* from, double tuning);
	double LocalNonMPIBranchLengthMove(const Link* from, double tuning);



	int GibbsSPR();
	// double GibbsSPR();
	void GlobalGibbsSPRScan(Link* down, Link* up, double* loglarray);
	void RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, double* loglarray, int& n);
	void RecursiveGibbsFillMap(Link* from, Link* fromup, map<pair<Link*,Link*>,double>& loglmap, double* loglarray, int& n);

	double OldMPIGibbsSPR();
	void RecursiveOldMPIGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap);

	double NonMPIGibbsSPR();
	void RecursiveNonMPIGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap);

	// various protected accessors
	// used for computation and maintenance from within PhyloProcess classes
	const int* GetData(int index)	{
		return data->GetState(index);
	}

	const int* GetData(const Link* from)	{
		if (! from->isLeaf())	{
			cerr << "error in PhyloProcess::GetData\n";
			exit(1);
		}
		return data->GetState(from->GetNode()->GetIndex());
	}

	/*
	int GetData(const Link* from, int site)	{
		return GetData(from)[site];
	}
	*/

	int* GetStates(const Node* node)	{
		return nodestate[GetNodeIndex(node)];
	}

	/*
	int GetState(const Node* node, int site)	{
		return nodestate[node][site];
	}

	void SetState(const Node* node, int site, int state)	{
		nodestate[node][site] = state;
	}

	BranchSitePath* GetBranchSitePath(const Branch* branch, int site)	{
		return submap[branch][site];
	}
	*/

	void CreateConditionalLikelihoods();
	void DeleteConditionalLikelihoods();
	virtual void UpdateConditionalLikelihoods();

	double*** GetConditionalLikelihoodVector(const Link* link)	{
		return condlmap[GetLinkIndex(link)];
	}

	void CreateMappings();
	void DeleteMappings();

	void CreateNodeStates();
	void DeleteNodeStates();
	
	void CreateMissingMap();
	void DeleteMissingMap();
	void FillMissingMap();
	void BackwardFillMissingMap(const Link* from);
	void ForwardFillMissingMap(const Link* from, const Link* up);

	int** missingmap;

	// sufficient statistics for rates and branch lengths (do not depend on the model)
	int GetSiteRateSuffStatCount(int site) {return siteratesuffstatcount[site];}
	double GetSiteRateSuffStatBeta(int site) {return siteratesuffstatbeta[site];}

	/*
	int GetBranchLengthSuffStatCount(const Branch* branch) {return branchlengthsuffstatcount[GetBranchIndex(branch)];}
	double GetBranchLengthSuffStatBeta(const Branch* branch) {return branchlengthsuffstatbeta[GetBranchIndex(branch)];}
	*/

	int GetBranchLengthSuffStatCount(int index) {return branchlengthsuffstatcount[index];}
	double GetBranchLengthSuffStatBeta(int index) {return branchlengthsuffstatbeta[index];}

	void GlobalUpdateSiteRateSuffStat();
	void GlobalUpdateBranchLengthSuffStat();

	void SlaveUpdateSiteRateSuffStat();
	void SlaveUpdateBranchLengthSuffStat();

	void GlobalGetMeanSiteRate();
	void SlaveSendMeanSiteRate();

	virtual int CountMapping();
	virtual int CountMapping(int site);
	virtual int GlobalCountMapping();
	void SlaveCountMapping();


	virtual double GetObservedCompositionalHeterogeneity(double* taxstat, double& meandist)	{
		return data->CompositionalHeterogeneity(taxstat,0,meandist);
	}

	virtual double GetCompositionalHeterogeneity(double* taxstat, double& meandist)	{
		return data->CompositionalHeterogeneity(taxstat,0,meandist);
	}

	virtual int GetNprocs() {
		return nprocs;
	}
	virtual int GetMyid() {
		return myid;
	}

	double*** sitecondlmap;
	double**** condlmap;
	BranchSitePath*** submap;
	int** nodestate;

	// sufficient statistics for rates and branch lengths (do not depend on the model)
	int* siteratesuffstatcount;
	double* siteratesuffstatbeta;
	// map<const Branch*, int> branchlengthsuffstatcount;
	// map<const Branch*, double> branchlengthsuffstatbeta;
	int* branchlengthsuffstatcount;
	double* branchlengthsuffstatbeta;

	bool condflag;

	SequenceAlignment* data;
	string datafile;
    SequenceAlignment* bkdata;
    int* steppingrank;
    double minsitecutoff;
    double maxsitecutoff;

	double* empfreq;

	Chrono propchrono;
	Chrono chronototal;
	/*
	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronocollapse;
	Chrono chronounfold;
	*/

	double* loglarray;

	// MPI
	int myid,nprocs;

	int size;
	void IncSize()	{size++;}
	int GetSize() {return size;}
	void SetSize(int insize) {size = insize;}
	int GetIndex() {
		// prior to version 1.7, samples are 1-indexed in trace files
		if (atof(version.substr(0,3).c_str()) < 1.7)	{
			return size;
		}
		// they are 0-indexed starting with version 1.7
		return size - 1;
	}

	void SetTopoBurnin(int intopoburnin)	{
		topoburnin = intopoburnin;
	}

	double GetNormFactor() {return GetNormalizationFactor();}

	string version;
	double totaltime;

	SequenceAlignment* testdata;
	int testnsite;
	int testsitemin;
	int testsitemax;
	int bksitemax;

	int dataclamped;
	int rateprior;
	int profileprior;
	int rootprior;

	int topoburnin;
	int fixbl;

	int sitesuffstat;

	int fixtopo;
	int NSPR;
	int NNNI;
	int dc;

};



#endif
