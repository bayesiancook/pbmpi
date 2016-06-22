
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AAMUTSELSITESPECIFICPHYLO_H
#define AAMUTSELSITESPECIFICPHYLO_H

//#include <cassert>
#include "AAMutSelSiteSpecificSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"
//#include "Parallel.h"

class AAMutSelSiteSpecificPhyloProcess : public virtual AAMutSelSiteSpecificSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	// s'inspirer de GeneralPathSuffStatGTRPhyloProcess
	// et GeneralPathSuffStatRASCATGTRPhyloProcess

	public:

	AAMutSelSiteSpecificPhyloProcess(string indatafile, string treefile, GeneticCodeType incodetype, int ncat, int me, int np)	{
		myid = me;
		nprocs = np;

		datafile = indatafile;
		codetype = incodetype;
		SequenceAlignment* nucdata = new FileSequenceAlignment(datafile,0,myid);
		CodonSequenceAlignment* codondata = new CodonSequenceAlignment(nucdata,true,codetype);
		CodonStateSpace* statespace = codondata->GetCodonStateSpace();
		const TaxonSet* taxonset = codondata->GetTaxonSet();

		if (treefile == "None")	{
			tree = new Tree(taxonset);
			if (myid == 0)	{
				tree->MakeRandomTree();
				GlobalBroadcastTree();
			}
			else	{
				SlaveBroadcastTree();
			}
		}
		else	{
			tree = new Tree(treefile);
		}
		tree->RegisterWith(taxonset,myid);
		
		int insitemin = -1,insitemax = -1;
		if (myid > 0) {
			int width = codondata->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = codondata->GetNsite();
			}
			else {
				insitemax = myid*width;
			}
		}

		//Create(tree,codondata,ncat,insitemin,insitemax,statespace);
		//
		//  The following line is the important difference that sets the number of components to the number of sites.
		Create(tree,codondata,codondata->GetNsite(),insitemin,insitemax,statespace);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	AAMutSelSiteSpecificPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		is >> datafile;
		is >> codetype;
		SequenceAlignment* nucdata = new FileSequenceAlignment(datafile,0,myid);
		CodonSequenceAlignment* codondata = new CodonSequenceAlignment(nucdata,true,codetype);
		CodonStateSpace* statespace = codondata->GetCodonStateSpace();
		const TaxonSet* taxonset = codondata->GetTaxonSet();

		int insitemin = -1,insitemax = -1;
		if (myid > 0) {
			int width = codondata->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = codondata->GetNsite();
			}
			else {
				insitemax = myid*width;
			}
		}

		tree = new Tree(taxonset);
		if (myid == 0)	{
			tree->ReadFromStream(is);
			GlobalBroadcastTree();
		}
		else	{
			SlaveBroadcastTree();
		}
		tree->RegisterWith(taxonset,0);

		Create(tree,codondata,1,insitemin,insitemax,statespace);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	virtual ~AAMutSelSiteSpecificPhyloProcess()	{
		Delete();
	}

	// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
	// are model dependent
	// should be implemented in .cpp file
        virtual void SlaveExecute(MESSAGE);
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		return 0;
	}

	double GetLogLikelihood()	{
		return logL;
	}

	void TraceHeader(ostream& os)	{
		os << "lnL\tlength\tNmode\tNocc\tnucsA\tnucsC\tnucsT\tnucsG\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT\tstatent";
		os << "\ttotaltime";
		os << "\tpruning\tsuffstat\tunfold\tcollapse";
		os << "\n";
	}

	void Trace(ostream& os)	{
		os << GetLogLikelihood();
		os << '\t' << GetTotalLength();
		os << '\t' << GetNcomponent();
		os << '\t' << GetNOccupiedComponent();
		os << '\t' << GetNucStat(0) << '\t' << GetNucStat(1) << '\t' << GetNucStat(2) << '\t' << GetNucStat(3);
		os << '\t' << GetNucRR(0) << '\t' << GetNucRR(1) << '\t' << GetNucRR(2) << '\t' << GetNucRR(3) << '\t' << GetNucRR(4) << '\t' << GetNucRR(5);
		os << '\t' << GetStatEnt();
		os << '\t' << ((int) (chronototal.GetTime() / 1000));

		if (chronototal.GetTime())	{
			os << '\t' << ((int) (100 * chronopruning.GetTime() /chronototal.GetTime()));
			os << '\t' << ((int) (100 * chronosuffstat.GetTime() /chronototal.GetTime()));
			os << '\t' << ((int) (100 * chronounfold.GetTime() /chronototal.GetTime()));
			os << '\t' << ((int) (100 * chronocollapse.GetTime() /chronototal.GetTime()));
		}
		else	{
			os << '\t' << '-';
			os << '\t' << '-';
			os << '\t' << '-';
			os << '\t' << '-';
		}
		os << '\n';
	}

	void ToStream(ostream& os)	{
		os << datafile << '\n';
		os << codetype << '\n';
		GetTree()->ToStream(os);
		GammaBranchProcess::ToStream(os);
		AAMutSelSiteSpecificProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AAMutSelSiteSpecificProfileProcess::FromStream(is);
	}

	// primary scheduler

	double Move(double tuning = 1.0)	{
		// cerr << "unfold\n";
		chronototal.Start();
		chronopruning.Start();
		//cerr << "bl\n";
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		//cerr << "gspr\n";
		//GibbsSPR(10);
		MoveTopo(5,0);
		//cerr << "collapse\n";
		chronopruning.Stop();

		chronosuffstat.Start();

		chronocollapse.Start();
		GlobalCollapse();
		chronocollapse.Stop();
		//cerr << "branch\n";
		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();
		AAMutSelSiteSpecificProfileProcess::Move(tuning,1,10);
		chronosuffstat.Stop();

		chronounfold.Start();
		GlobalUnfold();
		chronounfold.Stop();

		chronototal.Stop();
		//cerr << "ok\n";
		return 1;
	}


	// This function sets the allocation of one profile per site
	void SampleAlloc()	{
		Ncomponent = AAMutSelSiteSpecificProfileProcess::GetNsite();
		for (int i=0; i<AAMutSelSiteSpecificProfileProcess::GetNsite(); i++)	{
			CreateComponent(i);
			AddSite(i,i);

		}
	}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int sitemin, int sitemax)	{
		cerr << "In two-argument Create of AAMutSelSiteSpecificPhyloProcess.  Should not be here.\n";
		exit(1);
	}

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat, int sitemin, int sitemax, CodonStateSpace* instatespace)	{
		AAMutSelSiteSpecificSubstitutionProcess::Create(indata->GetNsite(),Naa,ncat,sitemin,sitemax,instatespace);
		GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,Naa,sitemin,sitemax); // added Naa here.
		GammaBranchProcess::Create(intree);
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AAMutSelSiteSpecificSubstitutionProcess::Delete();
		GammaBranchProcess::Delete();
	}

	GeneticCodeType codetype;
	CodonStateSpace* statespace;

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};

// enfin, le PhyloProcess ainsi construit peut etre instancie dans le main.cpp



#endif

