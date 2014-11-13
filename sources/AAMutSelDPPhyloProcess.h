
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AAMUTSELDPPHYLO_H
#define AAMUTSELDPPHYLO_H

//#include <cassert>
#include "AAMutSelDPSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"
//#include "Parallel.h"

class AAMutSelDPPhyloProcess : public virtual AAMutSelDPSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	// s'inspirer de GeneralPathSuffStatGTRPhyloProcess
	// et GeneralPathSuffStatRASCATGTRPhyloProcess

	public:

	AAMutSelDPPhyloProcess(string indatafile, string treefile, GeneticCodeType incodetype, int infixtopo, int infixbl, int inkappaprior, int indc, int me, int np)	{
	//AAMutSelDPPhyloProcess(Tree* intree, SequenceAlignment* indata, CodonStateSpace* instatespace)	{
		myid = me;
		nprocs = np;
		fixtopo = infixtopo;
		fixbl = infixbl;
		dc = indc;
		kappaprior = inkappaprior;
		datafile = indatafile;
		codetype = incodetype;
		SequenceAlignment* nucdata = new FileSequenceAlignment(datafile,0,myid);
		CodonSequenceAlignment* codondata = new CodonSequenceAlignment(nucdata,true,codetype);
		CodonStateSpace* statespace = codondata->GetCodonStateSpace();
		if (dc)	{
			codondata->DeleteAAConstantSites();
		}

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

		Create(tree,codondata,insitemin,insitemax,statespace);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	AAMutSelDPPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		is >> codetype;
		is >> fixtopo;
		is >> fixbl;
		is >> dc;
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

		Create(tree,codondata,insitemin,insitemax,statespace);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}

		//

		//Create(intree,indata,instatespace);
		//Unfold();
		//Collapse();
	}	
	virtual ~AAMutSelDPPhyloProcess()	{
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
		os << "lnL\tlength\tNmode\tnucsA\tnucsC\tnucsT\tnucsG\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT\tstatent";
		os << "\totaltime";
		os << "\tpruning\tsuffstat\tunfold\tcollapse";
		os << "\n";
	}

	void Trace(ostream& os)	{
		os << GetLogLikelihood();
		os << '\t' << GetTotalLength();
		os << '\t' << GetNcomponent();
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

	void ToStreamHeader(ostream& os)	{
		os << datafile << '\n';
		os << codetype << '\n';
		os << fixtopo << '\n';
		os << fixbl << '\n';
		os << dc << '\n';
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		AAMutSelDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AAMutSelDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	// primary scheduler

	double Move(double tuning = 1.0)	{
		// cerr << "unfold\n";
		chronototal.Start();
		chronounfold.Start();
		chronounfold.Stop();
		chronopruning.Start();
		//cerr << "bl\n";
		if (! fixbl)	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}
		//cerr << "gspr\n";
		if (! fixtopo)	{
			MoveTopo(50,0);
			//GibbsSPR(5);
		}
		chronopruning.Stop();
		chronosuffstat.Start();
		chronocollapse.Start();
		//cerr << "collapse\n";
		GlobalCollapse();
		chronocollapse.Stop();
		//cerr << "branch\n";
		if (! fixbl)	{
			GammaBranchProcess::Move(tuning,10);
		}
		//cerr << "GlobalUpdateParameters()\n";
		GlobalUpdateParameters();
		//cerr << "nucrr, nucstat, aaprofiles\n";
		AAMutSelDPProfileProcess::Move(tuning,1,10);
		chronototal.Stop();
		chronosuffstat.Stop();
		GlobalUnfold();
		//cerr << "ok\n";
		return 1;
	}


	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int sitemin, int sitemax)	{
		cerr << "In two-argument Create of AAMutSelFinitePhyloProcess.  Should not be here.\n";
		exit(1);
	}
	//virtual void Create(Tree* intree, SequenceAlignment* indata)	{
	//	cerr << "In two-argument Create of AAMutSelDPPhyloProcess.  Should not be here.\n";
	//	exit(1);
	//}

	virtual void Create(Tree* intree, SequenceAlignment* indata, int sitemin, int sitemax, CodonStateSpace* instatespace)	{
		AAMutSelDPSubstitutionProcess::Create(indata->GetNsite(),Naa,sitemin,sitemax,instatespace);
		GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,Naa,sitemin,sitemax); // added Naa here.
		GammaBranchProcess::Create(intree);
	}
	//virtual void Create(Tree* intree, SequenceAlignment* indata, CodonStateSpace* instatespace)	{
	//	AAMutSelDPSubstitutionProcess::Create(indata->GetNsite(),Naa,instatespace);
	//	GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,Naa); // added Naa here.
	//	GammaBranchProcess::Create(intree);
	//}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AAMutSelDPSubstitutionProcess::Delete();
		GammaBranchProcess::Delete();
	}

	GeneticCodeType codetype;
	CodonStateSpace* statespace;
	int dc;
	int fixtopo;
	int fixbl;

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};


#endif
