
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AACODONMUTSELSBDPPHYLO_H
#define AACODONMUTSELSBDPPHYLO_H

//#include <cassert>
#include "AACodonMutSelGeneSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"
//#include "Parallel.h"

class AACodonMutSelGenePhyloProcess : public virtual GeneralPathSuffStatMatrixPhyloProcess, public GenePhyloProcess, public virtual AACodonMutSelGeneSubstitutionProcess, public virtual GammaBranchProcess	{

	// s'inspirer de GeneralPathSuffStatGTRPhyloProcess
	// et GeneralPathSuffStatRASCATGTRPhyloProcess

	public:

	AACodonMutSelGenePhyloProcess(string indatafile, string treefile, GeneticCodeType incodetype, int infixtopo, int infixbl, int infixcodonprofile, int infixomega, int inkappaprior, double inmintotweight, int indc, AACodonMutSelProfileProcess* inhub, int inoffset)	{
		myid = 1;
		nprocs = 2;
		fixtopo = infixtopo;
		fixbl = infixbl;
		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		dc = indc;
		kappaprior = inkappaprior;
		SetMinTotWeight(inmintotweight);

		datafile = indatafile;
		codetype = incodetype;
		SequenceAlignment* nucdata = new FileSequenceAlignment(datafile,0,myid);
		CodonSequenceAlignment* codondata = new CodonSequenceAlignment(nucdata,true,codetype);
		CodonStateSpace* statespace = codondata->GetCodonStateSpace();
		if (dc)	{
			codondata->DeleteAAConstantSites();
		}
		TaxonSet* taxonset = codondata->GetTaxonSet();

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
		
		Create(tree,codondata,statespace,fixbl,fixcodonprofile,fixomega,inhub,inoffset);
		//if (myid == 0)	{
		//	Sample();
		//	GlobalUnfold();
		//}
	}


	virtual ~AACodonMutSelGenePhyloProcess()	{
		Delete();
	}

	StateSpace* GetStateSpace()	{
		return data->GetStateSpace();
	}

	virtual void GlobalUpdateSiteProfileSuffStat() {}
	virtual void SlaveUpdateSiteProfileSuffStat() {}


	void TraceHeader(ostream& os)	{
	}

	void Trace(ostream& os)	{
	}


	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		//AACodonMutSelSBDPProfileProcess::ToStream(os);
	}
	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		//AACodonMutSelSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}


	void SaveTree(string name)	{
		ofstream os((name + datafile + ".treelist").c_str(),ios_base::app);
		GetLengthTree()->ToStream(os);
	}

	// primary scheduler

	double Move(double tuning = 1.0)	{
		//cerr << "unfold\n";
		chronototal.Start();
		propchrono.Start();
		//chronopruning.Start();
		//cerr << "bl\n";
		NonMPIBranchLengthMove(tuning);
		NonMPIBranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			int nrep = (int) (GetGibbsFactor());
			for (int rep=0; rep<nrep; rep++)	{
				NonMPIGibbsSPR();
			}
		}
		//cerr << "collapse\n";
		//chronopruning.Stop();
		propchrono.Stop();

		/*
		if (! fixcodonprofile)	{
			AACodonMutSelProfileProcess::MoveCodonProfile(tuning,30,10);
			AACodonMutSelProfileProcess::MoveNucStatCodonProfile(tuning,30,10);
			AACodonMutSelProfileProcess::MoveNucRR(tuning,2);
		}
		else {
			AACodonMutSelProfileProcess::MoveNucRR(tuning,2);
			AACodonMutSelProfileProcess::MoveNucStat(tuning,2);
		}

		if (! fixomega)	{
			AACodonMutSelProfileProcess::MoveOmega(tuning);
			AACodonMutSelProfileProcess::MoveOmega(tuning*0.1);
			AACodonMutSelProfileProcess::MoveOmega(tuning*0.01);
		}
		*/


		chronocollapse.Start();
		Collapse();
		chronocollapse.Stop();
		//cerr << "branch\n";
		GammaBranchProcess::NonMPIMove(0.1 * tuning,10);
		GammaBranchProcess::NonMPIMove(tuning,10);

		chronototal.Stop();
		//cerr << "ok\n";
		return 1;
	}


	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, CodonStateSpace* instatespace, int infixbl, int infixcodonprofile, int infixomega, AACodonMutSelProfileProcess* inhub, int inoffset)	{
		AACodonMutSelGeneSubstitutionProcess::Create(indata->GetNsite(),Naa,instatespace,infixbl,infixcodonprofile,infixomega,inhub,inoffset);
		GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,Naa,0,indata->GetNsite());  
		GammaBranchProcess::Create(intree);
	}

	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AACodonMutSelGeneSubstitutionProcess::Delete();
	}

	GeneticCodeType codetype;
	CodonStateSpace* statespace;
	int fixtopo;
	int fixbl;
	int fixcodonprofile;
	int fixomega;
	int kappaprior;
	int dc;

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};



#endif

