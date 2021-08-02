
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
#include "AACodonMutSelSBDPSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"
//#include "Parallel.h"

class AACodonMutSelSBDPPhyloProcess : public virtual AACodonMutSelSBDPSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	public:

	AACodonMutSelSBDPPhyloProcess(string indatafile, string treefile, GeneticCodeType incodetype, int infixtopo, int infixbl, int inNSPR, int inNNNI, int infixcodonprofile, int infixomega, int inomegaprior, int inkappaprior, int indirweightprior, double inmintotweight, int indc, int me, int np)	{
		myid = me;
		nprocs = np;
		fixtopo = infixtopo;
		fixbl = infixbl;
		NSPR = inNSPR;
		NNNI = inNNNI;
		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		omegaprior = inomegaprior;
		dc = indc;
		kappaprior = inkappaprior;
		dirweightprior = indirweightprior;
		SetMinTotWeight(inmintotweight);

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

		Create(tree,codondata,insitemin,insitemax,statespace,fixcodonprofile,fixomega);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	AACodonMutSelSBDPPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		is >> codetype;
		if (atof(version.substr(0,3).c_str()) > 1.8)	{
            int nmax;
            is >> nmax;
            SetNmodeMax(nmax);
        }
		is >> fixtopo;
		is >> fixbl;
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
			is >> NSPR;
			is >> NNNI;
		}
		else	{
			NSPR = 10;
			NNNI = 0;
		}
		is >> fixcodonprofile;
		is >> fixomega;
		if (atof(version.substr(0,3).c_str()) > 1.5)	{
			is >> omegaprior;
		}
		else	{
			omegaprior = 0;
		}
		is >> kappaprior;
		is >> dirweightprior;
		is >> mintotweight;
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

		Create(tree,codondata,insitemin,insitemax,statespace,fixcodonprofile,fixomega);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	virtual ~AACodonMutSelSBDPPhyloProcess()	{
		Delete();
	}

	// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
	// are model dependent
	// should be implemented in .cpp file
        virtual void SlaveExecute(MESSAGE);
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();
	
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

    void GlobalSetSiteLogLCutoff();
    void SlaveSetSiteLogLCutoff();

    void GlobalSetTestData();
    void SlaveSetTestData();

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogLikelihood()	{
		return logL;
	}

	void TraceHeader(ostream& os)	{
		os << "iter\ttime\tpruning\tlnL\tlength\tcodonent\tomega\tNmode\tstatent\tstatalpha\tnucsA\tnucsC\tnucsG\tnucsT\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT";
		os << '\n'; 
		//os << "lnL\tlength\tNmode\tNocc\tnucsA\tnucsC\tnucsT\tnucsG\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT\tstatent";
		//os << "\ttotaltime";
		//os << "\tpruning\tsuffstat\tunfold\tcollapse";
		//os << "\n";
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

		//os << ((int) (chronototal.GetTime() / 1000));
		os << GetIndex();
		if (chronototal.GetTime())	{
			os << '\t' << chronototal.GetTime() / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
			//os << '\t' << ((double) ((int) (chronototal.GetTime() / (GetSize())))) / 1000;
			//os << '\t' << ((int) (100 * chronopruning.GetTime() /chronototal.GetTime()));
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' <<  GetLogLikelihood();
		os << '\t' << GetTotalLength();
		//os << '\t' << GetAlpha();
		os << '\t' << GetCodonProfileEntropy();
		os << '\t' << GetOmega();
		os << '\t' << GetNDisplayedComponent();
		//os << '\t' << GetNOccupiedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << GetNucStat(0) << '\t' << GetNucStat(1) << '\t' << GetNucStat(2) << '\t' << GetNucStat(3);
		os << '\t' << GetNucRR(0) << '\t' << GetNucRR(1) << '\t' << GetNucRR(2) << '\t' << GetNucRR(3) << '\t' << GetNucRR(4) << '\t' << GetNucRR(5);

		os << '\n';

		/*
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
		*/
	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << codetype << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.8)	{
            os << GetNmodeMax() << '\n';
        }
		os << fixtopo << '\n';
		os << fixbl << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
            os << NSPR << '\t' << NNNI << '\n';
        }
		os << fixcodonprofile << '\n';
		os << fixomega << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.5)	{
            os << omegaprior << '\n';
        }
		os << kappaprior << '\n';
		os << dirweightprior << '\n';
		os << mintotweight << '\n';
		os << dc << '\n';
		GetTree()->ToStream(os);
	}
	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		AACodonMutSelSBDPProfileProcess::ToStream(os);
	}
	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AACodonMutSelSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void Read(string name, int burnin, int every, int until);
	void ReadMapStats(string name, int burnin, int every, int until);
	int CountNonSynMapping(int i);
	int CountNonSynMapping();
	int GlobalNonSynMapping();
	virtual void SlaveNonSynMapping();

	// primary scheduler

	double Move(double tuning = 1.0)	{
		chronototal.Start();
		propchrono.Start();
		if (! fixbl)	{
			BranchLengthMove(0.1 * tuning);
			BranchLengthMove(tuning);
		}
		if (! fixtopo)	{
			//GibbsSPR(50);
			MoveTopo(NSPR,NNNI);
		}
		propchrono.Stop();

		chronosuffstat.Start();

		chronocollapse.Start();
		GlobalCollapse();
		chronocollapse.Stop();
		if (! fixbl)	{
			GammaBranchProcess::Move(0.1 * tuning,10);
			GammaBranchProcess::Move(tuning,10);
		}

		GlobalUpdateParameters();
		AACodonMutSelSBDPProfileProcess::Move(tuning,1,15);
		chronosuffstat.Stop();

		chronounfold.Start();
		GlobalUnfold();
		chronounfold.Stop();

		chronototal.Stop();
		return 1;
	}


	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int sitemin, int sitemax)	{
		cerr << "In two-argument Create of AACodonMutSelSBDPPhyloProcess.  Should not be here.\n";
		exit(1);
	}

	virtual void Create(Tree* intree, SequenceAlignment* indata, int sitemin, int sitemax, CodonStateSpace* instatespace, int infixcodonprofile, int infixomega)	{
		//cerr << "just before creating substitution process\n";
		//cerr.flush();
		AACodonMutSelSBDPSubstitutionProcess::Create(indata->GetNsite(),Naa,sitemin,sitemax,instatespace,infixcodonprofile,infixomega);
		//cerr << "just before creating phyloprocess\n";
		//cerr.flush();
		GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,Naa,sitemin,sitemax); 
		//cerr << "just before creating branch process\n";
		//cerr.flush();
		GammaBranchProcess::Create(intree);
		//cerr << "Done create\n";
		//cerr.flush();
	}

	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AACodonMutSelSBDPSubstitutionProcess::Delete();
	}

	GeneticCodeType codetype;
	CodonStateSpace* statespace;
	int fixcodonprofile;
	int fixomega;
    double siteloglcutoff;

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};



#endif

