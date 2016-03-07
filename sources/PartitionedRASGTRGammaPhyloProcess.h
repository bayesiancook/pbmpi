
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTRASGTR_H
#define PARTRASGTR_H

#include "PartitionedExpoConjugateGTRPhyloProcess.h"
#include "PartitionedDGamRateProcess.h"
#include "PartitionedExpoConjugateGTRPartitionedProfileProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"

// this is the final class implementing the PartitionedGTRPartitionedProfile phyloprocess

class PartitionedRASGTRSubstitutionProcess : public virtual PartitionedExpoConjugateGTRSubstitutionProcess, public virtual PartitionedDGamRateProcess, public virtual PartitionedExpoConjugateGTRPartitionedProfileProcess {

	public:

	PartitionedRASGTRSubstitutionProcess() {}
	virtual ~PartitionedRASGTRSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in PartitionedRASGTRSubstitutionProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int indim, int ncat, PartitionScheme rrscheme, PartitionScheme statscheme, PartitionScheme dgamscheme, int insitemin,int insitemax)	{
		PartitionedExpoConjugateGTRSubstitutionProcess::Create(indim, rrscheme, insitemin, insitemax);
		PartitionedDGamRateProcess::Create(ncat, dgamscheme);
		PartitionedExpoConjugateGTRPartitionedProfileProcess::Create(indim, rrscheme, statscheme);
	}

	virtual void Delete()	{
		PartitionedExpoConjugateGTRPartitionedProfileProcess::Delete();
		PartitionedDGamRateProcess::Delete();
		PartitionedExpoConjugateGTRSubstitutionProcess::Delete();
	}
};


class PartitionedRASGTRGammaPhyloProcess : public virtual PartitionedExpoConjugateGTRPhyloProcess, public virtual PartitionedRASGTRSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();


	PartitionedRASGTRGammaPhyloProcess(string indatafile, string treefile, string partfile, int nratecat, int iniscodon, GeneticCodeType incodetype, int innrrpart, string* inrrtype, int innstatpart, string* instattype, int infixtopo, int inNSPR, int inNNNI, double inmintotweight, int me, int np)	{
		myid = me;
		nprocs = np;

		fixtopo = infixtopo;
		NSPR = inNSPR;
		NNNI = inNNNI;
		iscodon = iniscodon;
		codetype = incodetype;
		SetMinTotWeight(inmintotweight);

		datafile = indatafile;
		SequenceAlignment* plaindata;
		if (iscodon)	{
			SequenceAlignment* tempdata = new FileSequenceAlignment(datafile,0,myid);
			plaindata = new CodonSequenceAlignment(tempdata,true,codetype);
		}
		else	{
			plaindata = new FileSequenceAlignment(datafile,0,myid);
		}
		const TaxonSet* taxonset = plaindata->GetTaxonSet();
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
			int width = plaindata->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = plaindata->GetNsite();
			}
			else {
				insitemax = myid*width;
			}
		}

		vector<PartitionScheme> schemes = ReadSchemes(schemefile);

		Create(tree,plaindata,nratecat,schemes[0],schemes[1],schemes[2],insitemin,insitemax);

		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	PartitionedRASGTRGammaPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		is >> schemefile;
		int nratecat;
		is >> nratecat;
		if (atof(version.substr(0,3).c_str()) > 1.3)	{
			is >> iscodon;
			is >> codetype;
			is >> mintotweight;
		}
		else	{
			iscodon = 0;
			codetype = Universal;
			mintotweight = -1;
		}
		string inrrtype;
		is >> inrrtype;
		is >> fixtopo;
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
			is >> NSPR;
			is >> NNNI;
		}
		else	{
			NSPR = 10;
			NNNI = 0;
		}
		//SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
		SequenceAlignment* plaindata;
		if (iscodon)	{
			SequenceAlignment* tempdata = new FileSequenceAlignment(datafile,0,myid);
			plaindata = new CodonSequenceAlignment(tempdata,true,codetype);
		}
		else	{
			plaindata = new FileSequenceAlignment(datafile,0,myid);
		}
		const TaxonSet* taxonset = plaindata->GetTaxonSet();

		int insitemin = -1,insitemax = -1;
		if (myid > 0) {
			int width = plaindata->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = plaindata->GetNsite();
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

		vector<PartitionScheme> schemes = ReadSchemes(schemefile);

		Create(tree,plaindata,nratecat,schemes[0],schemes[1],schemes[2],insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}

	}

	~PartitionedRASGTRGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "iter\ttime\ttopo\tloglik\tlength\talpha";

		if (nfreestat > 0)	{
			os << "\tstatent";
			if(nfreestat > 1)
				os << "\tstatalpha";
		}
		if (nfreerr > 0)
		{
			os << "\trrent\trrmean";
		}
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		os << GetSize();
		if (chronototal.GetTime())	{
			os << "\t" << chronototal.GetTime() / 1000;
			os << "\t" << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
		}
		else	{
			os << "\t" << 0;
			os << "\t" << 0;
		}

		os << "\t" << GetLogLikelihood();
		os << "\t" << GetRenormTotalLength();
		os << "\t" << GetAlpha();
		if (nfreestat > 0)
		{
			os << "\t" << GetStatEnt();

			if (nfreestat > 1)
				os << "\t" << GetMeanDirWeight();
		}
		if (nfreerr > 0)	{
			os << "\t" << GetRREntropy();
			os << "\t" << GetRRMean();
		}
		// os << '\t' << kappa << '\t' << GetAllocEntropy();
		os << "\n";

	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		// cerr << "BL move\n";
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		// cerr << "BL move ok\n";
		// cerr << "gibbs\n";
		if (! fixtopo)	{
			MoveTopo(NSPR,NNNI);
		}

		// cerr << "gibbs ok\n";
		propchrono.Stop();

		
		// MPI2: reactivate this in order to test the suff stat code
		// chronocollapse.Start();
		// cerr << "collapse\n";
		GlobalCollapse();
		// cerr << "collapse ok\n";
		// chronocollapse.Stop();

		// chronosuffstat.Start();
		// cerr << "branch process move\n";
		GammaBranchProcess::Move(tuning,10);
		// cerr << "branch process move ok\n";

		// cerr << "rate move\n";
		GlobalUpdateParameters();
		PartitionedDGamRateProcess::Move(0.3*tuning,10);
		PartitionedDGamRateProcess::Move(0.03*tuning,10);

		// cerr << "profile move\n";
		// is called inside ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		// GlobalUpdateParameters();
		PartitionedExpoConjugateGTRPartitionedProfileProcess::Move(1,1,10);
		if (iscodon){
			PartitionedExpoConjugateGTRPartitionedProfileProcess::Move(0.1,1,15);
			PartitionedExpoConjugateGTRPartitionedProfileProcess::Move(0.01,1,15);
		}

		if (! fixrr){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		// chronosuffstat.Stop();

		// chronounfold.Start();
		// cerr << "unfold\n";
		GlobalUnfold();
		// cerr << "unfold ok\n";
		// chronounfold.Stop();

		chronototal.Stop();

		// Trace(cerr);

		return 1;
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << schemefile << '\n';
		os << GetNcat() << '\n';
		os << iscodon << '\n';
		os << codetype << '\n';
		os << mintotweight << '\n';
		os << fixtopo << '\n';
		os << NSPR << '\t' << NNNI << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		PartitionedDGamRateProcess::ToStream(os);
		PartitionedExpoConjugateGTRPartitionedProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		PartitionedDGamRateProcess::FromStream(is);
		PartitionedExpoConjugateGTRPartitionedProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	vector<PartitionScheme> ReadSchemes(string partfile);

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat, PartitionScheme rrscheme, PartitionScheme statscheme, PartitionScheme dgamscheme, int insitemin,int insitemax)	{
		PartitionedExpoConjugateGTRPhyloProcess::Create(intree,indata,indata->GetNstate(),rrscheme,insitemin,insitemax);
		PartitionedRASGTRSubstitutionProcess::Create(indata->GetNstate(), ncat, rrscheme, statscheme, dgamscheme,insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		PartitionedRASGTRSubstitutionProcess::Delete();
		PartitionedExpoConjugateGTRPhyloProcess::Delete();
	}

	int fixtopo;
	int NSPR;
	int NNNI;
	int iscodon;
	GeneticCodeType codetype;

	string schemefile;
};

#endif

