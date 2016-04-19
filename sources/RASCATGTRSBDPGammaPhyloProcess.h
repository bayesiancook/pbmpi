
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATGTRSBDP_H
#define RASCATGTRSBDP_H

#include "ExpoConjugateGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "ExpoConjugateGTRSBDPProfileProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"

// this is the final class implementing the CATGTR phyloprocess

class RASCATGTRSBDPSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRSBDPProfileProcess {

	public:

	RASCATGTRSBDPSubstitutionProcess() {}
	virtual ~RASCATGTRSBDPSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate, string inrrtype, int insitemin,int insitemax)	{
		ExpoConjugateGTRSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		DGamRateProcess::Create(nsite,ncat);
		ExpoConjugateGTRSBDPProfileProcess::Create(nsite,nstate);
		SetRR(inrrtype);
	}

	virtual void Delete()	{
		ExpoConjugateGTRSBDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};


class RASCATGTRSBDPGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual RASCATGTRSBDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	using MixtureProfileProcess::LogStatPrior;

        virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();


	RASCATGTRSBDPGammaPhyloProcess(string indatafile, string treefile, int nratecat, int iniscodon, GeneticCodeType incodetype, string inrrtype, int infixtopo, int inNSPR, int inNNNI, int inkappaprior, double inmintotweight, int indc, int incinit, int me, int np)	{
		myid = me;
		nprocs = np;

		InitIncremental = incinit;

		fixtopo = infixtopo;
		NSPR = inNSPR;
		NNNI = inNNNI;
		iscodon = iniscodon;
		codetype = incodetype;
		dc = indc;
		kappaprior = inkappaprior;
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
		if (dc)	{
			plaindata->DeleteConstantSites();
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

		Create(tree,plaindata,nratecat,inrrtype,insitemin,insitemax);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	RASCATGTRSBDPGammaPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		int nratecat;
		is >> nratecat;
		if (atof(version.substr(0,3).c_str()) > 1.3)	{
			is >> iscodon;
			is >> codetype;
			is >> kappaprior;
			is >> mintotweight;
		}
		else	{
			iscodon = 0;
			codetype = Universal;
			kappaprior = 0;
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
		is >> dc;
		//SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
		SequenceAlignment* plaindata;
		if (iscodon)	{
			SequenceAlignment* tempdata = new FileSequenceAlignment(datafile,0,myid);
			plaindata = new CodonSequenceAlignment(tempdata,true,codetype);
		}
		else	{
			plaindata = new FileSequenceAlignment(datafile,0,myid);
		}
		if (dc)	{
			plaindata->DeleteConstantSites();
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

		Create(tree,plaindata,nratecat,inrrtype,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}

	}

	~RASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		// os << "\tkappa\tallocent";
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

		os << GetSize();
		if (chronototal.GetTime())	{
			os << '\t' << chronototal.GetTime() / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood();
		os << '\t' << GetRenormTotalLength();
		os << '\t' << GetAlpha();
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}
		// os << '\t' << kappa << '\t' << GetAllocEntropy();
		os << '\n';

	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
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
		// raph
		/*
		if (! fixtopo)	{
			cout << "SPR " << GibbsSPR(4)<<'\n';	
			for(int i=0; i<3; i++){
				cout << "NNI"<<i<<' '<< GibbsNNI(0.1,1)<<'\n';
			}

		}
		*/

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
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		// cerr << "profile move\n";
		// is called inside ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		// GlobalUpdateParameters();
		GlobalUpdateParameters();
		ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		if (iscodon){
			ExpoConjugateGTRSBDPProfileProcess::Move(0.1,1,15);
			ExpoConjugateGTRSBDPProfileProcess::Move(0.01,1,15);
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
		os << GetNcat() << '\n';
		os << iscodon << '\n';
		os << codetype << '\n';
		os << kappaprior << '\n';
		os << mintotweight << '\n';
		os << rrtype << '\n';
		os << fixtopo << '\n';
		os << NSPR << '\t' << NNNI << '\n';
		os << dc << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		ExpoConjugateGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		ExpoConjugateGTRSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
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

	virtual void ReadPB(int argc, char* argv[]);
	void ReadNocc(string name, int burnin, int every, int until);
	void ReadRelRates(string name, int burnin, int every, int until);
	void ReadSiteProfiles(string name, int burnin, int every, int until);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	double ProfileProposal(double* profile, int cat, double tuning);
	void ReadTestProfile(string name, int nrep, double tuning, int burnin, int every, int until);

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat, string inrrtype, int insitemin,int insitemax)	{
		ExpoConjugateGTRPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		RASCATGTRSBDPSubstitutionProcess::Create(indata->GetNsite(),ncat,indata->GetNstate(),inrrtype,insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATGTRSBDPSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

	int fixtopo;
	int NSPR;
	int NNNI;
	int iscodon;
	GeneticCodeType codetype;
	int dc;
};

#endif

