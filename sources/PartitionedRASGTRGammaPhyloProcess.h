
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

#include "PartitionedDGamRateProcess.h"
#include "PartitionedExpoConjugateGTRPartitionedProfileProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"
#include "PartitionedExpoConjugateGTRGammaPhyloProcess.h"


class PartitionedRASGTRGammaPhyloProcess : public virtual PartitionedExpoConjugateGTRPartitionedProfileProcess, public virtual PartitionedExpoConjugateGTRGammaPhyloProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();


	PartitionedRASGTRGammaPhyloProcess(string indatafile, string treefile, string inschemefile,bool inlinkgam,bool inunlinkgtr,string inrrtype,int nratecat, int iniscodon, GeneticCodeType incodetype, int infixtopo, int inNSPR, int inNNNI, double inmintotweight, int me, int np)	{
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

		schemefile = inschemefile;
		linkgam = inlinkgam;
		unlinkgtr = inunlinkgtr;
		rrtype = inrrtype;
		vector<PartitionScheme> schemes = PartitionedDGamRateProcess::ReadSchemes(schemefile, plaindata->GetNsite(), myid, linkgam, unlinkgtr, rrtype);

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
		is >> linkgam;
		is >> unlinkgtr;
		int nratecat;
		is >> nratecat;
		is >> rrtype;
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
		is >> fixtopo;
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
			is >> NSPR;
			is >> NNNI;
		}
		else	{
			NSPR = 10;
			NNNI = 0;
		}
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

		vector<PartitionScheme> schemes = PartitionedDGamRateProcess::ReadSchemes(schemefile, plaindata->GetNsite(), myid, linkgam, unlinkgtr, rrtype);

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

		if(PartitionedDGamRateProcess::GetNpart() > 1)
			os << "\tmultent\tmultalpha";

		if (nfreestat > 0)	{
			os << "\tstatent";
			if(nfreestat > 1)
				os << "\tstatalpha";
		}
		if (nfreerr > 0)
		{
			os << "\trrent\trrmean";
		}
		os << endl;
	}

	void Trace(ostream& hs)	{
		stringstream os;
		os << GetSize() - 1;
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

		if(PartitionedDGamRateProcess::GetNpart() > 1)
		{
			os << "\t" << GetMultiplierEntropy();
			os << "\t" << GetMultHyper();
		}

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
		os << endl;
		hs << os.str();

	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo(NSPR,NNNI);
		}

		propchrono.Stop();

		
		// MPI2: reactivate this in order to test the suff stat code
		// chronocollapse.Start();
		GlobalCollapse();
		// chronocollapse.Stop();

		// chronosuffstat.Start();
		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();
		PartitionedDGamRateProcess::Move(0.3*tuning,10);
		PartitionedDGamRateProcess::Move(0.03*tuning,10);

		// is called inside ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		// GlobalUpdateParameters();
		PartitionedExpoConjugateGTRPartitionedProfileProcess::Move(1,1,10);
		if (iscodon){
			PartitionedExpoConjugateGTRPartitionedProfileProcess::Move(0.1,1,15);
			PartitionedExpoConjugateGTRPartitionedProfileProcess::Move(0.01,1,15);
		}

		if (PartitionedGTRProfileProcess::GetNpart() == nfreerr){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		if(PartitionedDGamRateProcess::GetNpart() > 1)
		{
			LengthMultiplierMove(1,10);
			LengthMultiplierMove(0.1,10);
			LengthMultiplierMove(0.01,10);

			if(nfreerr > 0)
			{
				MultiplierRelRateMove(1,10);
				MultiplierRelRateMove(0.1,10);
				MultiplierRelRateMove(0.01,10);
			}
		}

		// chronosuffstat.Stop();

		// chronounfold.Start();
		GlobalUnfold();
		// chronounfold.Stop();

		chronototal.Stop();

		// Trace(cerr);

		return 1;
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << schemefile << '\n';
		os << linkgam << '\n';
		os << unlinkgtr << '\n';
		os << GetNcat() << '\n';
		os << rrtype << '\n';
		os << iscodon << '\n';
		os << codetype;
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

	// Importantly, this assumes that DGam partitions are always sub-partitions of GTR partitions
	virtual double GetNormalizationFactor()	{
		double total = 0;
		for (int dgampart=0; dgampart<PartitionedDGamRateProcess::GetNpart(); dgampart++)
		{
			vector<int> partsites = PartitionedDGamRateProcess::GetPartSites(dgampart);

			size_t rrpart = PartitionedGTRProfileProcess::GetSitePart(partsites.front());
			size_t statpart = PartitionedProfileProcess::GetSitePart(partsites.front());

			total += PartitionedGTRPartitionedProfileProcess::GetNormRate(rrpart,statpart) * PartitionedDGamRateProcess::GetRateMultiplier(dgampart) * partsites.size();
		}
		total /= GetNsite();
		return total;
	}

	virtual void ReadPB(int argc, char* argv[]);
	void GlobalSetTestData();
	void SlaveSetTestData();
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat, PartitionScheme rrscheme, PartitionScheme statscheme, PartitionScheme dgamscheme, int insitemin,int insitemax)	{
		PartitionedExpoConjugateGTRGammaPhyloProcess::Create(intree,indata,indata->GetNstate(),rrscheme,ncat,dgamscheme,insitemin,insitemax);
		PartitionedExpoConjugateGTRPartitionedProfileProcess::Create(indata->GetNstate(), rrscheme, statscheme);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		PartitionedExpoConjugateGTRPartitionedProfileProcess::Delete();
		PartitionedExpoConjugateGTRGammaPhyloProcess::Delete();
	}

	int fixtopo;
	int NSPR;
	int NNNI;
	int iscodon;
	GeneticCodeType codetype;

	string schemefile;
	string cvschemefile;
	bool linkgam;
	bool unlinkgtr;
	string rrtype;
};

#endif

