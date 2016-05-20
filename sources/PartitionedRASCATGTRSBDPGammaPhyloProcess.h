
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTRASCATGTRSBDP_H
#define PARTRASCATGTRSBDP_H

#include "PartitionedDGamRateProcess.h"
#include "PartitionedExpoConjugateGTRSBDPProfileProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"
#include "PartitionedExpoConjugateGTRGammaPhyloProcess.h"


class PartitionedRASCATGTRSBDPGammaPhyloProcess : public virtual PartitionedExpoConjugateGTRSBDPProfileProcess, public virtual PartitionedExpoConjugateGTRGammaPhyloProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();


	PartitionedRASCATGTRSBDPGammaPhyloProcess(string indatafile, string treefile, string inschemefile, bool inlinkgam,bool inunlinkgtr,bool inlinkmult,string inrrtype,int nratecat, int iniscodon, GeneticCodeType incodetype, int infixtopo, int inNSPR, int inNNNI, int inkappaprior, double inmintotweight, int me, int np)	{
		partoccupancy = 0;
		occupancyNeedsUpdating = true;

		myid = me;
		nprocs = np;

		fixtopo = infixtopo;
		NSPR = inNSPR;
		NNNI = inNNNI;
		iscodon = iniscodon;
		codetype = incodetype;
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

		Create(tree,plaindata,nratecat,inlinkmult,schemes[0],schemes[2],insitemin,insitemax);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	PartitionedRASCATGTRSBDPGammaPhyloProcess(istream& is, int me, int np)	{
		partoccupancy = 0;
		occupancyNeedsUpdating = true;

		myid = me;
		nprocs = np;

		FromStreamHeader(is);

		int insitemin = -1,insitemax = -1;
		if (myid > 0) {
			int width = data->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = data->GetNsite();
			}
			else {
				insitemax = myid*width;
			}
		}

		vector<PartitionScheme> schemes = PartitionedDGamRateProcess::ReadSchemes(schemefile, data->GetNsite(), myid, linkgam, unlinkgtr, rrtype);

		Create(tree,data,PartitionedDGamRateProcess::Ncat,linkmult,schemes[0],schemes[2],insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}

	}

	~PartitionedRASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "iter\ttime\ttopo\tloglik\tlength";
		if(PartitionedDGamRateProcess::GetNpart() > 1)
			os << "\tmeanalpha";
		else
			os << "\talpha";
		os << "\tNmode\tstatent\tstatalpha";

		if(PartitionedDGamRateProcess::GetNpart() > 1)
		{
			os << "\talphahyper";
			if(!LinkedMultipliers())
				os << "\tmultent\tmultalpha";
		}

		if (nfreerr > 0)
		{
			os << "\trrent\trrmean";
		}
		// os << "\tkappa\tallocent";
		os << endl;
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

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
		os << "\t" << GetNDisplayedComponent();
		os << "\t" << GetStatEnt();
		os << "\t" << GetMeanDirWeight();

		if(PartitionedDGamRateProcess::GetNpart() > 1)
		{
			os << "\t" << GetAlphaHyper();
			if(!LinkedMultipliers())
			{
				os << "\t" << GetMultiplierEntropy();
				os << "\t" << GetMultHyper();
			}
		}

		if (nfreerr > 0)
		{
			os << "\t" << GetRREntropy();
			os << "\t" << GetRRMean();
		}
		// os << '\t' << kappa << '\t' << GetAllocEntropy();
		os << endl;

	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
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
		PartitionedDGamRateProcess::Move(0.3*tuning,15);
		PartitionedDGamRateProcess::Move(0.03*tuning,15);

		// is called inside ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		// GlobalUpdateParameters();
		PartitionedExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		if (iscodon){
			PartitionedExpoConjugateGTRSBDPProfileProcess::Move(0.1,1,15);
			PartitionedExpoConjugateGTRSBDPProfileProcess::Move(0.01,1,15);
		}

		occupancyNeedsUpdating = true;

		if (PartitionedGTRProfileProcess::GetNpart() == nfreerr){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		if(PartitionedDGamRateProcess::GetNpart() > 1 && !LinkedMultipliers())
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
		os << PartitionedDGamRateProcess::LinkedMultipliers() << '\n';
		os << rrtype << '\n';
		os << iscodon << '\n';
		os << codetype;
		os << kappaprior << '\n';
		os << mintotweight << '\n';
		os << fixtopo << '\n';
		os << NSPR << '\t' << NNNI << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}
	void FromStreamHeader(istream& is)	{
		PhyloProcess::FromStreamHeader(is);
		is >> datafile;
		is >> schemefile;
		is >> linkgam;
		is >> unlinkgtr;
		is >> PartitionedDGamRateProcess::Ncat;
		is >> linkmult;
		is >> rrtype;
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
		is >> fixtopo;
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
			is >> NSPR;
			is >> NNNI;
		}
		else	{
			NSPR = 10;
			NNNI = 0;
		}
		if (iscodon)	{
			SequenceAlignment* tempdata = new FileSequenceAlignment(datafile,0,myid,false);
			data = new CodonSequenceAlignment(tempdata,true,codetype);
		}
		else	{
			data = new FileSequenceAlignment(datafile,0,myid,false);
		}
		const TaxonSet* taxonset = data->GetTaxonSet();
		tree = new Tree(taxonset);
		tree->ReadFromStream(is);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		PartitionedDGamRateProcess::ToStream(os);
		PartitionedExpoConjugateGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		PartitionedDGamRateProcess::FromStream(is);
		PartitionedExpoConjugateGTRSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadNocc(string name, int burnin, int every, int until);
	void ReadRelRates(string name, int burnin, int every, int until);
	void ReadSiteProfiles(string name, int burnin, int every, int until);
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat, bool inlinkmult, PartitionScheme rrscheme, PartitionScheme dgamscheme, int insitemin,int insitemax);
		
	virtual void Delete();


	// Importantly, this assumes that DGam partitions are always sub-partitions of GTR partitions
	double GetNormalizationFactor();
	double GetNormPartRate(int d, int p);

	void UpdatePartOccupancyNumbers();

	int NSPR;
	int NNNI;
	int iscodon;
	GeneticCodeType codetype;

	string schemefile;
	bool linkgam;
	bool unlinkgtr;
	string rrtype;

	double normFactor;
	bool occupancyNeedsUpdating;

	int*** partoccupancy;
};

#endif

