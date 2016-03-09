
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTRASCATGTRFINITE_H
#define PARTRASCATGTRFINITE_H

#include "PartitionedDGamRateProcess.h"
#include "PartitionedExpoConjugateGTRFiniteProfileProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"
#include "PartitionedExpoConjugateGTRGammaPhyloProcess.h"


class PartitionedRASCATGTRFiniteGammaPhyloProcess : public virtual PartitionedExpoConjugateGTRFiniteProfileProcess, public virtual PartitionedExpoConjugateGTRGammaPhyloProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();


	PartitionedRASCATGTRFiniteGammaPhyloProcess(string indatafile, string treefile, string inschemefile, int nratecat, int ncat, int infixncomp, int inempmix, string inmixtype, int infixtopo, int inNSPR, int inNNNI, int me, int np)	{
		partoccupancy = 0;

		myid = me;
		nprocs = np;

		fixtopo = infixtopo;
		NSPR = inNSPR;
		NNNI = inNNNI;

		datafile = indatafile;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);

		// data = new FileSequenceAlignment(datafile,0,myid);
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
		vector<PartitionScheme> schemes = PartitionedDGamRateProcess::ReadSchemes(schemefile, plaindata->GetNsite());

		if(myid == 0)
		{
			cerr << endl;
			cerr << "Read " << schemes[2].Npart << " partitions in scheme file '" << schemefile << "':\n";
			for(size_t i = 0; i < schemes[0].Npart; i++)
			{
				string t = schemes[0].partType[i] == "None" ? "GTR" : schemes[0].partType[i];

				cerr << t << "\t" << schemes[0].partSites[i].size() << " sites" << endl;
			}
			cerr << endl;
		}
		Create(tree,plaindata,nratecat,ncat,schemes[0],schemes[2],infixncomp,inempmix,inmixtype,insitemin,insitemax);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	PartitionedRASCATGTRFiniteGammaPhyloProcess(istream& is, int me, int np)	{
		partoccupancy = 0;

		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		is >> schemefile;
		int nratecat;
		is >> nratecat;
		int infixncomp;
		int inempmix;
		string inmixtype;
		string inrrtype;
		is >> infixncomp >> inempmix >> inmixtype;
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
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
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

		vector<PartitionScheme> schemes = PartitionedDGamRateProcess::ReadSchemes(schemefile, plaindata->GetNsite());

		Create(tree,plaindata,nratecat,1,schemes[0],schemes[2],infixncomp,inempmix,inmixtype,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}

	}

	~PartitionedRASCATGTRFiniteGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& hs)	{
		stringstream os;
		os << "iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha\tmultent\tmultalpha";
		if (nfreerr > 0)
		{
			os << "\trrent\trrmean";
		}
		// os << "\tkappa\tallocent";
		os << endl;
		hs << os.str();
	}

	void Trace(ostream& hs)	{

		UpdateOccupancyNumbers();

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
		os << "\t" << GetNDisplayedComponent();
		os << "\t" << GetStatEnt();
		os << "\t" << GetMeanDirWeight();
		os << "\t" << GetMultiplierEntropy();
		os << "\t" << GetMultHyper();
		if (nfreerr > 0)
		{
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
		PartitionedExpoConjugateGTRFiniteProfileProcess::Move(1,1,10);

		if (PartitionedGTRProfileProcess::GetNpart() == nfreerr){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}
		else
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
		os << GetNcat() << '\n';
		os << fixncomp << '\t' << empmix << '\t' << mixtype << '\n';
		os << fixtopo << '\n';
		os << NSPR << '\t' << NNNI << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		PartitionedDGamRateProcess::ToStream(os);
		PartitionedExpoConjugateGTRFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		PartitionedDGamRateProcess::FromStream(is);
		PartitionedExpoConjugateGTRFiniteProfileProcess::FromStream(is);
		//GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int nratecat,int ncat,PartitionScheme rrscheme, PartitionScheme dgamscheme, int infixncomp, int inempmix, string inmixtype, int insitemin,int insitemax);
		
	virtual void Delete();


	// Importantly, this assumes that DGam partitions are always sub-partitions of GTR partitions
	double GetNormalizationFactor();
	double GetNormPartRate(int d, int p);

	void UpdatePartOccupancyNumbers();

	int fixtopo;
	int NSPR;
	int NNNI;

	string schemefile;

	int*** partoccupancy;
};

#endif

