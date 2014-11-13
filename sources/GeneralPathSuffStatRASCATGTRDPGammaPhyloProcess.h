
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GPSSRASCATGTRDP_H
#define GPSSRASCATGTRDP_H

#include "GeneralPathSuffStatGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "GeneralPathSuffStatGTRDPProfileProcess.h"
#include "GammaBranchProcess.h"

// this is the final class implementing the CATGTR phyloprocess

class GeneralPathSuffStatRASCATGTRDPSubstitutionProcess : public virtual GeneralPathSuffStatGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual GeneralPathSuffStatGTRDPProfileProcess {

	public:

	GeneralPathSuffStatRASCATGTRDPSubstitutionProcess() {}
	virtual ~GeneralPathSuffStatRASCATGTRDPSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate, string inrrtype, int insitemin,int insitemax)	{
		GeneralPathSuffStatGTRSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		DGamRateProcess::Create(nsite,ncat);
		GeneralPathSuffStatGTRDPProfileProcess::Create(nsite,nstate);
		SetRR(inrrtype);
	}

	virtual void Delete()	{
		GeneralPathSuffStatGTRDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		GeneralPathSuffStatGTRSubstitutionProcess::Delete();
	}
};


class GeneralPathSuffStatRASCATGTRDPGammaPhyloProcess : public virtual GeneralPathSuffStatGTRPhyloProcess, public virtual GeneralPathSuffStatRASCATGTRDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();


	GeneralPathSuffStatRASCATGTRDPGammaPhyloProcess(string indatafile, string treefile, int nratecat, string inrrtype, int infixtopo, int inkappaprior, int indc, int me, int np)	{
		myid = me;
		nprocs = np;

		fixtopo = infixtopo;
		dc = indc;
		kappaprior = inkappaprior;

		datafile = indatafile;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
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

	GeneralPathSuffStatRASCATGTRDPGammaPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		is >> datafile;
		int nratecat;
		is >> nratecat;
		string inrrtype;
		is >> inrrtype;
		is >> fixtopo;
		is >> dc;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
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

	~GeneralPathSuffStatRASCATGTRDPGammaPhyloProcess() {
		Delete();
	}

	void UpdateRRSuffStat() {}
	void GlobalUpdateRRSuffStat() {}
	void SlaveUpdateRRSuffStat() {}

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
		os << "#time\ttimeperccyle\ttopo\tlnL\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

		os << ((int) (chronototal.GetTime() / 1000));
		if (chronototal.GetTime())	{
			os << '\t' << ((double) ((int) (chronototal.GetTime() / (1 + GetSize())))) / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			// os << '\t' << ((int) (chronosuffstat.GetTime() / chronototal.GetTime() * 100));
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
		os << '\n';
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
			MoveTopo(10,0);
		}
		// cerr << "gibbs ok\n";
		propchrono.Stop();

		// MPI2: reactivate this in order to test the suff stat code
		// cerr << "collapse\n";
		GlobalCollapse();
		// cerr << "collapse ok\n";

		// cerr << "branch process move\n";
		GammaBranchProcess::Move(tuning,10);
		// cerr << "branch process move ok\n";

		// cerr << "rate move\n";
		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		// cerr << "profile move\n";
		// is called inside GeneralPathSuffStatGTRDPProfileProcess::Move(1,1,10);
		// GlobalUpdateParameters();
		GeneralPathSuffStatGTRDPProfileProcess::Move(1,1,10);

		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}


		// cerr << "unfold\n";
		GlobalUnfold();
		// cerr << "unfold ok\n";

		chronototal.Stop();

		// Trace(cerr);

		return 1;
	}

	void ToStreamHeader(ostream& os)	{
		os << datafile << '\n';
		os << GetNcat() << '\n';
		os << rrtype << '\n';
		os << fixtopo << '\n';
		os << dc << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		GeneralPathSuffStatGTRDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		GeneralPathSuffStatGTRDPProfileProcess::FromStream(is);
	}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat, string inrrtype, int insitemin,int insitemax)	{
		GeneralPathSuffStatGTRPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		GeneralPathSuffStatRASCATGTRDPSubstitutionProcess::Create(indata->GetNsite(),ncat,indata->GetNstate(),inrrtype,insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatRASCATGTRDPSubstitutionProcess::Delete();
		GeneralPathSuffStatGTRPhyloProcess::Delete();
	}

	int fixtopo;
	int dc;
};

#endif

