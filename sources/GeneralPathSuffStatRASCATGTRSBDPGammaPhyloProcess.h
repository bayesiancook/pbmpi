
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSRASCATGTRSBDP_H
#define GENPATHSSRASCATGTRSBDP_H

#include "GeneralPathSuffStatGTRPhyloProcess.h"
#include "GeneralPathSuffStatGTRSubstitutionProcess.h"
#include "DGamRateProcess.h"
#include "GeneralPathSuffStatGTRSBDPProfileProcess.h"
#include "GammaBranchProcess.h"

class GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess : public virtual GeneralPathSuffStatGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual GeneralPathSuffStatGTRSBDPProfileProcess {

	public:

	GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess() {}
	virtual ~GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int nratecat, int nstate, string inrrtype, int insitemin,int insitemax)	{
		GeneralPathSuffStatGTRSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		DGamRateProcess::Create(nsite,nratecat);
		GeneralPathSuffStatGTRSBDPProfileProcess::Create(nsite,nstate);
		SetRR(inrrtype);
	}

	virtual void Delete()	{
		GeneralPathSuffStatGTRSBDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		GeneralPathSuffStatGTRSubstitutionProcess::Delete();
	}

};

class GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess : public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(string indatafile, string treefile, int nratecat, string inrrtype, int infixtopo, int inkappaprior, int indc, int me, int np)	{
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

	GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(istream& is, int me, int np)	{

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

		cerr << insitemin << '\t' << insitemax << '\n';
		cerr << nratecat << '\n';

		Create(tree,plaindata,nratecat,inrrtype,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	~GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void SlaveUpdateParameters();
	void GlobalUpdateParameters();

	virtual void SlaveExecute(MESSAGE);
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
		/*
		os << "\tprobinf";
		os << "\tstatinf";
		*/
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

		/*
		os << '\t' << GetInfProbCount();
		os << '\t' << GetStatInfCount();
		*/
		os << '\n';
	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		// cerr << "BL move\n";
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		// cerr << "gibbs\n";
		if (! fixtopo)	{
			MoveTopo(10,0);
		}
		propchrono.Stop();

		// MPI2: reactivate this in order to test the suff stat code
		// cerr << "collapse\n";
		GlobalCollapse();

		// cerr << "branch process move\n";
		GammaBranchProcess::Move(tuning,10);

		// cerr << "rates \n";
		GlobalUpdateParameters();
		// GlobalUpdateSiteRateSuffStat called inside DGamRateProcess::Move
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		// cerr << "profiles\n";
		// called inside GTRSBDPProfileProcess::Move
		// GlobalUpdateParameters();
		GeneralPathSuffStatGTRSBDPProfileProcess::Move(1,1,10);

		// GlobalUpdateParameters();

		/*
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		UpdateModeProfileSuffStat();
		*/

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
		GeneralPathSuffStatGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		GeneralPathSuffStatGTRSBDPProfileProcess::FromStream(is);
	}


	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int nratecat, string inrrtype, int insitemin,int insitemax)	{
		GeneralPathSuffStatMatrixPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Create(indata->GetNsite(),nratecat,indata->GetNstate(),inrrtype,insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
	}


	double LengthRelRateMove(double tuning, int nrep)	{

		double naccept = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogratio = - LogRRPrior() - LogLengthPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			int nbranch = MoveAllBranches(e);
			for (int i=0; i<GetNrr(); i++)	{
				rr[i] /= e;
			}
			deltalogratio += LogRRPrior() + LogLengthPrior();
			deltalogratio += (nbranch-GetNrr()) * m;

			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

			if (accepted)	{
				naccept++;
			}
			else	{
				MoveAllBranches(1.0/e);
				for (int i=0; i<GetNrr(); i++)	{
					rr[i] *= e;
				}
			}	
		}
		return naccept / nrep;
	}

	int fixtopo;
	int dc;
};

#endif

