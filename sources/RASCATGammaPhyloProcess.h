
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCAT_H
#define RASCAT_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "DGamRateProcess.h"
#include "PoissonDPProfileProcess.h"
#include "GammaBranchProcess.h"

class RASCATSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual DGamRateProcess, public virtual PoissonDPProfileProcess {

	public:

	RASCATSubstitutionProcess() {}
	virtual ~RASCATSubstitutionProcess() {}

	protected:

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate,int insitemin,int insitemax)	{
		PoissonSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		DGamRateProcess::Create(nsite,ncat);
		PoissonDPProfileProcess::Create(nsite,nstate);
	}

	virtual void Delete()	{
		PoissonDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class RASCATGammaPhyloProcess : public virtual PoissonPhyloProcess, public virtual RASCATSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

    void GlobalSetSiteLogLCutoff();
    void SlaveSetSiteLogLCutoff();

	RASCATGammaPhyloProcess() : empcount(0) {}

	RASCATGammaPhyloProcess(string indatafile, string treefile, int nratecat, int iniscodon, GeneticCodeType incodetype, int infixtopo, int inkappaprior, double inmintotweight, int indc, int me, int np)	{
		myid = me;
		nprocs = np;

		fixtopo = infixtopo;
		dc = indc;
		kappaprior = inkappaprior;
		SetMinTotWeight(inmintotweight);

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

		Create(tree,plaindata,nratecat,insitemin,insitemax);

		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	RASCATGammaPhyloProcess(istream& is, int me, int np)	{
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
        if (atof(version.substr(0,3).c_str()) > 1.7)	{
            is >> dirweightprior;
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

		Create(tree,plaindata,nratecat,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	~RASCATGammaPhyloProcess() {
		Delete();
	}

    void GlobalSetEmpiricalPrior(istream& is);
    void SlaveSetEmpiricalPrior();

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogLikelihood()	{
		return logL;
	}

	void TraceHeader(ostream& os)	{
		os << "iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha";
		// os << "\tkappa\tallocent";
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		os << GetIndex();
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

		os << '\t' << GetLogLikelihood() << '\t' << GetRenormTotalLength() << '\t' << GetAlpha();
		os << '\t' << GetNOccupiedComponent() << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		// os << '\t' << kappa << '\t' << GetAllocEntropy();

		os << '\n';
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();
		if (! fixbl)	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}
		if (! fixtopo)	{
			MoveTopo(10,0);
		}
		propchrono.Stop();

		GlobalCollapse();

		if (! fixbl)	{
			GammaBranchProcess::Move(tuning,10);
		}

		// this one is important 
		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);
		// RASCATSubstitutionProcess::MoveRate(tuning);

		// this one is not useful
		// because uniformized process:
		// conditional on discrete substitution mapping
		// profiles do not depend on branch lengths and site rates
		// GlobalUpdateParameters();

		PoissonDPProfileProcess::Move(1,1,5);

		GlobalUnfold();
		chronototal.Stop();

		return 1;
	
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadSiteProfiles(string name, int burnin, int every, int until);
	void ReadClusters(string name, int burnin, int every, int until);
	void ReadPostHyper(string name, int burnin, int every, int until);
	void ReadSiteProfileSuffStat(string name, int burnin, int every, int until);
    /*
    void ReadISSiteLogL(string name, string empname, int burnin, int nrep);
    void SlaveComputeISSiteLogL();
    */

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << GetNcat() << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.3)	{
            os << iscodon << '\n';
            os << codetype << '\n';
            os << kappaprior << '\n';
            os << mintotweight << '\n';
            if (atof(version.substr(0,3).c_str()) > 1.7)	{
                os << dirweightprior << '\n';
            }
		}
		os << fixtopo << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
            os << NSPR << '\t' << NNNI << '\n';
        }
		os << dc << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		PoissonDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		PoissonDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}


	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat,int insitemin,int insitemax)	{
		PoissonPhyloProcess::Create(intree,indata);
		// PoissonPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		RASCATSubstitutionProcess::Create(indata->GetNsite(),ncat,indata->GetNstate(),insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}

	int iscodon;
	GeneticCodeType codetype;
    double siteloglcutoff;
    double* empcount;
};

#endif

