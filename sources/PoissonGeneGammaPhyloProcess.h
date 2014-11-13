
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef POISSONGENE_H
#define POISSONGENE_H

#include "GeneProfileProcess.h"

#include "PoissonPhyloProcess.h"
#include "GammaBranchProcess.h"
#include "DGamRateProcess.h"

class PoissonGeneProfileProcess : public virtual PoissonProfileProcess, public virtual GeneProfileProcess	{

	public:

	PoissonGeneProfileProcess() : poissonhub(0) {}
	virtual ~PoissonGeneProfileProcess() {}

	protected:

	virtual void Create(int innsite, int innstate, PoissonProfileProcess* inhub, int inoffset)	{
		PoissonProfileProcess::Create(innsite,innstate);
		GeneProfileProcess::Create(innsite,innstate,inhub,inoffset);
		poissonhub = inhub;
	}

	virtual void Delete()	{
		GeneProfileProcess::Delete();
		PoissonProfileProcess::Delete();
	}

	PoissonProfileProcess* poissonhub;

};


class PoissonGeneSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual DGamRateProcess, public virtual PoissonGeneProfileProcess {

	public:

	PoissonGeneSubstitutionProcess() {}
	virtual ~PoissonGeneSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in Poisson Gene Sub Process::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate, PoissonProfileProcess* hub, int offset)	{
		PoissonSubstitutionProcess::Create(nsite,nstate,0,nsite);
		DGamRateProcess::Create(nsite,ncat);
		PoissonGeneProfileProcess::Create(nsite,nstate,hub,offset);
	}

	virtual void Delete()	{
		PoissonGeneProfileProcess::Delete();
		DGamRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}
};

class PoissonGeneGammaPhyloProcess : public virtual PoissonPhyloProcess, public GenePhyloProcess, public virtual PoissonGeneSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	PoissonGeneGammaPhyloProcess(string indatafile, string treefile, int nratecat, int infixtopo, int indc, PoissonProfileProcess* inhub, int inoffset)	{

		myid = 1;
		nprocs = 2;
		fixtopo = infixtopo;
		dc = indc;
		hub = inhub;
		offset = inoffset;

		datafile = indatafile;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
		if (dc)	{
			plaindata->DeleteConstantSites();
		}
		const TaxonSet* taxonset = plaindata->GetTaxonSet();
		if (treefile == "None")	{
			tree = new Tree(taxonset);
			tree->MakeRandomTree();
		}
		else	{
			tree = new Tree(treefile);
		}
		tree->RegisterWith(taxonset,0);
		
		// cerr << "gene create\n";
		Create(tree,plaindata,nratecat,inhub,inoffset);
		// cerr << "gene create ok\n";

		/*
		cerr << "sample\n";
		Sample();
		cerr << "unfold\n";
		Unfold();
		cerr << "unfold ok\n";
		*/
		datafile = indatafile;
	}

	~PoissonGeneGammaPhyloProcess() {
		Delete();
	}

	StateSpace* GetStateSpace()	{
		return data->GetStateSpace();
	}

	virtual void GlobalUpdateSiteProfileSuffStat() {}
	virtual void SlaveUpdateSiteProfileSuffStat() {}

	void SetMixtureParameters()	{
		UpdateZip();
	}

	void TraceHeader(ostream& os)	{
	}

	void Trace(ostream& os)	{
	}

	void SaveTree(string name)	{
		ofstream os((name + datafile + ".treelist").c_str(),ios_base::app);
		GetLengthTree()->ToStream(os);
	}

	double Move(double tuning = 1.0)	{

		// Unfold();
		chronototal.Start();
		propchrono.Start();
		// cerr << "branch\n";
		NonMPIBranchLengthMove(tuning);
		NonMPIBranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			int nrep = (int) (GetGibbsFactor());
			// int nrep = (int) (GetGibbsFactor() * data->GetNtaxa());
			for (int rep=0; rep<nrep; rep++)	{
				NonMPIGibbsSPR();
			}
		}
		// cerr << "gibbs ok\n";
		propchrono.Stop();
		
		// cerr << "collapse\n";
		Collapse();
		// cerr << "collapse ok\n";

		// cerr << "branch process move\n";
		GammaBranchProcess::NonMPIMove(tuning,10);
		// cerr << "branch process move ok\n";

		// cerr << "rate process move\n";
		DGamRateProcess::NonMPIMove(0.3*tuning,10);
		DGamRateProcess::NonMPIMove(0.03*tuning,10);
		// cerr << "rate process move ok\n";

		/*
		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}
		*/
		chronototal.Stop();

		// Trace(cerr);
		// cerr << "gene move ok\n";

		return 1;
	}

	void ToStreamHeader(ostream& os)	{
		/*
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << GetNcat() << '\n';
		os << rrtype << '\n';
		os << fixtopo << '\n';
		os << dc << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
		*/
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
	}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int nratecat, PoissonProfileProcess* inhub, int inoffset)	{
		PoissonPhyloProcess::Create(intree,indata);
		PoissonGeneSubstitutionProcess::Create(indata->GetNsite(),nratecat,indata->GetNstate(),inhub,inoffset);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		PoissonGeneSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}

	int fixtopo;
	int dc;
	PoissonProfileProcess* hub;
	string datafile;
};

#endif

