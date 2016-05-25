
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATSBDP_H
#define RASCATSBDP_H

#include "RASCATGammaPhyloProcess.h"
#include "PoissonSBDPProfileProcess.h"
#include "CodonSequenceAlignment.h"

/*
#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "DGamRateProcess.h"
#include "PoissonDPProfileProcess.h"
#include "GammaBranchProcess.h"
*/

class RASCATSBDPSubstitutionProcess : public virtual RASCATSubstitutionProcess, public virtual PoissonSBDPProfileProcess {

	public:

	RASCATSBDPSubstitutionProcess() {}
	virtual ~RASCATSBDPSubstitutionProcess() {}

	protected:

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate,int insitemin,int insitemax)	{
		RASCATSubstitutionProcess::Create(nsite,ncat,nstate,insitemin,insitemax);
		// PoissonSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		// DGamRateProcess::Create(nsite,ncat);
		PoissonSBDPProfileProcess::Create(nsite,nstate);
	}

	virtual void Delete()	{
		PoissonSBDPProfileProcess::Delete();
		// DGamRateProcess::Delete();
		// PoissonSubstitutionProcess::Delete();
		RASCATSubstitutionProcess::Delete();
	}

};

class RASCATSBDPGammaPhyloProcess : public virtual RASCATGammaPhyloProcess, public virtual RASCATSBDPSubstitutionProcess {
// PoissonPhyloProcess, public virtual RASCATSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	RASCATSBDPGammaPhyloProcess(string indatafile, string treefile, int nratecat, int iniscodon, GeneticCodeType incodetype, int infixtopo, int inNSPR, int inNNNI, int inkappaprior, double inmintotweight, int indc, int ininc, int me, int np)	{
		myid = me;
		nprocs = np;

		InitIncremental = ininc;

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

		Create(tree,plaindata,nratecat,insitemin,insitemax);

		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	RASCATSBDPGammaPhyloProcess(istream& is, int me, int np)	{
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

		Create(tree,plaindata,nratecat,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		PoissonSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual double Move(double tuning = 1.0)	{
		chronototal.Start();
		propchrono.Start();
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo(NSPR,NNNI);
		}
		//raph
		/*
		if (! fixtopo)	{
			cout << "SPR " << GibbsSPR(4)<<'\n';	
			for(int i=0; i<1; i++){
				cout << "NNI"<<i<<' '<< GibbsNNI(0.1,1)<<'\n';
			}

		}
		*/
		propchrono.Stop();


		GlobalCollapse();

		GammaBranchProcess::Move(tuning,10);

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

		PoissonSBDPProfileProcess::Move(1,1,5);
		if (iscodon)	{
			PoissonSBDPProfileProcess::Move(0.1,1,15);
			PoissonSBDPProfileProcess::Move(0.01,1,15);
		}

		GlobalUnfold();
		chronototal.Stop();

		// Trace(cerr);

		return 1;
	
	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
	}

	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat,int insitemin,int insitemax)	{
		/*
		PoissonPhyloProcess::Create(intree,indata);
		RASCATSubstitutionProcess::Create(indata->GetNsite(),ncat,indata->GetNstate(),insitemin,insitemax);
		GammaBranchProcess::Create(intree);
		*/
		RASCATSBDPSubstitutionProcess::Create(indata->GetNsite(),ncat,indata->GetNstate(),insitemin,insitemax);
		RASCATGammaPhyloProcess::Create(intree,indata,ncat,insitemin,insitemax);
	}
		
	virtual void Delete()	{
		/*
		GammaBranchProcess::Delete();
		RASCATSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
		*/
		RASCATGammaPhyloProcess::Delete();
		RASCATSBDPSubstitutionProcess::Delete();
	}

	void SlaveExecute(MESSAGE signal);
};

#endif

