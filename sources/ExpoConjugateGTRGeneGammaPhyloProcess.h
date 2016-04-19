
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef EXPOCONJGENE_H
#define EXPOCONJGENE_H

#include "GeneProfileProcess.h"

#include "ExpoConjugateGTRPhyloProcess.h"
#include "GammaBranchProcess.h"
#include "DGamRateProcess.h"

class ExpoConjugateGTRGeneProfileProcess : public virtual ExpoConjugateGTRProfileProcess, public virtual GeneProfileProcess	{

	public:

	ExpoConjugateGTRGeneProfileProcess() : expohub(0) {}
	virtual ~ExpoConjugateGTRGeneProfileProcess() {}

	void SetMixtureParameters()	{
		const double* tmprr = expohub->GetRR();
		for (int k=0; k<GetNrr(); k++)	{
			rr[k] = tmprr[k];
		}
		/*
		for (int i=0; i<GetNsite(); i++)	{
			alloc[i] = expohub->alloc[i+offset];
		}
		*/
	}

	protected:

	virtual void Create(int innsite, int innstate, ExpoConjugateGTRProfileProcess* inhub, int inoffset)	{
		ExpoConjugateGTRProfileProcess::Create(innsite,innstate);
		GeneProfileProcess::Create(innsite,innstate,inhub,inoffset);
		expohub = inhub;
	}

	virtual void Delete()	{
		GeneProfileProcess::Delete();
		ExpoConjugateGTRProfileProcess::Delete();
	}

	virtual SubMatrix* GetMatrix(int site)	{
		return expohub->GetMatrix(site+offset);
	}

	// get rr

	// should deactivate the move rr
	// should sync the rr between this and the hub

	virtual void UpdateMatrices() {
		cerr << "error: in gene update matrices\n";
		exit(1);
	}

	virtual void CreateMatrices()	{
		// cerr << "errror : in gene create matrices\n";
		// exit(1);
	}

	virtual void DeleteMatrices()	{
		// cerr << "errror : in gene delete matrices\n";
		// exit(1);
	}

	/*
	virtual const int* GetSiteProfileSuffStatCount(int site)	{
		return expohub->GetSiteProfileSuffStatCount(site+offset);
	}

	virtual const double* GetSiteProfileSuffStatBeta(int site)	{
		return expohub->GetSiteProfileSuffStatBeta(site+offset);
	}
	*/

	ExpoConjugateGTRProfileProcess* expohub;

};


class ExpoConjugateGTRGeneSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRGeneProfileProcess {

	public:

	ExpoConjugateGTRGeneSubstitutionProcess() {}
	virtual ~ExpoConjugateGTRGeneSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate, ExpoConjugateGTRProfileProcess* hub, int offset)	{
		ExpoConjugateGTRSubstitutionProcess::Create(nsite,nstate,0,nsite);
		DGamRateProcess::Create(nsite,ncat);
		ExpoConjugateGTRGeneProfileProcess::Create(nsite,nstate,hub,offset);
		// SetRR(inrrtype);
	}

	virtual void Delete()	{
		ExpoConjugateGTRGeneProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};


class ExpoConjugateGTRGeneGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public GenePhyloProcess, public virtual ExpoConjugateGTRGeneSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	ExpoConjugateGTRGeneGammaPhyloProcess(string indatafile, string treefile, int nratecat, int infixtopo, int indc, ExpoConjugateGTRProfileProcess* inhub, int inoffset)	{

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

	~ExpoConjugateGTRGeneGammaPhyloProcess() {
		Delete();
	}

	StateSpace* GetStateSpace()	{
		return data->GetStateSpace();
	}

	virtual void GlobalUpdateSiteProfileSuffStat() {}
	virtual void SlaveUpdateSiteProfileSuffStat() {}

	virtual void GlobalUpdateRRSuffStat() {}
	virtual void SlaveUpdateRRSuffStat() {}

	void TraceHeader(ostream& os)	{
	}

	void Trace(ostream& os)	{
	}

	/*
	void MakeTreeFiles(string name)	{
		ofstream ((name + datafile + ".treelist").c_str());
	}
	*/

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

	virtual void Create(Tree* intree, SequenceAlignment* indata, int nratecat, ExpoConjugateGTRProfileProcess* inhub, int inoffset)	{
		ExpoConjugateGTRPhyloProcess::Create(intree,indata,indata->GetNstate(),0,indata->GetNsite());
		ExpoConjugateGTRGeneSubstitutionProcess::Create(indata->GetNsite(),nratecat,indata->GetNstate(),inhub,inoffset);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		ExpoConjugateGTRGeneSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

	int fixtopo;
	int dc;
	ExpoConjugateGTRProfileProcess* hub;
	string datafile;
};

#endif

