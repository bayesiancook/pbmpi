
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATGTR_H
#define RASCATGTR_H

#include "ExpoConjugateGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "ExpoConjugateGTRFiniteProfileProcess.h"
#include "GammaBranchProcess.h"

// this is the final class implementing the CATGTR phyloprocess

class RASCATGTRSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRFiniteProfileProcess {

	public:

	RASCATGTRSubstitutionProcess() {}
	virtual ~RASCATGTRSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int nratecat, int ncat, int nstate, int infixncomp, int inempmix, string inmixtype, string rrtype, int insitemin,int insitemax)	{
		if (ncat == -1)	{
			ncat = nsite;
		}
	
		ExpoConjugateGTRSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		DGamRateProcess::Create(nsite,nratecat);
		ExpoConjugateGTRFiniteProfileProcess::Create(nsite,nstate,ncat,infixncomp,inempmix,inmixtype,rrtype);
	}

	virtual void Delete()	{
		ExpoConjugateGTRFiniteProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};

class RASCATGTRFiniteGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual RASCATGTRSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);

	void GlobalUpdateParameters();

	RASCATGTRFiniteGammaPhyloProcess(string indatafile, string treefile, int nratecat, int innmodemax, int ncat, int infixncomp, int inempmix, string inmixtype, string inrrtype, double indirweightprior, int infixtopo, int inNSPR, int inNNNI, int indc, int me, int np)	{

		myid = me;
		nprocs = np;

        withfulllogl = 0;

        dirweightprior = indirweightprior;

		fixtopo = infixtopo;
		NSPR = inNSPR;
		NNNI = inNNNI;
		dc = indc;

		datafile = indatafile;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
		if (dc)	{
			plaindata->DeleteConstantSites();
		}
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

        SetNmodeMax(innmodemax);

		Create(tree,plaindata,nratecat,ncat,infixncomp,inempmix,inmixtype,inrrtype,insitemin,insitemax);
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	RASCATGTRFiniteGammaPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

        withfulllogl = 0;

		FromStreamHeader(is);
		is >> datafile;
		int nratecat;
		is >> nratecat;
		if (atof(version.substr(0,3).c_str()) > 1.8)	{
            int nmax;
            is >> nmax;
            SetNmodeMax(nmax);
        }
		int infixncomp;
		int intmp;
		int inempmix;
		string inmixtype;
		string inrrtype;
		is >> intmp >> inempmix >> inmixtype;
		is >> inrrtype;
        int ncat = 1;
        if (intmp)  {
            infixncomp = 1;
            ncat = intmp;
        }
        else    {
            infixncomp = 0;
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

		Create(tree,plaindata,nratecat,ncat,infixncomp,inempmix,inmixtype,inrrtype,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	~RASCATGTRFiniteGammaPhyloProcess() {
		Delete();
	}

    void GlobalSetEmpiricalPrior(istream& is);
    void SlaveSetEmpiricalPrior();

    double GlobalGetSiteSteppingLogLikelihood(int site, int nrep, int restore)  {
        if (fixncomp && (GetNcomponent() == 1))    {
            return PhyloProcess::GlobalGetSiteSteppingLogLikelihood(site, nrep, restore);
        }
        return GlobalGetSiteSteppingLogLikelihoodNonIS(site, nrep, restore);
    }

    void SlaveGetSiteSteppingLogLikelihood()    {
        if (fixncomp && (GetNcomponent() == 1))    {
            PhyloProcess::SlaveGetSiteSteppingLogLikelihood();
        }
        else    {
            SlaveGetSiteSteppingLogLikelihoodNonIS();
        }
    }

    double GlobalGetSiteSteppingLogLikelihoodNonIS(int site, int nrep, int restore);
    void SlaveGetSiteSteppingLogLikelihoodNonIS();


	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogLikelihood()	{
		return logL;
	}

	void SlaveUpdateParameters();

	void TraceHeader(ostream& os)	{
		os << "iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
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

        if (withfulllogl)   {
            os << '\t' << GlobalGetFullLogLikelihood();
        }
        else    {
            os << '\t' << GetLogLikelihood();
        }
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
		if (! fixbl)	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}
		if (!fixtopo)	{
			MoveTopo(NSPR,NNNI);
		}
		propchrono.Stop();

		GlobalCollapse();

		if (! fixbl)	{
			GammaBranchProcess::Move(tuning,50);
			GammaBranchProcess::Move(0.1*tuning,50);
		}

		GlobalUpdateParameters();
		DGamRateProcess::Move(tuning,50);
		DGamRateProcess::Move(0.3*tuning,50);
		DGamRateProcess::Move(0.03*tuning,50);

		ExpoConjugateGTRFiniteProfileProcess::Move(1,1,10);

		if ((! fixrr) && (! fixbl))	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		GlobalUnfold();

		chronototal.Stop();
		return 1;
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << GetNcat() << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.8)	{
            os << GetNmodeMax() << '\n';
        }
        int tmp = 0;
        if (fixncomp)   {
            tmp = GetNcomponent();
        }
        else    {
            tmp = 0;
        }
		os << tmp << '\t' << empmix << '\t' << mixtype << '\n';
		// os << fixncomp << '\t' << empmix << '\t' << mixtype << '\n';
		os << rrtype << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.7)	{
            os << dirweightprior << '\n';
        }
		os << fixtopo << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
            os << NSPR << '\t' << NNNI << '\n';
        }
		os << dc << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	virtual void PrepareSiteLogLikelihood(int site) {
		int cat = ExpoConjugateGTRFiniteProfileProcess::alloc[site];
		if (! matrixarray[cat])	{
			cerr << "error in prepare site log likelihood: matrix is not allocated\n";
			exit(1);
			// CreateMatrix(cat);
		}
		UpdateMatrix(cat);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		ExpoConjugateGTRFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		ExpoConjugateGTRFiniteProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();
	void ReadPostHyper(string name, int burnin, int every, int until);

	void ReadRelRates(string name, int burnin, int every, int until, int verbose);
	void ReadSiteProfiles(string name, int burnin, int every, int until);

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int nratecat,int ncat,int infixncomp, int inempmix, string inmixtype, string inrrtype, int insitemin,int insitemax)	{
		ExpoConjugateGTRPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		RASCATGTRSubstitutionProcess::Create(indata->GetNsite(),nratecat,ncat,indata->GetNstate(),infixncomp, inempmix, inmixtype, inrrtype, insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATGTRSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

    int withfulllogl;
};

#endif

