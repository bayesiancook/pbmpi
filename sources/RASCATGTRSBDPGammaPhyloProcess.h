
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATGTRSBDP_H
#define RASCATGTRSBDP_H

#include "ExpoConjugateGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "ExpoConjugateGTRSBDPProfileProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"

// this is the final class implementing the CATGTR phyloprocess

class RASCATGTRSBDPSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRSBDPProfileProcess {

	public:

	RASCATGTRSBDPSubstitutionProcess() {}
	virtual ~RASCATGTRSBDPSubstitutionProcess() {}

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int nratecat, int nstate, string inrrtype, int insitemin,int insitemax)	{
		ExpoConjugateGTRSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		DGamRateProcess::Create(nsite,nratecat);
		ExpoConjugateGTRSBDPProfileProcess::Create(nsite,nstate);
		SetRR(inrrtype);
	}

	virtual void Delete()	{
		ExpoConjugateGTRSBDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};


class RASCATGTRSBDPGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual RASCATGTRSBDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	using MixtureProfileProcess::LogStatPrior;

    virtual void SlaveExecute(MESSAGE);
	void GlobalUpdateParameters();
	void SlaveUpdateParameters();

    void GlobalSetSiteLogLCutoff();
    void SlaveSetSiteLogLCutoff();

	RASCATGTRSBDPGammaPhyloProcess(string indatafile, string treefile, int nratecat, int innmodemax, int iniscodon, GeneticCodeType incodetype, string inrrtype, int infixtopo, int inNSPR, int inNNNI, int inkappaprior, double indirweightprior, double inmintotweight, int indc, int incinit, int me, int np)	{
		myid = me;
		nprocs = np;

		InitIncremental = incinit;

		fixtopo = infixtopo;
		NSPR = inNSPR;
		NNNI = inNNNI;
		iscodon = iniscodon;
		codetype = incodetype;
		dc = indc;
		kappaprior = inkappaprior;
        dirweightprior = indirweightprior;
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

        SetNmodeMax(innmodemax);

		Create(tree,plaindata,nratecat,inrrtype,insitemin,insitemax);
		/*
		if (fixbl)	{
			SetLengthsFromNames();
		}
		*/
		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	RASCATGTRSBDPGammaPhyloProcess(istream& is, int me, int np)	{

		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		int nratecat;
		is >> nratecat;
		if (atof(version.substr(0,3).c_str()) > 1.8)	{
            int nmax;
            is >> nmax;
            SetNmodeMax(nmax);
        }
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
		string inrrtype;
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

		Create(tree,plaindata,nratecat,inrrtype,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}

	}

	~RASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

    void GlobalSetEmpiricalPrior(istream& is);
    void SlaveSetEmpiricalPrior();

    double GlobalGetSiteSteppingLogLikelihood(int site, int nrep, int restore);
    void SlaveGetSiteSteppingLogLikelihood();

	void TraceHeader(ostream& os)	{
		os << "iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		// os << "\tkappa\tallocent";
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

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
		// os << '\t' << kappa << '\t' << GetAllocEntropy();
		os << '\n';

	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		if (! fixbl)	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}
		if (! fixtopo)	{
			MoveTopo(NSPR,NNNI);
		}
		propchrono.Stop();

		
		for (int rep=0; rep<5; rep++)	{
			GlobalCollapse();

			if (! fixbl)	{
				GammaBranchProcess::Move(tuning,10);
				GammaBranchProcess::Move(0.1*tuning,10);
			}

			GlobalUpdateParameters();
			DGamRateProcess::Move(tuning,10);
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);

			GlobalUpdateParameters();
			ExpoConjugateGTRSBDPProfileProcess::Move(1,1,2);
			if (iscodon){
				ExpoConjugateGTRSBDPProfileProcess::Move(0.1,1,3);
				ExpoConjugateGTRSBDPProfileProcess::Move(0.01,1,3);
			}
			GlobalUpdateParameters();

			if ((! fixrr) && (! fixbl)){
				LengthRelRateMove(1,10);
				LengthRelRateMove(0.1,10);
				LengthRelRateMove(0.01,10);
			}

			GlobalUnfold();
		}

		chronototal.Stop();

		return 1;
	}

	virtual void PrepareSiteLogLikelihood(int site) {
        /*
		int cat = ExpoConjugateGTRSBDPProfileProcess::alloc[site];
		if (! matrixarray[cat])	{
			cerr << "error in prepare site log likelihood: matrix is not allocated\n";
			exit(1);
			// CreateMatrix(cat);
		}
		UpdateMatrix(cat);
        */
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
        // number of rate categories
		os << GetNcat() << '\n';
		if (atof(version.substr(0,3).c_str()) > 1.8)	{
            os << GetNmodeMax() << '\n';
        }
		if (atof(version.substr(0,3).c_str()) > 1.3)	{
            os << iscodon << '\n';
            os << codetype << '\n';
            os << kappaprior << '\n';
            os << mintotweight << '\n';
        }
        if (atof(version.substr(0,3).c_str()) > 1.7)	{
            os << dirweightprior << '\n';
        }
		os << rrtype << '\n';
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
		ExpoConjugateGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		ExpoConjugateGTRSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadRelRates(string name, int burnin, int every, int until, int verbose);
	void ReadSiteProfiles(string name, int burnin, int every, int until);
	void ReadPostHyper(string name, int burnin, int every, int until);
	void ReadSiteProfileSuffStat(string name, int burnin, int every, int until);

	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int nratecat, string inrrtype, int insitemin,int insitemax)	{
		ExpoConjugateGTRPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		RASCATGTRSBDPSubstitutionProcess::Create(indata->GetNsite(),nratecat,indata->GetNstate(),inrrtype,insitemin,insitemax);
		GammaBranchProcess::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATGTRSBDPSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

	int iscodon;
	GeneticCodeType codetype;
    double siteloglcutoff;
    /*
    double* empcount;
    double* empbeta;
    */
};

#endif

