
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENE_H
#define MULTIGENE_H

#include "GeneProfileProcess.h"
#include "SBDPProfileProcess.h"

class MultiGeneMixture : public virtual SBDPProfileProcess	{

	public:

	MultiGeneMixture() : gibbsfactor(0.1) {}
	virtual ~MultiGeneMixture() {}

	double GetLogLikelihood();
	double GetMeanLength();
	double GetMeanAlpha();

	virtual void Trace(ostream& os) = 0;
	virtual void TraceHeader(ostream& os) = 0;

	void MakeFiles();

	StateSpace* GetStateSpace()	{
		return statespace;
	}
		
	double GetGibbsFactor()	{
		return gibbsfactor;
	}

	string GetVersion()	{
		return "1.3";
	}

	int GetNprocs()	{
		return nprocs;
	}

	int GetMyid()	{
		return myid;
	}

	int GetSiteMin()	{
		return 0;
	}

	int GetSiteMax()	{
		return GetGlobalNsite();
	}

	int GetGlobalNsite(int proc = -1)	{
		if (proc == -1)	{
			proc = myid;
		}
		return globalnsite[proc];
	}

	void ToStreamHeader(ostream& os)	{
		os << GetVersion() << '\n';
		propchrono.ToStream(os);
		chronototal.ToStream(os);
		os << name << '\n';
		os << datafile << '\n';
		os << nratecat << '\n';
		os << fixtopo << '\n';
		os << dc << '\n';
		os << kappaprior << '\n';
		os << gibbsfactor << '\n';
	}

	void FromStreamHeader(istream& is)	{
		string version;
		is >> version;
		if ((version.substr(0,3) != "1.3") && (version.substr(0,3) != "1.2"))	{
			cerr << "error: version too old: " << version << '\n';
			exit(1);
		}
		propchrono.FromStream(is);
		chronototal.FromStream(is);

		is >> name;
		is >> datafile;
		is >> nratecat;
		is >> fixtopo;
		is >> dc;
		is >> kappaprior;
		is >> gibbsfactor;
	}

	virtual void ToStream(ostream& os) = 0;

	virtual void FromStream(istream& is) = 0;

	void GlobalToStream(ostream& os);
	void GlobalFromStream(istream& is);
	void SlaveToStream();
	void SlaveFromStream();

	virtual double Move(double tuning = 1.0) = 0;

	void WaitLoop();
	void SaveTrees();
	void SlaveSaveTrees();

	protected:

        virtual void SlaveExecute(MESSAGE);

	// virtual double GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep) = 0;
	virtual void SlaveMixMove() = 0;

	void GlobalSample();
	void GlobalUnfold();
	void GlobalCollapse();
	void GlobalGeneMove();

	void SlaveSample();
	void SlaveUnfold();
	void SlaveCollapse();
	void SlaveGeneMove();

	void AllocateAlignments(string datafile, string treefile, int dc);
	virtual void Create(int inNsite, int Nstate, int nratecat, int dc) = 0;
	virtual void CreateSuffStat() = 0;
	virtual void DeleteSuffStat() = 0;

	virtual void GlobalUpdateParameters() = 0;
	virtual void SlaveUpdateParameters() = 0;

	void GlobalCollectGeneLikelihoods();
	void SlaveSendGeneLikelihoods();

	void GlobalCollectGeneLengths();
	void SlaveSendGeneLengths();

	void GlobalCollectGeneAlphas();
	void SlaveSendGeneAlphas();

	virtual void GlobalUpdateSiteProfileSuffStat() = 0;
	virtual void SlaveUpdateSiteProfileSuffStat() = 0;
	void UpdateSiteProfileSuffStat() {
		cerr << "error: in update site profile suff stat of global mixture\n";
		exit(1);
	}

	/*
	void GlobalUpdateRRSuffStat();
	void SlaveUpdateRRSuffStat();
	void UpdateRRSuffStat() {
		cerr << "error: in update site rr suff stat of global mixture\n";
		exit(1);
	}
	*/

	/*
	int* GetSiteProfileSuffStatCount(int site)	{
		return siteprofilesuffstatcount[site];
	}

	double* GetSiteProfileSuffStatBeta(int site)	{
		return siteprofilesuffstatbeta[site];
	}
	*/

	int myid;
	int nprocs;

	int GlobalNsite;
	int* globalnsite;
	int Nstate;
	StateSpace* statespace;
	int Ngene;
	int* genealloc;
	int* genesize;
	string* genename;
	string* treename;

	double lnL;
	double* genelnL;
	double* genealpha;
	double* genelength;
	double* tmpgenelnL;

	int fixtopo;
	int dc;

	string name;

	Chrono propchrono;
	Chrono chronototal;
	
	GenePhyloProcess** process;

	string datafile;
	int nratecat;
	double gibbsfactor;
};

#endif
