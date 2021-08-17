
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef POISSONSUB_H
#define POISSONSUB_H

#include "SubstitutionProcess.h"
#include "PoissonProfileProcess.h"
#include "RateProcess.h"

class PoissonSubstitutionProcess : public virtual SubstitutionProcess, public virtual PoissonProfileProcess, public virtual RateProcess {

	public:

	PoissonSubstitutionProcess() : zipstat(0) {}
	virtual ~PoissonSubstitutionProcess() {}


	const double* GetStationary(int site) {return zipstat[site];}
	int GetNstate(int site) {return GetZipSize(site);}

	virtual void ConditionalLikelihoodsToStatePostProbs(double*** aux,double*** statepostprob, int nodelabel, bool condalloc = false);

	protected:

	// CPU Level 3: implementations of likelihood propagation and substitution mapping methods
	void Propagate(double*** from, double*** to, double time, bool condalloc = false);

	// CPU Level 3: implementations of likelihood propagation and substitution mapping methods
	void SitePropagate(int site, double** from, double** to, double time, bool condalloc = false);

	BranchSitePath** SamplePaths(int* stateup, int* statedown, double time);
	BranchSitePath** SampleRootPaths(int* rootstate);

	void SimuPropagate(int* stateup, int* statedown, double time);
	void SimuPropagateZip(int* stateup, int* statedown, double time);

	// CPU Level 1: gathering sufficient statistics from substitution mappings
	void AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, double branchlength, BranchSitePath** patharray, int* nonmissing);
	void AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray, int* nonmissing);
	void AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root);

	/*
	void AddSiteRateSuffStat(int* siteratesuffstatcount, BranchSitePath** patharray);
	void AddBranchLengthSuffStat(int& count, BranchSitePath** patharray);
	void AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root);
	*/

	virtual void Delete() {
		DeleteZip();
		SubstitutionProcess::Delete();
	}

	void CreateZip();
	void DeleteZip();
	void UpdateZip();
	void UpdateZip(int site);
	void UnzipBranchSitePath(BranchSitePath** patharray, int* nodestateup, int* nodestatedown);
	int GetRandomStateFromZip(int site, int state);

	// will be implemented in PoissonPhyloProcess
	virtual int GetZipSize(int site) = 0;
	virtual int GetOrbitSize(int site) = 0;
	virtual int GetStateFromZip(int site, int state) = 0;
	virtual bool InOrbit(int site, int state) = 0;

	void ChooseTrueStates(BranchSitePath** patharray, int* nodestateup, int* nodestatedown, bool root);
	void ChooseRootTrueStates(int* nodestate);
	void ZipToTruePostProbs(double*** statepostprob, int nodelabel);

	private:

	double** zipstat;
};

#endif

