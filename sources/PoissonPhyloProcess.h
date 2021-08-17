
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef POISSONPHYLO_H
#define POISSONPHYLO_H

#include "PhyloProcess.h"
#include "PoissonSubstitutionProcess.h"

class PoissonPhyloProcess : public virtual PhyloProcess, public virtual PoissonSubstitutionProcess	{

	public:

	PoissonPhyloProcess() : siteprofilesuffstatcount(0), allocsiteprofilesuffstatcount(0), zipdata(0) {}
	virtual ~PoissonPhyloProcess() {}

	void Unfold()	{
		UpdateZip();
		PhyloProcess::Unfold();
	}

	virtual void PrepareSiteLogLikelihood(int site) {
		UpdateZip(site);
	}

	void Collapse();

	// protected:

	// true data here !
	virtual void Create(Tree* intree, SequenceAlignment* indata);
	virtual void Delete();

	// in fact, same object as GetData, but now with its true type
	ZippedSequenceAlignment* GetZipData()	{
		if (! zipdata)	{
			cerr << "null zip\n";
			exit(1);
		}
		return zipdata;
	}

	virtual int GetNstate() {return truedata->GetNstate();}
	int GetZipSize(int site) {return GetZipData()->GetZipSize(site);}
	int GetOrbitSize(int site) {return GetZipData()->GetOrbitSize(site);}
	int GetStateFromZip(int site, int state) {return GetZipData()->GetStateFromZip(site,state);}
	bool InOrbit(int site, int state) {return GetZipData()->InOrbit(site,state);}
	
	void CreateSuffStat();
	void DeleteSuffStat();

	void UpdateSiteRateSuffStat();
	void UpdateSiteProfileSuffStat();
	void PoissonUpdateSiteProfileSuffStat();
	void UpdateBranchLengthSuffStat();

	void GlobalUpdateSiteProfileSuffStat();
	void SlaveUpdateSiteProfileSuffStat();

	int RecursiveUpdateSiteProfileSuffStat(const Link* from, int site);

	const int* GetSiteProfileSuffStatCount(int site) {return siteprofilesuffstatcount[site];}

	void SetDataFromLeaves()	{
		SampleTrueNodeStates(GetRoot());
		PhyloProcess::SetDataFromLeaves();
	}

	// virtual void RecursiveSimulateForward(const Link* from);

	void SampleTrueNodeStates(const Link* from);

	virtual double GetObservedCompositionalHeterogeneity(double* taxstat, double& meandist)	{
		return truedata->CompositionalHeterogeneity(taxstat,0,meandist);
	}

	void RecursiveUnzipBranchSitePath(const Link* from);
	void SlaveWriteMappings();

	void GlobalSetTestData();
	void SlaveSetTestData();

	// private:

	int** siteprofilesuffstatcount;
	int* allocsiteprofilesuffstatcount;

	ZippedSequenceAlignment* zipdata;
	SequenceAlignment* truedata;
};

#endif
