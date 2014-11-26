
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSMATRIXPHYLO_H
#define GENPATHSSMATRIXPHYLO_H

#include "MatrixPhyloProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"

class GeneralPathSuffStatMatrixPhyloProcess : public virtual MatrixPhyloProcess, public virtual GeneralPathSuffStatMatrixSubstitutionProcess {

	public:

	GeneralPathSuffStatMatrixPhyloProcess() : siterootstate(0), sitepaircount(0), sitewaitingtime(0) {}
	virtual ~GeneralPathSuffStatMatrixPhyloProcess() {}

	// this is the log of the site likelihood?
	double logSiteProbPath(int site, SubMatrix* mat);

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata, int indim,int insitemin,int insitemax)	{

		GeneralPathSuffStatMatrixSubstitutionProcess::Create(indata->GetNsite(),indim,insitemin,insitemax);
		MatrixPhyloProcess::Create(intree,indata,indim);
	}

	virtual void Delete()	{
		MatrixPhyloProcess::Delete();
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
	}

	void Unfold();
	void Collapse();

	map<pair<int,int>,int>& GetSitePairCount(int site) {return sitepaircount[site];}
	int GetSiteRootState(int site) {return siterootstate[site];}
	map<int,double>& GetSiteWaitingTime(int site) {return sitewaitingtime[site];}

	// should also create the matrices
	void GlobalUnfold();
	// void GlobalCollapse();

	void CreateSuffStat();
	void DeleteSuffStat();

	void GlobalUpdateSiteProfileSuffStat();
	void SlaveUpdateSiteProfileSuffStat();

	void UpdateSiteRateSuffStat();
	void UpdateBranchLengthSuffStat();
	void UpdateSiteProfileSuffStat();

	int CountMapping(int site);
	// int CountMapping();
	// int GlobalCountMapping();

	int* siterootstate;
	map< pair<int,int>, int>* sitepaircount;
	map<int,double>* sitewaitingtime;
	
};

#endif

