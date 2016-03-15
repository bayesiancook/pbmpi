
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTEXPCONGTRPHYLO_H
#define PARTEXPCONGTRPHYLO_H

#include "PartitionedGTRPhyloProcess.h"
#include "PartitionedExpoConjugateGTRSubstitutionProcess.h"
#include "RateProcess.h"

// this is the Exponential Conjugate version
class PartitionedExpoConjugateGTRGammaPhyloProcess : public virtual PartitionedGTRGammaPhyloProcess, public virtual PartitionedExpoConjugateGTRSubstitutionProcess {

	public:

	PartitionedExpoConjugateGTRGammaPhyloProcess() : siteprofilesuffstatcount(0), siteprofilesuffstatbeta(0), allocsiteprofilesuffstatcount(0), allocsiteprofilesuffstatbeta(0) {}
	virtual ~PartitionedExpoConjugateGTRGammaPhyloProcess() {}

	protected:

	virtual void Create(Tree* intree, SequenceAlignment* indata,int indim, PartitionScheme rrscheme, int ncat, bool inlinkmult, PartitionScheme dgamscheme, int insitemin,int insitemax)	{
		PartitionedExpoConjugateGTRSubstitutionProcess::Create(indata->GetNstate(), rrscheme, insitemin,insitemax);
		PartitionedGTRGammaPhyloProcess::Create(intree,indata, indata->GetNstate(), rrscheme, ncat, inlinkmult, dgamscheme, insitemin,insitemax);
	}

	virtual void Delete()	{
		PartitionedGTRGammaPhyloProcess::Delete();
		PartitionedExpoConjugateGTRSubstitutionProcess::Delete();
	}

	/*
	void Unfold();
	void Collapse();
	*/

	const int* GetSiteProfileSuffStatCount(int site) {return siteprofilesuffstatcount[site];}
	const double* GetSiteProfileSuffStatBeta(int site) {return siteprofilesuffstatbeta[site];}

	// protected:
	public:

	void CreateSuffStat();
	void DeleteSuffStat();

	void GlobalUpdateRRSuffStat();
	void GlobalUpdateSiteProfileSuffStat();

	void SlaveUpdateRRSuffStat();
	void SlaveUpdateSiteProfileSuffStat();

	void UpdateSiteRateSuffStat();
	void UpdateBranchLengthSuffStat();
	void UpdateRRSuffStat();
	void UpdateSiteProfileSuffStat();

	int** siteprofilesuffstatcount;
	double** siteprofilesuffstatbeta;

	int* allocsiteprofilesuffstatcount;
	double* allocsiteprofilesuffstatbeta;

	int* tmpcount;
	double* tmpbeta;
};

#endif

