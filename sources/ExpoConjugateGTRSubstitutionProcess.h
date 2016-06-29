
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef EXPCONGTRSUB_H
#define EXPCONGTRSUB_H

#include "GTRSubstitutionProcess.h"
#include "ExpoConjugateGTRProfileProcess.h"
#include "RateProcess.h"

class ExpoConjugateGTRSubstitutionProcess : public virtual GTRSubstitutionProcess, public virtual ExpoConjugateGTRProfileProcess, public virtual RateProcess	{

	public:

	ExpoConjugateGTRSubstitutionProcess() {}
	virtual ~ExpoConjugateGTRSubstitutionProcess() {}

	virtual int GetNstate(int site) {return GetDim();}

	protected:

	void Create(int innsite, int indim,int insitemin,int insitemax)	{
		ExpoConjugateGTRProfileProcess::Create(innsite,indim);
		SubstitutionProcess::Create(innsite,indim,insitemin,insitemax);
	}

	void Delete() {
		SubstitutionProcess::Delete();
		ExpoConjugateGTRProfileProcess::Delete();
	}

	// CPU : level 1
	// gathering sufficient statistics from substitution mappings at each site, for a given branch
	void AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing);
	void AddRRSuffStat(int* rrsuffstatcount, double* rrsuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing);

	void AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray, int* nonmissing);
	void AddSiteProfileSuffStat(int** siteprofilesuffstatcount, double** siteprofilesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing);
	/*
	void AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength);
	void AddRRSuffStat(int* rrsuffstatcount, double* rrsuffstatbeta, BranchSitePath** patharray, double branchlength);

	void AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray);
	void AddSiteProfileSuffStat(int** siteprofilesuffstatcount, double** siteprofilesuffstatbeta, BranchSitePath** patharray, double branchlength, bool isroot);
	*/

};

#endif

