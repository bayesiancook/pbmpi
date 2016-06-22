
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSSUB_H
#define GENPATHSSSUB_H

#include "MatrixSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"

class GeneralPathSuffStatMatrixSubstitutionProcess : public virtual MatrixSubstitutionProcess, public virtual GeneralPathSuffStatMatrixProfileProcess {

	public:
	GeneralPathSuffStatMatrixSubstitutionProcess() {}
	virtual ~GeneralPathSuffStatMatrixSubstitutionProcess() {}

	protected:

	void AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray, int* nonmissing);
	void AddSiteRateSuffStat(int* count, double* beta, BranchSitePath** patharray, double length, int* nonmissing);
	void AddSiteProfileSuffStat(int* siterootstate, map<pair<int,int>, int>* sitepaircount, map<int,double>* sitewaitingtime, BranchSitePath** patharray, double branchlength, int* nonmissing);
	/*
	void AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray);
	void AddSiteRateSuffStat(int* count, double* beta, BranchSitePath** patharray, double length);
	// void AddSiteProfileSuffStat(int& siterootstate, map<pair<int,int>, int>& sitepaircount, map<int,double>& sitewaitingtime, BranchSitePath* path, double efflength, bool isroot);
	void AddSiteProfileSuffStat(int* siterootstate, map<pair<int,int>, int>* sitepaircount, map<int,double>* sitewaitingtime, BranchSitePath** patharray, double branchlength, bool isroot);
	*/

};

#endif

