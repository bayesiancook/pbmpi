
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* General Path SuffStat Matrix Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixSubstitutionProcess::AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddGeneralPathRateSuffStat(count,beta,GetRate(i),GetMatrix(i));
		}
	}
}

void GeneralPathSuffStatMatrixSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddGeneralPathRateSuffStat(siteratesuffstatcount[i],siteratesuffstatbeta[i],branchlength,GetMatrix(i));
		}
	}
}

void GeneralPathSuffStatMatrixSubstitutionProcess::AddSiteProfileSuffStat(int* siterootstate, map<pair<int,int>, int>* sitepaircount, map<int,double>* sitewaitingtime, BranchSitePath** patharray, double branchlength, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddGeneralPathSuffStat(sitepaircount[i],sitewaitingtime[i],GetRate(i)*branchlength);
		}
		else if (nonmissing[i] == 2)	{
			siterootstate[i] = patharray[i]->GetFinalState();
		}
	}
}

/*
void GeneralPathSuffStatMatrixSubstitutionProcess::AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		patharray[i]->AddGeneralPathRateSuffStat(count,beta,GetRate(i),GetMatrix(i));
	}
}

void GeneralPathSuffStatMatrixSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		patharray[i]->AddGeneralPathRateSuffStat(siteratesuffstatcount[i],siteratesuffstatbeta[i],branchlength,GetMatrix(i));
	}
}

void GeneralPathSuffStatMatrixSubstitutionProcess::AddSiteProfileSuffStat(int* siterootstate, map<pair<int,int>, int>* sitepaircount, map<int,double>* sitewaitingtime, BranchSitePath** patharray, double branchlength, bool isroot)	{
	if (!isroot)	{
		// non root case
		for (int i=sitemin; i<sitemax; i++)	{
		// for (int i=0; i<GetNsite(); i++)	{
			// patharray[i]->AddGeneralPathSuffStat(sitepaircount[i], sitewaitingtime[i], efflength);
			patharray[i]->AddGeneralPathSuffStat(sitepaircount[i],sitewaitingtime[i],GetRate(i)*branchlength);
		}
	}
	else	{
		// root case
		for (int i=sitemin; i<sitemax; i++)	{
		// for (int i=0; i<GetNsite(); i++)	{
			siterootstate[i] = patharray[i]->GetInitState();
		}
	}
}

*/
