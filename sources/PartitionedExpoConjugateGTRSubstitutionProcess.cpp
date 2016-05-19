
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionedExpoConjugateGTRSubstitutionProcess.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugate GTR Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//	* Gather sufficient statistics for resampling continuous parameters of the model
//	(branch lengths, relative rates, and site-specific rates and profiles)
//	just scanning through substitution mappings (singly linked lists of substitutions and waiting times)
//	counting them and summing over their exponential rates
//	(CPU level 1)
//-------------------------------------------------------------------------


void PartitionedExpoConjugateGTRSubstitutionProcess::AddRRSuffStat(int** rrsuffstatcount, double** rrsuffstatbeta, BranchSitePath** patharray, double branchlength)	{
	for (int i=sitemin; i<sitemax; i++)	{
	    if(!fixrr[GetSitePart(i)])
	        patharray[i]->AddRRSuffStat(rrsuffstatcount[GetSitePart(i)],rrsuffstatbeta[GetSitePart(i)],GetRate(i)*branchlength,GetProfile(i),GetNstate(i));
	}
}

void PartitionedExpoConjugateGTRSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength)	{
	for (int i=sitemin; i<sitemax; i++)	{
	    patharray[i]->AddRateSuffStat(siteratesuffstatcount[i],siteratesuffstatbeta[i],branchlength,GetRR(GetSitePart(i)),GetProfile(i),GetNstate(i));
	}
}

void PartitionedExpoConjugateGTRSubstitutionProcess::AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	    patharray[i]->AddRateSuffStat(count,beta,GetRate(i),GetRR(GetSitePart(i)),GetProfile(i),GetNstate(i));
	}
}

void PartitionedExpoConjugateGTRSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, double** siteprofilesuffstatbeta, BranchSitePath** patharray, double branchlength, bool isroot)	{
	if (!isroot)	{
		// non root case
		for (int i=sitemin; i<sitemax; i++)	{
		    patharray[i]->AddProfileSuffStat(siteprofilesuffstatcount[i],siteprofilesuffstatbeta[i],GetRate(i)*branchlength,GetRR(GetSitePart(i)),GetNstate(i));
		}
	}
	else	{
		// root case
		for (int i=sitemin; i<sitemax; i++)	{
		    siteprofilesuffstatcount[i][patharray[i]->GetInitState()]++;
		}
	}
}


