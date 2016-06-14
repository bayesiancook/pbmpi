
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugateGTRSubstitutionProcess.h"


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


void ExpoConjugateGTRSubstitutionProcess::AddRRSuffStat(int* rrsuffstatcount, double* rrsuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,GetRate(i)*branchlength,GetProfile(i),GetNstate(i));
		}
	}
}

void ExpoConjugateGTRSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddRateSuffStat(siteratesuffstatcount[i],siteratesuffstatbeta[i],branchlength,GetRR(),GetProfile(i),GetNstate(i));
		}
	}
}

void ExpoConjugateGTRSubstitutionProcess::AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddRateSuffStat(count,beta,GetRate(i),GetRR(),GetProfile(i),GetNstate(i));
		}
	}
}

void ExpoConjugateGTRSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, double** siteprofilesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			patharray[i]->AddProfileSuffStat(siteprofilesuffstatcount[i],siteprofilesuffstatbeta[i],GetRate(i)*branchlength,GetRR(),GetNstate(i));
		}
		else if (nonmissing[i] == 2)	{
			siteprofilesuffstatcount[i][patharray[i]->GetFinalState()]++;
		}
	}
}

/*
void ExpoConjugateGTRSubstitutionProcess::AddRRSuffStat(int* rrsuffstatcount, double* rrsuffstatbeta, BranchSitePath** patharray, double branchlength)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		patharray[i]->AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,GetRate(i)*branchlength,GetProfile(i),GetNstate(i));
	}
}

void ExpoConjugateGTRSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		patharray[i]->AddRateSuffStat(siteratesuffstatcount[i],siteratesuffstatbeta[i],branchlength,GetRR(),GetProfile(i),GetNstate(i));
	}
}

void ExpoConjugateGTRSubstitutionProcess::AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		patharray[i]->AddRateSuffStat(count,beta,GetRate(i),GetRR(),GetProfile(i),GetNstate(i));
	}
}

void ExpoConjugateGTRSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, double** siteprofilesuffstatbeta, BranchSitePath** patharray, double branchlength, bool isroot)	{
	if (!isroot)	{
		// non root case
		for (int i=sitemin; i<sitemax; i++)	{
		// for (int i=0; i<GetNsite(); i++)	{
			patharray[i]->AddProfileSuffStat(siteprofilesuffstatcount[i],siteprofilesuffstatbeta[i],GetRate(i)*branchlength,GetRR(),GetNstate(i));
		}
	}
	else	{
		// root case
		for (int i=sitemin; i<sitemax; i++)	{
		// for (int i=0; i<GetNsite(); i++)	{
			siteprofilesuffstatcount[i][patharray[i]->GetInitState()]++;
		}
	}
}

*/
