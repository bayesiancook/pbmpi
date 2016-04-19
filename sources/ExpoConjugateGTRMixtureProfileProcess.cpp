
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugateGTRMixtureProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugateGTRMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void ExpoConjugateGTRMixtureProfileProcess::Create(int innsite, int indim)	{
	ExpoConjugateGTRProfileProcess::Create(innsite,indim);
	if (! profilesuffstatcount)	{
		GTRMixtureProfileProcess::Create(innsite,indim);
		profilesuffstatcount = new int*[GetNmodeMax()];
		profilesuffstatbeta = new double*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = new int[GetDim()];
			profilesuffstatbeta[i] = new double[GetDim()];
		}
	}
}

void ExpoConjugateGTRMixtureProfileProcess::Delete() {
	if (profilesuffstatcount)	{
		for (int i=0; i<GetNmodeMax(); i++)	{
			delete[] profilesuffstatcount[i];
			delete[] profilesuffstatbeta[i];
		}
		delete[] profilesuffstatcount;
		delete[] profilesuffstatbeta;
		profilesuffstatcount = 0;
		profilesuffstatbeta = 0;
		GTRMixtureProfileProcess::Delete();
	}
}


void ExpoConjugateGTRMixtureProfileProcess::UpdateModeProfileSuffStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
			profilesuffstatbeta[i][k] = 0;
		}
	}
	for (int i=0; i<GetNsite(); i++)	{
		const int* count = GetSiteProfileSuffStatCount(i);
		const double* beta = GetSiteProfileSuffStatBeta(i);
		int cat = alloc[i];
		/*
		if (cat < 0)	{
			cerr << "negative alloc\n";
			exit(1);
		}
		if (cat >= Ncomponent)	{
			cerr << "overflow\n";
			cerr << cat << '\t' << Ncomponent << '\n';
			exit(1);
		}
		*/
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[cat][k] += count[k];
			profilesuffstatbeta[cat][k] += beta[k];
		}
	}
}

double ExpoConjugateGTRMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += profilesuffstatcount[cat][k] * log(profile[cat][k]) - profilesuffstatbeta[cat][k] * profile[cat][k];
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

void ExpoConjugateGTRMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	for (int k=0; k<GetDim(); k++)	{

		int tmp = profilesuffstatcount[cat1][k];
		profilesuffstatcount[cat1][k] = profilesuffstatcount[cat2][k];
		profilesuffstatcount[cat2][k] = tmp;

		double temp = profilesuffstatbeta[cat1][k];
		profilesuffstatbeta[cat1][k] = profilesuffstatbeta[cat2][k];
		profilesuffstatbeta[cat2][k] = temp;
	}
}


//-------------------------------------------------------------------------
//	* udpate eq. frequency profiles based on sufficient statistics
//	(CPU level 3)
//-------------------------------------------------------------------------

double ExpoConjugateGTRMixtureProfileProcess::LogStatProb(int site, int cat)	{
	const int* count = GetSiteProfileSuffStatCount(site);
	const double* beta = GetSiteProfileSuffStatBeta(site);
	double* pi = profile[cat];
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += count[k] * log(pi[k]) - beta[k] * pi[k];
	}
	return total;
}



double ExpoConjugateGTRMixtureProfileProcess::PoissonDiffLogSampling(int cat, int site)	{

	const int* nsub = GetSiteProfileSuffStatCount(site);
	int* catnsub = profilesuffstatcount[cat];
	int totalsub = 0;
	double priorweight = 0;
	int grandtotal = 0;
	for (int k=0; k<GetDim(); k++)	{
		totalsub += nsub[k];
		priorweight += dirweight[k];
		grandtotal += catnsub[k];
	}
	
	double total = 0;
	for (int j=0; j< totalsub; j++)	{
		total -= log(priorweight + grandtotal + j);
	}
	for (int k=0; k<GetDim(); k++)	{
		for (int j=0; j< nsub[k]; j++)	{
			total += log(dirweight[k] + catnsub[k] + j);
		}
	}
	return total;
}


