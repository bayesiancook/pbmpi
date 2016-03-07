
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionedExpoConjugateGTRPartitionedProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionedExpoConjugateGTRPartitionedProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PartitionedExpoConjugateGTRPartitionedProfileProcess::Create(int indim, PartitionScheme rrscheme, PartitionScheme statscheme)	{
	if (! profilesuffstatcount)	{
		PartitionedExpoConjugateGTRProfileProcess::Create(indim, rrscheme);
		PartitionedGTRPartitionedProfileProcess::Create(indim, rrscheme, statscheme);

		int Np = PartitionedProfileProcess::GetNpart();
		profilesuffstatcount = new int*[Np];
		profilesuffstatbeta = new double*[Np];
		for (int i=0; i< Np; i++)	{
			profilesuffstatcount[i] = new int[GetDim()];
			profilesuffstatbeta[i] = new double[GetDim()];
		}
	}
}

void PartitionedExpoConjugateGTRPartitionedProfileProcess::Delete() {
	if (profilesuffstatcount)	{
		for (int i=0; i<PartitionedProfileProcess::GetNpart(); i++)	{
			delete[] profilesuffstatcount[i];
			delete[] profilesuffstatbeta[i];
		}
		delete[] profilesuffstatcount;
		delete[] profilesuffstatbeta;
		profilesuffstatcount = 0;
		profilesuffstatbeta = 0;
		PartitionedGTRPartitionedProfileProcess::Delete();
		PartitionedExpoConjugateGTRProfileProcess::Delete();
	}
}


void PartitionedExpoConjugateGTRPartitionedProfileProcess::UpdateModeProfileSuffStat()	{
	for (int i=0; i<PartitionedProfileProcess::GetNpart(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
			profilesuffstatbeta[i][k] = 0;
		}
	}
	for (int i=0; i<GetNsite(); i++)	{
		const int* count = GetSiteProfileSuffStatCount(i);
		const double* beta = GetSiteProfileSuffStatBeta(i);
		int cat = PartitionedProfileProcess::GetSitePart(i);

		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[cat][k] += count[k];
			profilesuffstatbeta[cat][k] += beta[k];
		}
	}
}

double PartitionedExpoConjugateGTRPartitionedProfileProcess::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += profilesuffstatcount[cat][k] * log(profile[cat][k]) - profilesuffstatbeta[cat][k] * profile[cat][k];
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}


//-------------------------------------------------------------------------
//	* udpate eq. frequency profiles based on sufficient statistics
//	(CPU level 3)
//-------------------------------------------------------------------------

double PartitionedExpoConjugateGTRPartitionedProfileProcess::LogStatProb(int site, int cat)	{
	const int* count = GetSiteProfileSuffStatCount(site);
	const double* beta = GetSiteProfileSuffStatBeta(site);
	double* pi = profile[cat];
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += count[k] * log(pi[k]) - beta[k] * pi[k];
	}
	return total;
}

void PartitionedExpoConjugateGTRPartitionedProfileProcess::ToStream(ostream& os)	{

	if(nfreestat > 0)
	{
		if(nfreestat > 1)
			for (int j=0; j<GetDim(); j++)	{
				os << dirweight[j] << '\t';
			}
			os << '\n';

		for(int p = 0; p < PartitionedProfileProcess::GetNpart(); p++)
		{
			if(!fixstat[p])
			{
				for (int j=0; j<GetDim(); j++)	{
					os << profile[p][j] << '\t';
				}
				os << '\n';
			}
		}
	}

	for(int p = 0; p < PartitionedGTRProfileProcess::GetNpart(); p++)
	{
		if(!fixrr[p])
		{
			for (int i=0; i<GetNrr(); i++)	{
				os << rr[p][i] << '\t';
			}
			os << '\n';
		}
	}
}

void PartitionedExpoConjugateGTRPartitionedProfileProcess::FromStream(istream& is)	{

	if(nfreestat > 0)
	{
		if(nfreestat > 1)
			for (int i=0; i<GetDim(); i++)	{
				is >> dirweight[i];
			}

		for(int p = 0; p < PartitionedProfileProcess::GetNpart(); p++)
		{
			if(!fixstat[p])
				for (int j=0; j<GetDim(); j++)
					is >> profile[p][j];
		}
	}

	for(int p = 0; p < PartitionedGTRProfileProcess::GetNpart(); p++)
	{
		if(!fixrr[p])
			for (int i=0; i<GetNrr(); i++)
				is >> rr[p][i];
	}
}
