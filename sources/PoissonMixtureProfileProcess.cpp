
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonMixtureProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PoissonMixtureProfileProcess::Create(int innsite, int indim)	{
	if (! profilesuffstatcount)	{
		PoissonProfileProcess::Create(innsite,indim);
		MixtureProfileProcess::Create(innsite,indim);
		profilesuffstatcount  = new int*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = new int[GetDim()];
		}
		// SampleProfile();
	}
}

void PoissonMixtureProfileProcess::Delete() {
	if (profilesuffstatcount)	{
		for (int i=0; i<GetNmodeMax(); i++)	{
			delete[] profilesuffstatcount[i];
		}
		delete[] profilesuffstatcount;
		profilesuffstatcount = 0;
		PoissonProfileProcess::Delete();
		MixtureProfileProcess::Delete();
	}
}

void PoissonMixtureProfileProcess::UpdateModeProfileSuffStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
		}
	}
	for (int i=0; i<GetNsite(); i++)	{
		const int* count = GetSiteProfileSuffStatCount(i);
		int cat = alloc[i];
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[cat][k] += count[k];
		}
	}
}

double PoissonMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	double priorweight = 0;
	double postweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[k] + profilesuffstatcount[cat][k]) - rnd::GetRandom().logGamma(dirweight[k]);
		priorweight += dirweight[k];
		postweight += dirweight[k] + profilesuffstatcount[cat][k];
	}
	total += rnd::GetRandom().logGamma(priorweight) - rnd::GetRandom().logGamma(postweight);
	return total;
}

double PoissonMixtureProfileProcess::DiffLogSampling(int cat, int site)	{

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

double PoissonMixtureProfileProcess::LogStatProb(int site, int cat)	{
	const int* nsub = GetSiteProfileSuffStatCount(site);
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += nsub[k] * log(profile[cat][k]);
	}
	return total;
}

double PoissonMixtureProfileProcess::MoveProfile()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		MoveProfile(i);
	}
	return 1;
}

double PoissonMixtureProfileProcess::MoveProfile(int cat)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[cat][k] = rnd::GetRandom().sGamma(dirweight[k] + profilesuffstatcount[cat][k]);
		if (profile[cat][k] < stateps)	{
			profile[cat][k] = stateps;
		}
		total += profile[cat][k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[cat][k] /= total;
	}
	return 1;
}

double PoissonMixtureProfileProcess::LogStatIntPrior()	{
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		if (occupancy[k])	{
			total += LogStatIntPrior(k);
		}
	}
	return total;
}

double PoissonMixtureProfileProcess::LogStatIntPrior(int cat)	{

	double total = 0;
	int tot = 0;
	double totweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[k] + profilesuffstatcount[cat][k]);
		total -= rnd::GetRandom().logGamma(dirweight[k]);
		totweight += dirweight[k];
		tot += profilesuffstatcount[cat][k];
	}
	total -= rnd::GetRandom().logGamma(totweight + tot);
	total += rnd::GetRandom().logGamma(totweight);
	return total;
}
	
double PoissonMixtureProfileProcess::MoveDirWeights(double tuning, int nrep)	{

	UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatIntPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[k] *= e;
			deltalogprob += LogHyperPrior() + LogStatIntPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				dirweight[k] /= e;
			}
		}
	}
	SampleStat();
	return naccepted / nrep / GetDim();
}
void PoissonMixtureProfileProcess::RemoveSite(int site, int cat)	{
	occupancy[cat] --;
	if (activesuffstat)	{
		const int* nsub = GetSiteProfileSuffStatCount(site);
		int* catnsub = profilesuffstatcount[cat];
		for (int k=0; k<GetDim(); k++)	{
			catnsub[k] -= nsub[k];
		}
	}
}

void PoissonMixtureProfileProcess::AddSite(int site, int cat)	{
	alloc[site] = cat;
	occupancy[cat] ++;
	UpdateZip(site);
	if (activesuffstat)	{
		const int* nsub = GetSiteProfileSuffStatCount(site);
		int* catnsub = profilesuffstatcount[cat];
		for (int k=0; k<GetDim(); k++)	{
			catnsub[k] += nsub[k];
		}
	}
}

void PoissonMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	for (int k=0; k<GetDim(); k++)	{
		int tmp = profilesuffstatcount[cat1][k];
		profilesuffstatcount[cat1][k] = profilesuffstatcount[cat2][k];
		profilesuffstatcount[cat2][k] = tmp;
	}
	/*
	int* tmp = profilesuffstatcount[cat1];
	profilesuffstatcount[cat1] = profilesuffstatcount[cat2];
	profilesuffstatcount[cat2] = tmp;
	*/
}

