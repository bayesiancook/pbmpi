
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MixtureProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void MixtureProfileProcess::Create(int innsite, int indim)	{
	if (! profile)	{
		ProfileProcess::Create(innsite,indim);
		allocprofile = new double[GetNmodeMax() * GetDim()];
		profile = new double*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profile[i] = allocprofile + i*GetDim();
		}
		alloc = new int[GetNsite()];
		occupancy = new int[GetNmodeMax()];
		dirweight = new double[GetDim()];
		logstatprior = new double[GetNmodeMax()];
		profilesuffstatlogprob = new double[GetNmodeMax()];
        empdirweightalpha = new double[GetDim()];
        empdirweightbeta = new double[GetDim()];
        for (int i=0; i<GetDim(); i++)  {
            empdirweightalpha[i] = 1.0;
            empdirweightbeta[i] = 1.0;
            dirweight[i] = 1.0;
        }
	}
}

void MixtureProfileProcess::Delete()	{
	if (profile)	{
        delete[] empdirweightalpha;
        delete[] empdirweightbeta;
		delete[] profilesuffstatlogprob;
		delete[] logstatprior;
		delete[] allocprofile;
		delete[] profile;
		profile = 0;
		ProfileProcess::Delete();
	}
}

void MixtureProfileProcess::SetEmpiricalDirWeightPrior(double* inalpha, double* inbeta) {
    for (int k=0; k<GetDim(); k++)  {
        empdirweightalpha[k] = inalpha[k];
        empdirweightbeta[k] = inbeta[k];
    }
}

double MixtureProfileProcess::GetMeanDirWeight()	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += dirweight[k];
	}
	return total;
}

double MixtureProfileProcess::GetStatEnt()	{
	double total = 0;
	UpdateOccupancyNumbers();
	for (int k=0; k<GetNcomponent(); k++)	{
		total += occupancy[k] * GetStatEnt(k);
	}
	return total / GetNsite();
}

double MixtureProfileProcess::GetStatEnt(int k)	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		if (profile[k][i] <= 0)	{
			cerr << "error: 0 entry in profile\n";
			cerr << profile[k][i] << '\n';
			exit(1);
		}
		total -= profile[k][i] * log(profile[k][i]);
	}
	if (std::isnan(total))	{
		cerr << "entropy is nan\n";
		exit(1);
	}
	return  total;
}

double MixtureProfileProcess::GetCenterStatEnt()	{
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		totalweight += dirweight[k];
	}
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		double w = dirweight[k] / totalweight;
		total -= w * log(w);
	}
	return total;
}

void MixtureProfileProcess::RenormalizeProfiles()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += profile[i][k];
		}
		for (int k=0; k<GetDim(); k++)	{
			profile[i][k] /= total;
		}
	}
}

void MixtureProfileProcess::SampleProfile()	{
	SampleHyper();
	SampleAlloc();
	SampleStat();
	// UpdateComponents();
}

void MixtureProfileProcess::SampleStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		SampleStat(profile[i]);
	}
}

void MixtureProfileProcess::SampleStat(int i)	{
	SampleStat(profile[i]);
}

double MixtureProfileProcess::ResampleEmptyProfiles()	{

	UpdateOccupancyNumbers();
	for (int i=0; i<GetNcomponent(); i++)	{
		if (! occupancy[i])	{
			SampleStat(i);
		}
	}
	return 1.0;
}
/*
void MixtureProfileProcess::SampleStat(double* prof)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = rnd::GetRandom().sGamma(dirweight[k]);
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
	}
	total = 0;
	for (int k=0; k<GetDim(); k++)	{
		if (prof[k] < stateps)	{
			prof[k] = stateps;
		}
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
		if (isnan(prof[k]))	{
			cerr << "nan in sample stat\n";
			cerr << dirweight[k] << '\n';
			exit(1);
		}
	}
	// UpdateComponent(i);
}
*/
void MixtureProfileProcess::SampleStat(double* prof, double statmin)	{
	if (! statmin)	{
		statmin = stateps;
	}
	double total = 0;
	int infreached = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = rnd::GetRandom().sGamma(dirweight[k]);
		if (prof[k] < statmin)	{
			prof[k] = statmin;
			infreached = 1;
		}
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
	}
	if (infreached)	{
		statinfcount++;
	}
	totstatcount++;
}

void MixtureProfileProcess::SampleEmpiricalStat(double* prof, const double* count, double statmin)	{
	if (! statmin)	{
		statmin = stateps;
	}
	double total = 0;
	int infreached = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = rnd::GetRandom().sGamma(dirweight[k] + count[k]);
		if (prof[k] < statmin)	{
			prof[k] = statmin;
			infreached = 1;
		}
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
	}
	if (infreached)	{
		statinfcount++;
	}
	totstatcount++;
}

double MixtureProfileProcess::LogProfilePrior()	{
	double total = 0;
	total += LogHyperPrior();
	total += LogAllocPrior();
	total += LogStatPrior();
	return total;
}

void MixtureProfileProcess::UpdateOccupancyNumbers()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		occupancy[i] = 0;
	}
	for (int i=0; i<GetNsite(); i++)	{
		occupancy[alloc[i]]++;
	}
}

double MixtureProfileProcess::LogStatPrior()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		LogStatPrior(i);
	}
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += logstatprior[i];
	}
	return total;
}

double MixtureProfileProcess::LogStatPrior(int cat)	{
	double total = 0;
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += (dirweight[k] - 1) * log(profile[cat][k]) - rnd::GetRandom().logGamma(dirweight[k]);
		totalweight += dirweight[k];
	}
	total += rnd::GetRandom().logGamma(totalweight);
	logstatprior[cat] = total;
	return total;
}

double MixtureProfileProcess::LogStatPrior(const double* profile) {
	double total = 0;
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += (dirweight[k] - 1) * log(profile[k]) - rnd::GetRandom().logGamma(dirweight[k]);
		totalweight += dirweight[k];
	}
	total += rnd::GetRandom().logGamma(totalweight);
	return total;
}

double MixtureProfileProcess::EmpiricalLogStatPrior(const double* profile, const double* count) {
	double total = 0;
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += (dirweight[k] + count[k] - 1) * log(profile[k]) - rnd::GetRandom().logGamma(dirweight[k] + count[k]);
		totalweight += dirweight[k] + count[k];
	}
	total += rnd::GetRandom().logGamma(totalweight);
	return total;
}

double MixtureProfileProcess::ProfileSuffStatLogProb()	{
	// simply, sum over all components
	for (int i=0; i<GetNcomponent(); i++)	{
		ProfileSuffStatLogProb(i);
	}
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += profilesuffstatlogprob[i];
	}
	return total;
}

double MixtureProfileProcess::MoveDirWeights(double tuning, int nrep)	{
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[k] *= e;
			deltalogprob += LogHyperPrior() + LogStatPrior();
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
	ResampleEmptyProfiles();
	return naccepted / nrep / GetDim();
}

/*
double MixtureProfileProcess::MoveDirWeights(double tuning, int nrep)	{
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double bk = dirweight[k];
			double min = GetMinTotWeight() - GetMeanDirWeight() + dirweight[k];
			if (min < 0)	{
				min = 0;
			}
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			dirweight[k] += m;
			if (dirweight[k] < min)	{
				dirweight[k] = 2*min - dirweight[k];
			}
			deltalogprob += LogHyperPrior() + LogStatPrior();
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				dirweight[k] = bk;
			}
		}
	}
	return naccepted / nrep / GetDim();
}
*/

void MixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	int tmp = occupancy[cat1];
	occupancy[cat1] = occupancy[cat2];
	occupancy[cat2] = tmp;

	for (int k=0; k<GetDim(); k++)	{
		double tmp = profile[cat1][k];
		profile[cat1][k] = profile[cat2][k];
		profile[cat2][k] = tmp;
	}
	/*
	double* temp = profile[cat1];
	profile[cat1] = profile[cat2];
	profile[cat2] = temp;
	*/

	for (int i=0; i<GetNsite(); i++)	{
		if (alloc[i] == cat1)	{
			alloc[i] = cat2;
		}
		else if (alloc[i] == cat2)	{
			alloc[i] = cat1;
		}
	}
}

