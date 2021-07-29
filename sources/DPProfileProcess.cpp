
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "DPProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* DPProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

const double meandir[] = {0.499737,0.171262,0.183399,0.225593,0.197453,0.211819,0.173191,0.175454,0.3181,0.240008,0.187577,0.324778,0.205587,0.395097,0.162356,0.519427,0.526213,0.349177,0.0511527,0.130222};

void DPProfileProcess::SampleHyper()	{
	kappa = 1;
	// kappa = GetNsite() / 5 ;
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = 1;
		// dirweight[i] = meandir[i];
	}
}
	
void DPProfileProcess::PriorSampleHyper()	{
	kappa = 1;
	if (kappaprior == 0)	{
        double a = profilefrac + (1-profilefrac)*empkappaalpha;
        double b = profilefrac*0.1 + (1-profilefrac)*empkappabeta;
        kappa = rnd::GetRandom().Gamma(a,b);
        // kappa = 10.0 * rnd::GetRandom().sExpo();
	}
	else 	{
        double u = 8 * (rnd::GetRandom().Uniform() - 4.0);
        kappa = exp(log(10.0) * u);
	}
	for (int k=0; k<GetDim(); k++)	{
        double a = profilefrac + (1-profilefrac)*empdirweightalpha[k];
        double b = profilefrac + (1-profilefrac)*empdirweightbeta[k];
        dirweight[k] = rnd::GetRandom().Gamma(a,b);
		// dirweight[k] = rnd::GetRandom().sExpo();
	}
}
	
// 1 component
/*
void DPProfileProcess::SampleAlloc()	{
	CreateComponent(0);
	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{
		AddSite(i,0);
	}
}
*/
	
// Nsite compoents
/*
void DPProfileProcess::SampleAlloc()	{
	Ncomponent = GetNsite();
	for (int i=0; i<GetNsite(); i++)	{
		CreateComponent(i);
		AddSite(i,i);
	}
}
*/

double DPProfileProcess::LogProxy(int site, int cat)	{
	return 0;
}

void DPProfileProcess::SampleAlloc()	{

	CreateComponent(0);
	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{
		double* p = new double[Ncomponent+1];
		double total = 0;
		double max = 0;
		for (int k=0; k<Ncomponent; k++)	{
			double tmp = log(occupancy[k]) * LogProxy(i,k);
			if ((!k) || (max < tmp))	{
				max = tmp;
			}
			p[k] = tmp;
		}
		p[Ncomponent] = log(kappa) + LogProxy(i,Ncomponent);
		if (max < p[Ncomponent])	{
			max = p[Ncomponent];
		}
		for (int k=0; k<=Ncomponent; k++)	{
			double tmp = exp(p[k] - max);
			total += tmp;
			p[k] = total;
		}
		double q = total * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<=Ncomponent) && (q > p[k])) k++;
		if (k == Ncomponent+1)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		if (k==Ncomponent)	{
			CreateComponent(k);
			Ncomponent++;
		}
		AddSite(i,k);
		delete[] p;
	}
}

double DPProfileProcess::LogHyperPrior()	{
	double total = 0;
	if (kappaprior == 0)	{
        if (profilefrac == 1.0) {
            total = -kappa / 10.0;
        }
        else    {
            double a = profilefrac + (1-profilefrac)*empkappaalpha;
            double b = profilefrac*0.1 + (1-profilefrac)*empkappabeta;
            total += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(kappa) - b*kappa;
        }
	}
	else 	{
		total = -log(kappa);
		if ((kappa < 1e-4) || (kappa > 1e4))	{
			total = -std::numeric_limits<double>::infinity();
		}
	}
	double sum = 0;
	for (int k=0; k<GetDim(); k++)	{
        if (profilefrac == 1.0) {
            total -= dirweight[k];
        }
        else    {
            double a = profilefrac + (1-profilefrac)*empdirweightalpha[k];
            double b = profilefrac + (1-profilefrac)*empdirweightbeta[k];
            total += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(dirweight[k]) - b*dirweight[k];
        }
		sum += dirweight[k];
	}
	if (sum < GetMinTotWeight())	{
        total = -std::numeric_limits<double>::infinity();
	}
	return total;
}

double DPProfileProcess::LogAllocPrior()	{
	double total = GetNOccupiedComponent() * log(kappa);
	// double total = GetNcomponent() * log(kappa);
	for (int i=0; i<GetNsite(); i++)	{
		total -= log(kappa + i);
	}
	return total;
}

double DPProfileProcess::MoveHyper(double tuning, int nrep)	{
	double total = 0;
	if (!burnin)	{
		total += MoveKappa(tuning,nrep);
	}
    if (! dirweightprior)   {
        total += MoveDirWeights(tuning,nrep);
    }
	return total;
}

double DPProfileProcess::MoveKappa(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogAllocPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		kappa *= e;
		deltalogprob += LogHyperPrior() + LogAllocPrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		/*
		if (kappa < 1)	{
			accepted = 0;
		}
		*/
		if (accepted)	{
			naccepted++;
		}
		else	{
			kappa /= e;
		}
	}
	return naccepted / nrep;
}



