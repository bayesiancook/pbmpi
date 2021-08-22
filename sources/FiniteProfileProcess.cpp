
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "BiologicalSequences.h"
#include "FiniteProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* FiniteProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void FiniteProfileProcess::Create(int innsite, int indim, int ncat, int infixncomp,int inempmix, string inmixtype)	{
	if (! weight)	{
		empmix = inempmix;
		mixtype = inmixtype;
        fixncomp = infixncomp;
        if (empmix) {
            ReadNcomponent(mixtype);
            Ncomponent = Ncat;
        }
        else    {
            Ncomponent = ncat;
        }
        if (Ncomponent > GetNmodeMax()) {
            cerr << "error in FiniteProfileProcess::Create: max number of components is " << GetNmodeMax() << '\n';
            exit(1);
        }
		MixtureProfileProcess::Create(innsite,indim);
		weight = new double[GetNmodeMax()];
        empdirweight = new double[GetDim()];
        for (int k=0; k<GetDim(); k++)  {
            empdirweight[k] = 1.0;
        }
	}
}

void FiniteProfileProcess::Delete()	{
	if (weight)	{
        delete[] empdirweight;
		if (statfix)	{
			for (int i=0; i<Ncat; i++)	{
				delete[] statfix[i];
			}
			delete[] statfix;
			delete[] empweight;
		}
		delete[] weight;
		MixtureProfileProcess::Delete();
	}
}

void FiniteProfileProcess::SampleWeights()	{
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		weight[k] = rnd::GetRandom().sGamma(weightalpha);
		total += weight[k];
	}
	for (int k=0; k<GetNcomponent(); k++)	{
		weight[k] /= total;
	}
}

void FiniteProfileProcess::ResampleWeights()	{
	UpdateOccupancyNumbers();
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		weight[k] = rnd::GetRandom().sGamma(weightalpha + occupancy[k]);
		if (weight[k] < 1e-10)	{
			weight[k] = 1e-10;
		}
		total += weight[k];
	}
	for (int k=0; k<GetNcomponent(); k++)	{
		weight[k] /= total;
	}
}

void FiniteProfileProcess::SampleHyper()	{
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = 1.0;
	}
}
	
void FiniteProfileProcess::PriorSampleHyper()   {
    if (! dirweightprior)   {
        for (int k=0; k<GetDim(); k++)	{
            double a = profilefrac + (1-profilefrac) * empdirweightalpha[k];
            double b = profilefrac + (1-profilefrac) * empdirweightbeta[k];
            dirweight[k] = rnd::GetRandom().Gamma(a,b);
        }
    }
    else    {
        if (fixncomp && (GetNcomponent() == 1))   {
            for (int k=0; k<GetDim(); k++)  {
                dirweight[k] = profilefrac + (1.0-profilefrac) * empdirweight[k];
            }
        }
        else    {
            for (int k=0; k<GetDim(); k++)  {
                dirweight[k] = 1.0;
            }
        }
    }
}

void FiniteProfileProcess::SampleAlloc()	{
	if (!GetNcomponent())	{
		cerr << "error in sample alloc: " << GetNcomponent() << '\n';
		exit(1);
	}
	if (empmix)	{
		ReadStatFix(mixtype);
		SetStatFix();
	}
	SampleWeights();
	if (GetNcomponent() == GetNsite())	{
		for (int i=0; i<GetNsite(); i++)	{
			AddSite(i,i);
		}	
	}
	else	{
		for (int i=0; i<GetNsite(); i++)	{
			int choose = rnd::GetRandom().FiniteDiscrete(GetNcomponent(),weight);
			AddSite(i,choose);
		}
	}
	ResampleWeights();
}

void FiniteProfileProcess::DrawProfileFromPrior()	{

	if (! GetMyid())	{
		cerr << "error: in master DrawProfileFromPrior\n";
		exit(1);
	}

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		RemoveSite(i,alloc[i]);
		int choose = rnd::GetRandom().FiniteDiscrete(GetNcomponent(),weight);
		AddSite(i,choose);
	}
}


void FiniteProfileProcess::SampleStat()	{
	if (! empmix)	{
		MixtureProfileProcess::SampleStat();
	}
}

double FiniteProfileProcess::LogWeightPrior()	{

	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		total += rnd::GetRandom().logGamma(weightalpha + occupancy[k]);
	}
	total -= rnd::GetRandom().logGamma(GetNcomponent() * weightalpha + GetNsite());
	total += rnd::GetRandom().logGamma(GetNcomponent()*weightalpha) - GetNcomponent() * rnd::GetRandom().logGamma(weightalpha);
	return total;
}

double FiniteProfileProcess::LogHyperPrior()	{
	double total = 0;
    if (! dirweightprior)   {
        double sum = 0;
        for (int k=0; k<GetDim(); k++)	{
            if (profilefrac == 1.0) {
                total -= dirweight[k];
            }
            else    {
                double a = profilefrac + (1-profilefrac) * empdirweightalpha[k];
                double b = profilefrac + (1-profilefrac) * empdirweightbeta[k];
                total += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(dirweight[k]) - b*dirweight[k];
            }
            sum += dirweight[k];
        }
        if (sum < GetMinTotWeight())	{
            total = -std::numeric_limits<double>::infinity();
        }
    }
	total -= weightalpha;
	return total;
}

double FiniteProfileProcess::MoveHyper(double tuning, int nrep)	{
	double total = 0;
    if (! dirweightprior)   {
		total += MoveDirWeights(tuning,nrep);
	}
	return total;
}

double FiniteProfileProcess::MoveWeightAlpha(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogWeightPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		weightalpha *= e;
		deltalogprob += LogHyperPrior() + LogWeightPrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			weightalpha /= e;
		}
	}
	return naccepted / nrep;
}

double FiniteProfileProcess::MoveNcomponent(int nrep)	{

    UpdateOccupancyNumbers();
	int Kmax = GetNmodeMax();
    int K = GetNcomponent();
    int N = GetNsite();

    int k0 = 0;
    for (int k=0; k<K; k++) {
        if (! occupancy[k])   {
            k0++;
        }
    }

    for (int rep=0; rep<nrep; rep++)    {
        if (k0 == 0)    {
            double ratio = 0.5 * (K+1) / (k0+1) * K / (N+K);
            if (K == Kmax)  {
                ratio = 0;
            }
            if (rnd::GetRandom().Uniform() < ratio)  {
                occupancy[K] = 0;
                K++;
                k0++;
            }
        }
        else    {
            if (rnd::GetRandom().Uniform() < 0.5)    {
                double ratio = (K+1.0) / (k0+1.0) * K / (N+K);
                if (K == Kmax)  {
                    ratio = 0;
                }
                if (rnd::GetRandom().Uniform() < ratio)  {
                    occupancy[K] = 0;
                    K++;
                    k0++;
                }
            }
            else    {
                double ratio = 1.0 * k0 / K * (N+K-1) / (K-1);
                if (k0 == 1)    {
                    ratio *= 2;
                }
                if (rnd::GetRandom().Uniform() < ratio)  {
                    int k = K-1;
                    while ((k >= 0) && occupancy[k]) k--;
                    if (k < 0) {
                        cerr << "error when scanning for empty component\n";
                        exit(1);
                    }
                    K--;
                    k0--;
                    if (k != K) {
                        SwapComponents(k,K);
                    }
                }
            }
        }
    }

	if (K > Ncomponent)	{
		for (int k=Ncomponent; k<K; k++)	{
			CreateComponent(k);
		}
	}
	else if (K < Ncomponent)	{
		for (int k=K; k<Ncomponent; k++)	{
			DeleteComponent(k);
		}
	}
	Ncomponent = K;
    ResampleEmptyProfiles();
    ResampleWeights();
    return 1.0;
}

/* old version
double FiniteProfileProcess::MoveNcomponent(int nrep)	{

	int Kmax = GetNmodeMax();

	int kmax = Ncomponent-1;
	int kmin = 0;
	do	{
		while (!occupancy[kmax])	{
			kmax--;
		}
		while (occupancy[kmin])	{
			kmin++;
		}
		if (kmin < kmax)	{
			SwapComponents(kmin,kmax);
		}
	} while (kmin < kmax);
	int Kmin = kmin;
	if (Kmin > Kmax)	{
		Kmin--;
	}
	int K = Ncomponent;
	int nacc = 0;

	for (int rep=0; rep<nrep; rep++)	{
		double ratio = 1;
		int BK = K;
		if (K == Kmin)	{
			ratio *= 0.5 * ((double) K+1);
			ratio *= (K - Kmin + 1);
			ratio /= (GetNsite() + K * weightalpha);
			ratio *= K * weightalpha;
			K++;
		}
		else	{
			if (rnd::GetRandom().Uniform() < 0.5)	{
				ratio *= ((double) (K+1)) / (K - Kmin + 1);
				ratio /= (GetNsite() + K * weightalpha);
				ratio *= K * weightalpha;
				K++;
			}
			else	{
				if (K == Kmin + 1)	{
					ratio *= 2.0 / K;
					ratio *= (GetNsite() + (K - 1) * weightalpha);
					ratio /= (K-1) * weightalpha;
					K--;
				}
				else	{
					ratio *= ((double) (K - Kmin)) / K;
					ratio *= (GetNsite() + (K - 1) * weightalpha);
					ratio /= (K-1) * weightalpha;
					K--;
				}
			}
		}
		bool accept = ((rnd::GetRandom().Uniform() < ratio) && (K<=Kmax));
		if (! accept)	{
			K = BK;
		}
		else	{
			nacc ++;
		}
	}
	if (K > Ncomponent)	{
		for (int k=Ncomponent; k<K; k++)	{
			CreateComponent(k);
		}
	}
	else if (K < Ncomponent)	{
		for (int k=K; k<Ncomponent; k++)	{
			DeleteComponent(k);
		}
	}
	Ncomponent = K;
	return ((double) nacc) / nrep;
}
*/


void FiniteProfileProcess::ReadNcomponent(string filename)	{
	mixtype = filename;
	if ((filename == "WLSR5") || (filename == "wlsr5"))	{
		Ncat = WLSR5N;
	}
	else if ((filename == "CG6") || (filename == "cg6") || (filename == "c6") || (filename == "C6"))	{
		Ncat = 6;
	}
	else if ((filename == "CG10") || (filename == "cg10"))	{
		Ncat = 10;
	}
	else if ((filename == "CG20") || (filename == "cg20"))	{
		Ncat = 20;
	}
	else if ((filename == "CG30") || (filename == "cg30"))	{
		Ncat = 30;
	}
	else if ((filename == "CG40") || (filename == "cg40"))	{
		Ncat = 40;
	}
	else if ((filename == "CG50") || (filename == "cg50"))	{
		Ncat = 50;
	}
	else if ((filename == "CG60") || (filename == "cg60"))	{
		Ncat = 60;
	}
	else if ((filename == "c10") || (filename == "C10"))	{
		Ncat = C10N;
	}
	else if ((filename == "c20") || (filename == "C20"))	{
		Ncat = C20N;
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		Ncat = C30N;
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		Ncat = C30N;
	}
	else if ((filename == "c40") || (filename == "C40"))	{
		Ncat = C40N;
	}
	else if ((filename == "c50") || (filename == "C50"))	{
		Ncat = C50N;
	}
	else if ((filename == "c60") || (filename == "C60"))	{
		Ncat = C60N;
	}
	else if ((filename == "lg") || (filename == "LG"))	{
		Ncat = 1;
	}
	else if ((filename == "uniform") || (filename == "Uniform"))	{
		Ncat = 1;
	}
	else	{
		ifstream is(filename.c_str());
		if (!is)	{
			cerr << "error : unrecognized file for empirical mixture : " << filename << '\n';
			exit(1);
		}
		int tmpNstate;
		is >> tmpNstate;
		for (int k=0; k<tmpNstate; k++)	{
			string c;
			is >> c;
		}
		is >> Ncat;
	}
}

void FiniteProfileProcess::ReadStatFix(string filename)	{
	mixtype = filename;
	int Nstate = GetDim();
	if ((filename == "WLSR5") || (filename == "wlsr5"))	{
		if (Nstate != 20)	{
			cerr << "error: WLRS5 is for aminoacids\n";
		}
		int ncat = WLSR5N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat-1; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = WLSR5StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
		}
		statfix[Ncat-1] = new double[Nstate];	
		double* tmp = GetEmpiricalFreq();
		for (int k=0; k<Nstate; k++)	{
			statfix[Ncat-1][k] = tmp[k];
		}
		for (int i=0; i<Ncat; i++)	{
			empweight[i] = 1.0 / WLSR5N;
		}
	}
	else if ((filename == "CG6") || (filename == "cg6") || (filename == "c6") || (filename == "C6"))	{
		if (Nstate != 20)	{
			cerr << "error: CG6 is for aminoacids\n";
		}
		int ncat = 6;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG6StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "CG10") || (filename == "cg10"))	{
		if (Nstate != 20)	{
			cerr << "error: CG10 is for aminoacids\n";
		}
		int ncat = 10;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG10StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "CG20") || (filename == "cg20"))	{
		if (Nstate != 20)	{
			cerr << "error: CG20 is for aminoacids\n";
		}
		int ncat = 20;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG20StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "CG30") || (filename == "cg30"))	{
		if (Nstate != 20)	{
			cerr << "error: CG30 is for aminoacids\n";
		}
		int ncat = 30;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG30StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "CG40") || (filename == "cg40"))	{
		if (Nstate != 20)	{
			cerr << "error: CG40is for aminoacids\n";
		}
		int ncat = 40;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG40StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "CG50") || (filename == "cg50"))	{
		if (Nstate != 20)	{
			cerr << "error: CG50 is for aminoacids\n";
		}
		int ncat = 50;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG50StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "CG60") || (filename == "cg60"))	{
		if (Nstate != 20)	{
			cerr << "error: CG60 is for aminoacids\n";
		}
		int ncat = 60;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG60StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "c10") || (filename == "C10"))	{
		if (Nstate != 20)	{
			cerr << "error: C10 is for aminoacids\n";
		}
		int ncat = C10N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C10StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C10StatWeight[i];
		}
	}
	else if ((filename == "c20") || (filename == "C20"))	{
		if (Nstate != 20)	{
			cerr << "error: C20 is for aminoacids\n";
		}
		int ncat = C20N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C20StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C20StatWeight[i];
		}
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		if (Nstate != 20)	{
			cerr << "error: C30 is for aminoacids\n";
		}
		int ncat = C30N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C30StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C30StatWeight[i];
		}
	}
	else if ((filename == "c40") || (filename == "C40"))	{
		if (Nstate != 20)	{
			cerr << "error: C40 is for aminoacids\n";
		}
		int ncat = C40N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C40StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C40StatWeight[i];
		}
	}
	else if ((filename == "c50") || (filename == "C50"))	{
		if (Nstate != 20)	{
			cerr << "error: C50 is for aminoacids\n";
		}
		int ncat = C50N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C50StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C50StatWeight[i];
		}
	}
	else if ((filename == "c60") || (filename == "C60"))	{
		if (Nstate != 20)	{
			cerr << "error: C60 is for aminoacids\n";
		}
		int ncat = C60N;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C60StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C60StatWeight[i];
		}
	}
	else if ((filename == "lg") || (filename == "LG"))	{
		int ncat = 1;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = LG_Stat[k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0;
		}
	}
	else if ((filename == "uniform") || (filename == "Uniform"))	{
		int ncat = 1;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = 1.0 / Nstate;
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0;
		}
	}
	else	{
		ifstream is(filename.c_str());
		if (!is)	{
			cerr << "error : unrecognized file for empirical mixture : " << filename << '\n';
			exit(1);
		}
		int tmp;
		is >> tmp;
		if (tmp != Nstate)	{
			cerr << "error when reading empirical mixture : bad number of states\n";
			exit(1);
		}
		// read alphabet
		int permut[Nstate];
		for (int k=0; k<Nstate; k++)	{
			string c;
			is >> c;
			int l=0;
			while ((l<Nstate) && (c != GetStateSpace()->GetState(l))) l++;
			// while ((l<Nstate) && (c != AminoAcids[l])) l++;
			if (l == Nstate)	{
				cerr << "error when reading empirical mixture in " << filename << ": does not recognise letter " << c << '\n';
				cerr << "file should be formatted as follows\n";
				cerr << "list of amino acids (e.g. A C D E ...)\n";
				cerr << "Ncat\n";
				cerr << "empweight_1 (1 number) profile_1 (20 numbers)\n";
				cerr << "empweight_2 (1 number) profile_2 (20 numbers)\n";
				cerr << "...\n";
				cerr << "empweight_Ncat (1 number) profile_Ncat (20 numbers)\n";
				for (int a=0; a<Nstate; a++)	{
					cerr << "::" << GetStateSpace()->GetState(a) << "::";
				}
				cerr << '\n';
				exit(1); 
			}
			permut[k] = l;
		}
        int ncat;
		is >> ncat;
        if (ncat != Ncat)   {
            cerr << "error: incorrect number of components\n";
            cerr << "read: " << ncat << "; allocated: " << Ncat << '\n';
            exit(1);
        }
		statfix = new double*[Ncat];
		empweight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			is >> empweight[i];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				is >> statfix[i][permut[k]];
				if (statfix[i][permut[k]]<stateps)	{
					statfix[i][permut[k]] = stateps;
				}
				total += statfix[i][permut[k]];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
		}
		double total = 0;
		for (int i=0; i<Ncat; i++)	{
			total += empweight[i];
		}
		for (int i=0; i<Ncat; i++)	{
			empweight[i] /= total;
		}
	}
}

void FiniteProfileProcess::SetStatFix()	{

	if (! Ncat)	{
		cerr << "error in set stat fix\n";
		exit(1);
	}

	for (int k=0; k<Ncomponent; k++)	{
		DeleteComponent(k);
	}
	Ncomponent = Ncat;
	for (int k=0; k<Ncomponent; k++)	{
		CreateComponent(k);
	}
	
	for (int k=0; k<Ncomponent; k++)	{
		weight[k] = empweight[k];
		for (int i=0; i<GetDim(); i++)	{
			profile[k][i] = statfix[k][i];
		}
	}
	fixncomp = true;
	empmix = true;
}

