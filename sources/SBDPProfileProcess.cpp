
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "SBDPProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* SBDPProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void SBDPProfileProcess::Create(int innsite, int indim)	{

	if (! V)	{
		DPProfileProcess::Create(innsite,indim);
		V = new double[GetNmodeMax()];
		weight = new double[GetNmodeMax()];
	}
}

void SBDPProfileProcess::Delete()	{

	if (V)	{
		delete[] V;
		delete[] weight;
		DPProfileProcess::Delete();
	}
}

void SBDPProfileProcess::SampleAlloc()	{

	for (int k=0; k<GetNmodeMax(); k++)	{
		CreateComponent(k);
	}
	Ncomponent = GetNmodeMax();

	SampleWeights();
	for (int i=0; i<GetNsite(); i++)	{
		double U = rnd::GetRandom().Uniform();
		double total = weight[0];
		int k = 0;
		while ((k<GetNmodeMax()) && (total < U))	{
			k++;
			total += weight[k];
		}
		if (k == GetNmodeMax())	{
			cerr << "error in SBDPProfileProcess::SampleAlloc: overflow\n";
			exit(1);
		}
		AddSite(i,k);
	}
}

/*
void SBDPProfileProcess::IncrementalSampleAlloc()	{

	kappa = 0.1;
	int successful = 0;

	while (! successful)	{

		successful = 1;

		for (int i=0; i<GetNsite(); i++)	{
			RemoveSite(i,alloc[i]);
		}

		AddSite(0,0);
		Ncomponent = 1;
		
		for (int i=0; i<GetNsite(); i++)	{

			if (!successful)	{
				AddSite(i,0);
			}
			else	{
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
					Ncomponent++;
					if (Ncomponent == GetNmodeMax())	{
						cerr << "not successful : " << kappa << '\t' << Ncomponent << '\n';
						successful = 0;
					}
				}
				AddSite(i,k);
				delete[] p;
			}
		}
		if (! successful)	{
			kappa /= 10;
		}
		else	{
			cerr << "successful : " << kappa << '\t' << Ncomponent << '\n';
		}
	}

	Ncomponent = GetNmodeMax();
	ResampleWeights();
	cerr << "init incremental ok\n";
}
*/

void SBDPProfileProcess::IncrementalSampleAlloc()	{

	kappa = 0.1;

	for (int i=0; i<GetNsite(); i++)	{
		RemoveSite(i,alloc[i]);
	}

	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{

		int K = Ncomponent + 1;
		if (K > GetNmodeMax())	{
			K--;
		}
		double* p = new double[K];
		double total = 0;
		double max = 0;
		for (int k=0; k<K; k++)	{
			double w = occupancy[k];
			if (! w)	{
				w = kappa;
			}
			double tmp = log(w) * LogProxy(i,k);
			if ((!k) || (max < tmp))	{
				max = tmp;
			}
			p[k] = tmp;
		}
		for (int k=0; k<K; k++)	{
			double tmp = exp(p[k] - max);
			total += tmp;
			p[k] = total;
		}
		double q = total * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<K) && (q > p[k])) k++;
		if (k == K)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		if (k==Ncomponent)	{
			if (Ncomponent <= GetNmodeMax())	{
				Ncomponent++;
			}
		}
		AddSite(i,k);
		delete[] p;
	}

	Ncomponent = GetNmodeMax();
	ResampleWeights();
	cerr << "init incremental ok\n";
}

void SBDPProfileProcess::SampleWeights()	{

	double cumulProduct = 1.0;
	double totweight = 0;
	double v, x, y;
	for (int k=0; k<GetNcomponent(); k++)	{
		x = rnd::GetRandom().sGamma(1.0);
		y = rnd::GetRandom().sGamma(kappa);
		v = x / (x+y);
		V[k] = v;
		if (k == GetNcomponent() - 1)	{
			V[k] = 1;
			v = 1;
		}
		weight[k] = v * cumulProduct;
		cumulProduct *= (1 - v);	
		totweight += weight[k];
	}
}

void SBDPProfileProcess::ResampleWeights()	{

	UpdateOccupancyNumbers();
	// ???
	int remainingOcc = GetNsite();
	double cumulProduct = 1.0;
	double totweight = 0;
	double v, x, y;
	for (int k=0; k<GetNcomponent(); k++)	{
		remainingOcc -= occupancy[k];
		x = rnd::GetRandom().sGamma(1 + occupancy[k]);
		y = rnd::GetRandom().sGamma(kappa + remainingOcc);
		v = x / (x+y);
		V[k] = v;
		if (k == GetNcomponent() - 1)	{
			double tmp = cumulProduct * (1 - v);
			if (maxweighterror < tmp)	{
				maxweighterror = tmp;
			}
			V[k] = 1;
			v = 1;
		}
		weight[k] = v * cumulProduct;
		cumulProduct *= (1 - v);
		totweight += weight[k];
	}
}

double SBDPProfileProcess::MoveOccupiedCompAlloc(int k0)	{

	int nrep = (int) (k0 * kappa);
	UpdateOccupancyNumbers();
	ResampleWeights();
	double total = 0.0;
	int Nocc = GetNOccupiedComponent();
	if (Nocc != 1)	{
		for (int i=0; i<nrep; i++)	{
			int* occupiedComponentIndices = new int[Nocc];
			int j=0;
			for (int k=0; k<GetNcomponent(); k++)	{
				if (occupancy[k] != 0)	{
					occupiedComponentIndices[j] = k;
					j++;
				}
			}
			if (j != Nocc)	{
				cerr << "error in MoveOccupiedCompAlloc.\n";
				exit(1);
			}
			int* indices = new int[2];
			rnd::GetRandom().DrawFromUrn(indices,2,Nocc);
			int cat1 = occupiedComponentIndices[indices[0]];
			int cat2 = occupiedComponentIndices[indices[1]];
			double logMetropolis = (occupancy[cat2] - occupancy[cat1]) * log(weight[cat1] / weight[cat2]);
			int accepted = (log(rnd::GetRandom().Uniform()) < logMetropolis);
			if (accepted)	{
				total += 1.0;
				// SwapComponents(cat1, cat2);
				MixtureProfileProcess::SwapComponents(cat1,cat2);
			
			}
			delete[] occupiedComponentIndices;
			delete[] indices; 
		}
		return total /= nrep;
	}
	return 0;
}

double SBDPProfileProcess::LogIntegratedAllocProb()	{
	int remainingOcc = GetNsite();
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		if (remainingOcc)	{
			remainingOcc -= occupancy[k];
			total += log(kappa) + rnd::GetRandom().logGamma(1 + occupancy[k]) + rnd::GetRandom().logGamma(kappa + remainingOcc) - rnd::GetRandom().logGamma(1 + kappa + occupancy[k] + remainingOcc);
		}
	}
	if (remainingOcc)	{
		cerr << "error in allocation count\n";
		exit(1);
	}
	return total;
}

double SBDPProfileProcess::MoveKappa(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogIntegratedAllocProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		kappa *= e;
		deltalogprob += LogHyperPrior() + LogIntegratedAllocProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			kappa /= e;
		}
	}
	// not useful
	// ResampleWeights();
	return naccepted / nrep;
}

double SBDPProfileProcess::LogStatPrior()	{

	UpdateOccupancyNumbers();
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		if (occupancy[i])	{
			total += DPProfileProcess::LogStatPrior(i);
		}
	}
	return total;
}

void SBDPProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	double tempv = V[cat1];
	V[cat1] = V[cat2];
	V[cat2] = tempv;
	double tempw = weight[cat1];
	weight[cat1] = weight[cat2];
	weight[cat2] = tempw;
}



double SBDPProfileProcess::MoveAdjacentCompAlloc(int k0)	{

	int nrep = (int) (k0 * kappa);
	ResampleWeights();
	
	double total = 0;

	for (int i=0; i<nrep; i++)	{
		//int cat1 = (int)(rnd::GetRandom().Uniform() * (GetNcomponent()-1));
		int cat1 = (int)(rnd::GetRandom().Uniform() * (GetNcomponent()-2));  
		int cat2 = cat1 + 1;
		double logMetropolis = (occupancy[cat1] * log(1 - V[cat2])) - (occupancy[cat2] * log(1-V[cat1]));
		int accepted = (log(rnd::GetRandom().Uniform()) < logMetropolis);
		if (accepted)	{
			total += 1.0;
			SwapComponents(cat1,cat2);
		}
	}

	return total /= nrep;
}


