
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonSubstitutionProcess.h"

#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Poisson Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// same methods as above, but specialized for the CAT poisson model
// probably less interesting to parallelize, at least in a first step


//-------------------------------------------------------------------------
//	* conditional likelihood propagation
//	(CPU level 3)
//-------------------------------------------------------------------------

void PoissonSubstitutionProcess::Propagate(double*** from, double*** to, double time, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            const double* stat = GetStationary(i);
            for (int j=0; j<GetNrate(i); j++)	{
                if ((! condalloc) || (ratealloc[i] == j))	{
                    double* tmpfrom = from[i][j];
                    double* tmpto = to[i][j];
                    double expo = exp(-GetRate(i,j) * time);
                    double tot = 0;
                    int nstate = GetNstate(i);
                    for (int k=0; k<nstate; k++)	{
                        tot += (*tmpfrom++) * (*stat++);
                        // tot += tmpfrom[k] * stat[k];
                    }
                    tmpfrom -= nstate;
                    stat -= nstate;
                    tot *= (1-expo);
                    for (int k=0; k<nstate; k++)	{	
                        (*tmpto++) = expo * (*tmpfrom++) + tot;
                        // tmpto[k] = expo * tmpfrom[k] + tot;
                    }
                    (*tmpto) = (*tmpfrom);
                    tmpto -= nstate;
                    tmpfrom -= nstate;
                    // tmpto[GetNstate(i)] = tmpfrom[GetNstate(i)];
                }
            }
        }
    }
}

/*
// version directly unzipped
void PoissonSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{

	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetProfile(i);
		int j = ratealloc[i];
		double expo = exp(-GetRate(i,j) * time);
		if (rnd::GetRandom().Uniform() < expo)	{
			statedown[i] = stateup[i];
		}
		else	{
			double u = rnd::GetRandom().Uniform();
			int k = 0;
			double cumul = stat[k];
			while ((k < GetDim()) && (u > cumul))	{
				k++;
				if (k == GetDim())	{
					cerr << "error in PoissonSubstitutionProcess::SimuPropagate: overflow\n";
					exit(1);
				}
				cumul += stat[k];
			}
			statedown[i] = k;
		}
	}
}
*/

void PoissonSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{

	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		int nstate = GetNstate(i);
		int j = ratealloc[i];
		double expo = exp(-GetRate(i,j) * time);
		if (rnd::GetRandom().Uniform() < expo)	{
			statedown[i] = stateup[i];
		}
		else	{
			double u = rnd::GetRandom().Uniform();
			int k = 0;
			double cumul = stat[k];
			while ((k < nstate) && (u > cumul))	{
				k++;
				if (k == GetDim())	{
					cerr << "error in PoissonSubstitutionProcess::SimuPropagate: overflow\n";
					exit(1);
				}
				cumul += stat[k];
			}
			statedown[i] = k;
		}
	}
}

//-------------------------------------------------------------------------
//	* sample substitution mappings conditional on states at nodes 
//	(CPU level 3)
//-------------------------------------------------------------------------

// root version
BranchSitePath** PoissonSubstitutionProcess::SampleRootPaths(int* state)	{
	// BranchSitePath** patharray = new BranchSitePath*[sitemax - sitemin];
	BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		patharray[i] = new BranchSitePath(0,state[i]);
	}
	return patharray;
}

// general version
BranchSitePath** PoissonSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time) 	{
	BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		double rate = GetRate(i);
		double l = rate * time;
		int dup = stateup[i];
		int ddown = statedown[i];
		double pi = stat[ddown];

		int m = 0;
		int mmax = 1000;
		
		if (dup == ddown)	{
			double fact = pi * exp(-l);
			double total = exp(-l);
			double q = rnd::GetRandom().Uniform() * (exp(-l) * (1 - pi) + pi);
			while ((m<mmax) && (total < q))	{
				m++;
				fact *= l / m;
				total += fact;
			}
			if (m == mmax)	{
				suboverflowcount ++;
			}
		}
		else	{
			double fact = pi * exp(-l);
			double total = 0;
			double q = rnd::GetRandom().Uniform() * (1 - exp(-l)) * pi;
			while ((m<mmax) && (total < q))	{
				m++;
				fact *= l / m;
				total += fact;
			}
			if (m == mmax)	{
				suboverflowcount ++;
			}
		}
		patharray[i] = new BranchSitePath(m,ddown);
	}
	return patharray;
}

//-------------------------------------------------------------------------
//	* gather sufficient statistics 
//	(CPU level 3)
//-------------------------------------------------------------------------


void PoissonSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, double* siteratesuffstatbeta, double branchlength, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			siteratesuffstatcount[i] += patharray[i]->GetNsub();
			siteratesuffstatbeta[i] += branchlength;
		}
	}
}


void PoissonSubstitutionProcess::AddBranchLengthSuffStat(int& count, double& beta, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (nonmissing[i] == 1)	{
			count += patharray[i]->GetNsub();
			beta += GetRate(i);
		}
	}
}


void PoissonSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root)	{
	cerr << "error: in PoissonSubstitutionProcess::AddSiteProfileSuffStat: deprecated\n";
	exit(1);
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (root || patharray[i]->GetNsub())	{
			siteprofilesuffstatcount[i][GetRandomStateFromZip(i,patharray[i]->GetFinalState())]++;
		}
	}
}

/*
void PoissonSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		siteratesuffstatcount[i] += patharray[i]->GetNsub();
	}
}


void PoissonSubstitutionProcess::AddBranchLengthSuffStat(int& count, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		count += patharray[i]->GetNsub();
	}
}


void PoissonSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		if (root || patharray[i]->GetNsub())	{
			siteprofilesuffstatcount[i][GetRandomStateFromZip(i,patharray[i]->GetFinalState())]++;
		}
	}
}
*/

void PoissonSubstitutionProcess::ChooseTrueStates(BranchSitePath** patharray, int* nodestateup, int* nodestatedown, bool root)	{
	for (int i=sitemin; i<sitemax; i++)	{
		int tmp = nodestateup[i];
		if (root || patharray[i]->GetNsub())	{
			tmp = GetRandomStateFromZip(i,patharray[i]->GetFinalState());
		}
		//cerr  << ' '<< i  << ' '<< patharray[i]->GetNsub() << ' '<< tmp  << '\n';
		nodestatedown[i] = tmp;
	}
}

void PoissonSubstitutionProcess::ChooseRootTrueStates(int* nodestate)	{

	for (int i=sitemin; i<sitemax; i++)	{
		int tmp = GetRandomStateFromZip(i,nodestate[i]);
		nodestate[i] = tmp;
	}
}

void PoissonSubstitutionProcess::ConditionalLikelihoodsToStatePostProbs(double*** aux,double*** statepostprob, int nodelabel, bool condalloc)	{

	SubstitutionProcess::ConditionalLikelihoodsToStatePostProbs(aux,statepostprob,nodelabel,condalloc);
	ZipToTruePostProbs(statepostprob,nodelabel);
}

void PoissonSubstitutionProcess::ZipToTruePostProbs(double*** statepostprob, int nodelabel)	{

	double tmp[GetGlobalNstate()];
	for (int i=sitemin; i<sitemax; i++)	{
		double* s = statepostprob[i][nodelabel];
		double* pi = GetProfile(i);
	
		double total = 0;
		for (int k=0; k<GetOrbitSize(i); k++)	{
			tmp[GetStateFromZip(i,k)] = s[k];
			total += zipstat[i][k];
		}
		if (GetZipSize(i) > GetOrbitSize(i))	{
			double s0 = s[GetOrbitSize(i)];
			for (int k=0; k<GetGlobalNstate(); k++)	{
				if (! InOrbit(i,k))	{
					tmp[k] = s0 * pi[k] / (1 - total);
				}
			}
		}
		double totprob = 0;
		for (int k=0; k<GetGlobalNstate(); k++)	{
			s[k] = tmp[k];
			totprob += tmp[k];
		}
		if (fabs(totprob - 1) > 1e-7)	{
			cerr << "error in PoissonSubstitutionProcess::ZipToTruePostProbs: total is not 1\n";
			cerr << totprob << '\n';
			exit(1);
		}
	}
}

//-------------------------------------------------------------------------
//	* recomputing the equilibrium frequency profiles of the recoded process
//
//-------------------------------------------------------------------------


void PoissonSubstitutionProcess::CreateZip()	{
	if (! zipstat)	{
		zipstat = new double*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			zipstat[i] = new double[GetDim()];
			// zipstat[i] = new double[GetNstate(i)];
		}
	}
}

void PoissonSubstitutionProcess::DeleteZip()	{
	for (int i=0; i<GetNsite(); i++)	{
		delete[] zipstat[i];
	}
	delete[] zipstat;
	zipstat = 0;
}

void PoissonSubstitutionProcess::UpdateZip()	{
	
	for (int i=0; i<GetNsite(); i++)	{
		UpdateZip(i);
	}
}

void PoissonSubstitutionProcess::UpdateZip(int i)	{
		double total = 0;
		double* pi = GetProfile(i);
		for (int k=0; k<GetOrbitSize(i); k++)	{
			zipstat[i][k] = pi[GetStateFromZip(i,k)];
			total += zipstat[i][k];
		}
		if (GetZipSize(i) > GetOrbitSize(i))	{
			zipstat[i][GetOrbitSize(i)] = 1-total;
		}
}

int PoissonSubstitutionProcess::GetRandomStateFromZip(int site, int zipstate)	{
	int truestate = 0;
	if ((GetZipSize(site) != GetOrbitSize(site)) && (zipstate == GetOrbitSize(site)))	{
		double v = rnd::GetRandom().Uniform();
		double u = zipstat[site][GetOrbitSize(site)] * v;
		double total = 0;
		double* pi = GetProfile(site);

		/*
		int myid;
		MPI_Comm_rank(MPI_COMM_WORLD,&myid);
		if (myid == 4)	{
			cerr << myid << "??" << '\t' << site << '\t' << v << '\t' << u << '\t' << GetOrbitSize(site) << '\t' << zipstat[site][GetOrbitSize(site)] << '\n';
			for (int k=0; k<GetDim(); k++)	{
				cerr << pi[k] << '\t';
			}
			cerr << '\n';
		}
		*/

		int choose = -1;
		while ((choose < GetDim()) && (total < u))	{
			choose ++;
			if (choose == GetDim())	{
				cerr << "error in getstatefromzip\n";
				cerr << choose << '\t' << GetDim() << '\n';
				cerr << total << '\n';
				cerr << v << '\t' << zipstat[site][GetOrbitSize(site)] << '\t' << u << '\n';
				cerr << total - zipstat[site][GetOrbitSize(site)] << '\n';
				double newtotal = 0;
				cerr << '\n';
				for (int k=0; k<GetDim(); k++)	{
					cerr << pi[k] << '\t';
					cerr << InOrbit(site,k) << '\n';
					if (! InOrbit(site,k))	{
						newtotal += pi[k];
					}
				}
				cerr << '\n';
				for (int k=0; k<=GetOrbitSize(site); k++)	{
					cerr << zipstat[site][k] << '\n';
				}
				cerr << "new total : " << newtotal << '\t' << newtotal - total << '\t' << total - choose << '\n';
				exit(1);
			}
			if (!InOrbit(site,choose))	{
				total += pi[choose];
			}
		}
		truestate = choose;
	}
	else	{
		truestate = GetStateFromZip(site,zipstate);
	}
	return truestate;
}


void PoissonSubstitutionProcess::UnzipBranchSitePath(BranchSitePath** patharray, int* nodestateup, int* nodestatedown){
	for (int i=sitemin; i<sitemax; i++)	{
		int nsub = patharray[i]->GetNsub();
		patharray[i]->nsub=0;
		double* times = new double[nsub+1];
		for(int j = 0; j < nsub; j++){
			times[j] = rnd::GetRandom().Uniform();
			for(int k = 0; k < j; k++){
				if(times[k]>times[j]	){
					times[nsub] = times[k];
					times[k] = times[j];
					times[j] = times[nsub];
				}
			}
		}
		double mem = 0;
		for(int j = 0; j < nsub; j++){
			times[nsub] = times[j];
			times[j] = times[j] - mem;
			mem = times[nsub];
		}
		times[nsub]=1-mem;

		int previousstate = nodestateup[i];
		patharray[i]->Init()->SetState(previousstate);
		double* pi = GetProfile(i);
		for(int j = 0; j < nsub-1; j++){
			int newstate = rnd::GetRandom().DrawFromDiscreteDistribution(pi, GetDim());
			if(newstate != previousstate){
			      patharray[i]->Append(newstate, times[j]);
			      previousstate=newstate;
			}
			else{
			      times[j+1]+=times[j];
			}
		}
		if(previousstate == nodestatedown[i]){
			times[nsub] += times[nsub-1];
		}
		else{
			patharray[i]->Append(nodestatedown[i], times[nsub-1]);
		}
		patharray[i]->Last()->SetRelativeTime(times[nsub]);
		delete[] times;
	}
}


void PoissonSubstitutionProcess::SitePropagate(int i, double** from, double** to, double time, bool condalloc)	{

	const double* stat = GetStationary(i);
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmpfrom = from[j];
			double* tmpto = to[j];
			double expo = exp(-GetRate(i,j) * time);
			double tot = 0;
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				tot += (*tmpfrom++) * (*stat++);
				// tot += tmpfrom[k] * stat[k];
			}
			tmpfrom -= nstate;
			stat -= nstate;
			tot *= (1-expo);
			for (int k=0; k<nstate; k++)	{	
				(*tmpto++) = expo * (*tmpfrom++) + tot;
				// tmpto[k] = expo * tmpfrom[k] + tot;
			}
			(*tmpto) = (*tmpfrom);
			tmpto -= nstate;
			tmpfrom -= nstate;
			// tmpto[GetNstate(i)] = tmpfrom[GetNstate(i)];
		}
	}
}

