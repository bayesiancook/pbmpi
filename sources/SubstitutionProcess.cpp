
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "SubstitutionProcess.h"
#include "Random.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* allocations / deallocations
//-------------------------------------------------------------------------

void SubstitutionProcess::Create(int site, int dim, int insitemin, int insitemax)	{
	sitemin = insitemin;
	sitemax = insitemax;
	if (! ratealloc)	{
		RateProcess::Create(site);
		ProfileProcess::Create(site, dim);
		// ratealloc = new int[sitemax - sitemin];
		ratealloc = new int[GetNsite()];
	}
}

void SubstitutionProcess::Delete() {
	if (ratealloc)	{
		delete[] ratealloc;
		ratealloc = 0;
		ProfileProcess::Delete();
		RateProcess::Delete();
	}
};

void SubstitutionProcess::CreateCondSiteLogL()	{
	if (condsitelogL)	{
		cerr << "error in SubstitutionProcess::CreateSiteLogL\n";
		exit(1);
	}
	// sitelogL = new double[sitemax - sitemin];
	// condsitelogL = new double*[sitemax - sitemin];
	sitelogL = new double[GetNsite()];
	meansiterate = new double[GetNsite()];
	condsitelogL = new double*[GetNsite()];
	// for (int i=sitemin; i<sitemax; i++)	{
	for (int i=0; i<GetNsite(); i++)	{
		condsitelogL[i] = new double[GetNrate(i)];
	}
}

void SubstitutionProcess::DeleteCondSiteLogL()	{
	if (condsitelogL)	{
		// for (int i=sitemin; i<sitemax; i++)	{
		for (int i=0; i<GetNsite(); i++)	{
			delete[] condsitelogL[i];
		}
		delete[] condsitelogL;
		delete[] sitelogL;
		delete[] meansiterate;
		condsitelogL = 0;
		sitelogL = 0;
	}
}

double*** SubstitutionProcess::CreateConditionalLikelihoodVector()	{
	//cout << "VECTOR ALLOCATION: " << sitemax << "  " << sitemin << endl;
	// double*** condl = new double**[sitemax - sitemin];
	double*** condl = new double**[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		condl[i] = new double*[GetNrate(i)];
		for (int j=0; j<GetNrate(i); j++)	{
			condl[i][j] = new double[GetNstate(i) + 1];
			double* tmp = condl[i][j];
			for (int k=0; k<GetNstate(i); k++)	{
				tmp[k] = 1.0;
			}
			tmp[GetNstate(i)] = 0;
		}
	}
	//cout << "Test element " << condl[0][1][0] << endl;
	return condl;
}

void SubstitutionProcess::DeleteConditionalLikelihoodVector(double*** condl)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			delete[] condl[i][j];
		}
		delete[] condl[i];
	}
	delete[] condl;
	condl = 0;
}

//-------------------------------------------------------------------------
//	* elementary computations on conditional likelihood vectors 
//	(CPU level 1)
//-------------------------------------------------------------------------

// set the vector uniformly to 1 
void SubstitutionProcess::Reset(double*** t, bool condalloc, bool all)	{
	for (int i=sitemin; i<sitemax; i++)	{
        // if (ActiveSite(i))  {
        if (all || ActiveSite(i))  {
            for (int j=0; j<GetNrate(i); j++)	{
                if ((! condalloc) || (ratealloc[i] == j))	{
                    double* tmp = t[i][j];
                    int nstate = GetNstate(i);
                    for (int k=0; k<nstate; k++)	{
                        (*tmp++) = 1.0;
                        // tmp[k] = 1.0;
                    }
                    *tmp = 0;
                    tmp -= nstate;
                    // tmp[GetNstate(i)] = 0;
                }
            }
        }
	}
}
	
// initialize the vector according to the data observed at a given leaf of the tree (contained in const int* state)
// steta[i] == -1 means 'missing data'. in that case, conditional likelihoods are all 1
void SubstitutionProcess::Initialize(double*** t, const int* state, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            for (int j=0; j<GetNrate(i); j++)	{
                if ((! condalloc) || (ratealloc[i] == j))	{
                    double* tmp = t[i][j];
                    int nstate = GetNstate(i);
                    tmp[nstate] = 0;
                    if (state[i] == -1)	{
                        for (int k=0; k<nstate; k++)	{
                            (*tmp++) = 1.0;
                            // tmp[k] = 1.0;
                        }
                        tmp -= nstate;
                    }
                    else	{
                        for (int k=0; k<nstate; k++)	{
                            (*tmp++) = 0;
                            // tmp[k] = 0;
                        }
                        tmp -= nstate;
                        tmp[state[i]] = 1.0;
                    }
                }
            }
        }
    }
}

// multiply two conditional likelihood vectors, term by term
void SubstitutionProcess::Multiply(double*** from, double*** to, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            for (int j=0; j<GetNrate(i); j++)	{
                if ((! condalloc) || (ratealloc[i] == j))	{
                    double* tmpfrom = from[i][j];
                    double* tmpto = to[i][j];
                    int nstate = GetNstate(i);
                    for (int k=0; k<nstate; k++)	{
                        (*tmpto++) *= (*tmpfrom++);
                        // tmpto[k] *= tmpfrom[k];
                    }
                    *tmpto += *tmpfrom;
                    tmpto -= nstate;
                    tmpfrom -= nstate;
                    // tmpto[GetNstate(i)] += tmpfrom[GetNstate(i)];
                }
            }
        }
    }
}

// multiply a conditional likelihood vector by the (possibly site-specific) stationary probabilities of the process
void SubstitutionProcess::MultiplyByStationaries(double*** to, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            const double* stat = GetStationary(i);
            for (int j=0; j<GetNrate(i); j++)	{
                if ((! condalloc) || (ratealloc[i] == j))	{
                    double* tmpto = to[i][j];
                    int nstate = GetNstate(i);
                    for (int k=0; k<nstate; k++)	{	
                        (*tmpto++) *= (*stat++);
                        // tmpto[k] *= stat[k];
                    }
                    tmpto -= nstate;
                    stat -= nstate;
                }
            }
        }
    }
}

// to avoid numerical errors: all entries for a given site and a given rate
// are divided by the largest among them
// and the residual is stored in the last entry of the vector
void SubstitutionProcess::Offset(double*** t, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            for (int j=0; j<GetNrate(i); j++)	{
                if ((! condalloc) || (ratealloc[i] == j))	{
                    double* tmp = t[i][j];
                    double max = 0;
                    for (int k=0; k<GetNstate(i); k++)	{
                        if (tmp[k] <0)	{
                            cerr << "error in pruning: negative prob : " << tmp[k] << "\n";
                            exit(1);
                            // tmp[k] = 0;
                        }
                        if (max < tmp[k])	{
                            max = tmp[k];
                        }
                    }
                    /*
                    if (max == 0)	{
                        max = 1e-12;
                    }
                    */
                    if (max < 0)	{
                        cerr << "error in pruning (offset function): null likelihood\n";
                        exit(1);
                    }
                    if (max > 0)	{
                        for (int k=0; k<GetNstate(i); k++)	{
                            tmp[k] /= max;
                        }
                    }
                    tmp[GetNstate(i)] += log(max);
                }
            }
        }
    }
}

//-------------------------------------------------------------------------
//	* likelihood computation (last step, once the recursion has proceeded throughout the entire tree) 
//	(CPU level 2)
//-------------------------------------------------------------------------

double SubstitutionProcess::ComputeLikelihood(double*** aux, bool condalloc)	{

	for (int i=sitemin; i<sitemax; i++)	{
        sitelogL[i] = 0;
        for (int j=0; j<GetNrate(i); j++)	{
            condsitelogL[i][j] = 0;
        }
    }

	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            if (condalloc)	{
                int j = ratealloc[i];
                double* t = aux[i][j];
                double tot = 0;
                int nstate = GetNstate(i);
                for (int k=0; k<nstate; k++)	{
                    tot += (*t++);
                }
                if (tot == 0)	{
                    // dirty !
                    tot = 1e-12;
                }
                sitelogL[i] = log(tot) + (*t);
                t -= nstate;
            }
            else	{
                double max = 0;
                double* logl = condsitelogL[i];
                for (int j=0; j<GetNrate(i); j++)	{
                    double* t = aux[i][j];
                    double tot = 0;
                    int nstate = GetNstate(i);
                    for (int k=0; k<nstate; k++)	{
                        tot += (*t++);
                    }
                    if (tot < 0)	{
                        cerr << "error in SubstitutionProcess::ComputeLikelihood: negative prob\n";
                        cerr << tot << '\n';
                        exit(1);
                    }
                    if (std::isnan(tot))	{
                        cerr << "error in SubstitutionProcess::ComputeLikelihood: tot is nan\n";
                        for (int k=0; k<nstate; k++)	{
                            cerr << GetStationary(i)[k] << '\t';
                        }
                        cerr << '\n';
                        exit(1);
                    }
                    if (tot == 0)	{
                        // dirty !
                        tot = 1e-12;
                    }
                    logl[j] = log(tot) + (*t);
                    t -= nstate;
                    if ((!j) || (max < logl[j]))	{
                        max = logl[j];
                    }
                }
                double total = 0;
                double meanrate = 0;
                if (std::isinf(max))	{
                    meanrate = 1.0;
                    sitelogL[i] = max;
                }
                else	{
                    for (int j=0; j<GetNrate(i); j++)	{
                        double tmp = GetRateWeight(i,j) * exp(logl[j] - max);
                        total += tmp;
                        meanrate += tmp * GetRate(i,j);
                    }
                    sitelogL[i] = log(total) + max;
                }
                if (std::isnan(sitelogL[i]))	{
                    cerr << "error in SubstitutionProcess::ComputeNodeLikelihood: nan\n";
                    cerr << total << '\t' << max << '\n';
                    for (int j=0; j<GetNrate(i); j++)	{
                        cerr << logl[j] << '\t';
                    }
                    cerr << '\n';
                    exit(1);
                }
                meanrate /= total;
                meansiterate[i] = meanrate;
            }
        }

    }
    logL = 0;
    for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            logL += sitelogL[i];
        }
    }
    return logL;
}
	

void SubstitutionProcess::ConditionalLikelihoodsToStatePostProbs(double*** aux,double*** statepostprob, int nodelabel, bool condalloc)	{

	for (int i=sitemin; i<sitemax; i++)	{
        if (condalloc)	{
            int j = ratealloc[i];
            double* t = aux[i][j];
            double* s = statepostprob[i][nodelabel];
            /*
            if (GetNstate(i) != GetNstate())	{
                cerr << "error in SubstitutionProcess::ConditionalLikelihoodsToStatePostProbs: non matching number of states\n";
                exit(1);
            }
            */
            double total = 0;
            for (int k=0; k<GetNstate(i); k++)	{
                s[k] = t[k];
                total += s[k];
            }
            for (int k=0; k<GetNstate(i); k++)	{
                s[k] /= total;
            }
        }
        else	{
            double* s = statepostprob[i][nodelabel];
            for (int k=0; k<GetNstate(i); k++)	{
                s[k] =0;
            }
            for (int j=0; j<GetNrate(i); j++)	{
                double* t = aux[i][j];
                for (int k=0; k<GetNstate(i); k++)	{
                    s[k] += t[k];
                }
            }
            double total = 0;
            for (int k=0; k<GetNstate(i); k++)	{
                total += s[k];
            }
            for (int k=0; k<GetNstate(i); k++)	{
                s[k] /= total;
            }
        }
    }
}

//-------------------------------------------------------------------------
//	* sample the allocation of each site to one of the available rate categories
// 	for each site, the rate caegory is chosen with probability proportional
//	to the likelihood for that site, under that rate category
//	is aux == 0, this likelihood is assumed to be already computed and stored in int** condsitelogL (a member of SubstitutionProcess)
//	otherwise, aux is assumed to be a conditional likelihood vector multiplied by stationaries, and thus can be used to compute those likelihoods 
//	(CPU level 2)
//-------------------------------------------------------------------------

void SubstitutionProcess::DrawAllocations(double*** aux)	{

	if (aux)	{
		ComputeLikelihood(aux);
	}
	for (int i=sitemin; i<sitemax; i++)	{
        if (GetNrate(i) == 1)	{
            ratealloc[i] = 0;
        }
        else	{
            double max = 0;
            double* logl = condsitelogL[i];
            for (int j=0; j<GetNrate(i); j++)	{
                if ((!j) || (max < logl[j]))	{
                    max = logl[j];
                }
            }
            double cumul = 0;
            double p[GetNrate(i)];
            for (int j=0; j<GetNrate(i); j++)	{
                cumul += GetRateWeight(i,j) * exp(logl[j] - max);
                p[j] = cumul;
            }
            double u = rnd::GetRandom().Uniform() * cumul;
            int j = 0;
            while ((j<GetNrate(i)) && (p[j] < u)) j++;
            if (j == GetNrate(i))	{
                cerr << "error in SubstitutionProcess::SampleAlloc\n";
                exit(1);
            }
            ratealloc[i] = j;
        }
    }
}

void SubstitutionProcess::DrawAllocationsFromPrior()	{

	for (int i=sitemin; i<sitemax; i++)	{
        int k = (int) (GetNrate(i) * rnd::GetRandom().Uniform());
        ratealloc[i] = k;
	}
}

//-------------------------------------------------------------------------
//	* sample states for each site, based on vectors of posterior probability stored in double*** t
//	this is assumed to be conditional on rate allocations 
//	(CPU level 2)
//-------------------------------------------------------------------------

void SubstitutionProcess::ChooseStates(double*** t, int* states)	{
	for (int i=sitemin; i<sitemax; i++)	{
		int j = ratealloc[i];
		double* tmp = t[i][j];
		double total = 0;
		for (int k=0; k<GetNstate(i); k++)	{
			total += tmp[k];
		}
		double u = rnd::GetRandom().Uniform() * total;
		double tot = tmp[0];
		int k = 0;
		while ((k<GetNstate(i)) && (tot < u))	{
			k++;
			if (k==GetNstate(i))	{
				cerr << "error in SubstitutionProcess::ChooseState\n";
				exit(1);
			}
			tot += tmp[k];
		}
		states[i] = k;
		for (int l=0; l<GetNstate(i); l++)	{
			tmp[l] = 0;
		}
		tmp[k] =  1;
	}
}

void SubstitutionProcess::ChooseStatesAtEquilibrium(int* states)	{
	
	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		int nstate = GetNstate(i);
		double cumul[nstate];
		double tot = 0;
		for (int k=0; k<nstate; k++)	{
			tot += stat[k];
			cumul[k] = tot;
		}
		if (fabs(tot-1) > 1e-7)	{
			cerr << "error in SubstitutionProcess::ChooseStatesAtEquilibrium: does not sum to 1: " << tot << '\n';
			exit(1);
		}
		double u = rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<nstate) && (u > cumul[k]))	{
			k++;
		}
		if (k == nstate)	{
			cerr << "error in SubstitutionProcess::ChooseStatesAtEquilibrium: overflow\n";
			exit(1);
		}
		states[i] = k;
	}
}

void SubstitutionProcess::SetCondToStates(double*** t, int* states)	{
	for (int i=sitemin; i<sitemax; i++)	{
		int j = ratealloc[i];
		double* tmp = t[i][j];
		for (int l=0; l<GetNstate(i); l++)	{
			tmp[l] = 0;
		}
		if ((states[i] < 0) || (states[i] >= GetNstate(i)))	{
			cerr << "error: state out of bound\n";
			cerr << states[i] << '\n';
			exit(1);
		}
		tmp[states[i]] =  1;
	}
}


// single-site versions

//-------------------------------------------------------------------------
//	* elementary computations on conditional likelihood vectors 
//	(CPU level 1)
//-------------------------------------------------------------------------

// set the vector uniformly to 1 
void SubstitutionProcess::SiteReset(int i, double** t, bool condalloc)	{
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmp = t[j];
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				(*tmp++) = 1.0;
			}
			*tmp = 0;
			tmp -= nstate;
		}
	}
}
	
// initialize the vector according to the data observed at a given leaf of the tree (contained in const int* state)
// steta[i] == -1 means 'missing data'. in that case, conditional likelihoods are all 1
void SubstitutionProcess::SiteInitialize(int i, double** t, const int state, bool condalloc)	{

	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmp = t[j];
			int nstate = GetNstate(i);
			tmp[nstate] = 0;
			if (state == -1)	{
				for (int k=0; k<nstate; k++)	{
					(*tmp++) = 1.0;
				}
				tmp -= nstate;
			}
			else	{
				for (int k=0; k<nstate; k++)	{
					(*tmp++) = 0;
				}
				tmp -= nstate;
				tmp[state] = 1.0;
			}
		}
	}
}

// multiply two conditional likelihood vectors, term by term
void SubstitutionProcess::SiteMultiply(int i, double** from, double** to, bool condalloc)	{
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmpfrom = from[j];
			double* tmpto = to[j];
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				(*tmpto++) *= (*tmpfrom++);
			}
			*tmpto += *tmpfrom;
			tmpto -= nstate;
			tmpfrom -= nstate;
		}
	}
}

// multiply a conditional likelihood vector by the (possibly site-specific) stationary probabilities of the process
void SubstitutionProcess::SiteMultiplyByStationaries(int i, double** to, bool condalloc)	{
	const double* stat = GetStationary(i);
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmpto = to[j];
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{	
				(*tmpto++) *= (*stat++);
			}
			tmpto -= nstate;
			stat -= nstate;
		}
	}
}

// to avoid numerical errors: all entries for a given site and a given rate
// are divided by the largest among them
// and the residual is stored in the last entry of the vector
void SubstitutionProcess::SiteOffset(int i, double** t, bool condalloc)	{
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmp = t[j];
			double max = 0;
			for (int k=0; k<GetNstate(i); k++)	{
				if (tmp[k] <0)	{
					cerr << "error in pruning: negative prob : " << tmp[k] << "\n";
					exit(1);
					tmp[k] = 0;
				}
				if (max < tmp[k])	{
					max = tmp[k];
				}
			}
			if (max > 0)	{
				for (int k=0; k<GetNstate(i); k++)	{
					tmp[k] /= max;
				}
			}
			tmp[GetNstate(i)] += log(max);
		}
	}
}

//-------------------------------------------------------------------------
//	* likelihood computation (last step, once the recursion has proceeded throughout the entire tree) 
//	(CPU level 2)
//-------------------------------------------------------------------------

double SubstitutionProcess::SiteComputeLikelihood(int i, double** aux, bool condalloc)	{

    double sitelogl = 0;
	if (condalloc)	{
		int j = ratealloc[i];
		double* t = aux[j];
		double tot = 0;
		int nstate = GetNstate(i);
		for (int k=0; k<nstate; k++)	{
			tot += (*t++);
		}
		if (tot == 0)	{
			// dirty !
			tot = 1e-12;
		}
		sitelogl = log(tot) + (*t);
		t -= nstate;
	}
	else	{
		double max = 0;
		double logl[GetNrate(i)];
		for (int j=0; j<GetNrate(i); j++)	{
			double* t = aux[j];
			double tot = 0;
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				tot += (*t++);
			}
			if (tot < 0)	{
				cerr << "negative total in site likelihood\n";
				exit(1);
				// dirty !
				tot = 1e-12;
			}
			logl[j] = log(tot) + (*t);
			t -= nstate;
			if ((!j) || (max < logl[j]))	{
				max = logl[j];
			}
		}
		double total = 0;
		double meanrate = 0;
		for (int j=0; j<GetNrate(i); j++)	{
			double tmp = 0;
			if (! isinf(logl[j]))	{
				tmp = GetRateWeight(i,j) * exp(logl[j] - max);
			}
			total += tmp;
			meanrate += tmp * GetRate(i,j);
		}
		sitelogl = log(total) + max;
		if (isnan(sitelogl))    {
			cerr << "nan site logl\n";
			cerr << max << '\t' << total << '\n';
			for (int j=0; j<GetNrate(i); j++)	{
				cerr << '\t' << logl[j] << '\t' << exp(logl[j] - max) << '\n';
			}
			exit(1);
		}
		meanrate /= total;
		// meansiterate[i] = meanrate;
	}
	return sitelogl;
}
	

//-------------------------------------------------------------------------
//	* sample the allocation of each site to one of the available rate categories
// 	for each site, the rate caegory is chosen with probability proportional
//	to the likelihood for that site, under that rate category
//	is aux == 0, this likelihood is assumed to be already computed and stored in int** condsitelogL (a member of SubstitutionProcess)
//	otherwise, aux is assumed to be a conditional likelihood vector multiplied by stationaries, and thus can be used to compute those likelihoods 
//	(CPU level 2)
//-------------------------------------------------------------------------
/*
void SubstitutionProcess::SiteDrawAllocations(int i, double** aux)	{

	if (aux)	{
		SiteComputeLikelihood(i,aux);
	}

	if (GetNrate(i) == 1)	{
		ratealloc[i] = 0;
	}
	else	{
		double max = 0;
		double* logl = condsitelogL[i];
		for (int j=0; j<GetNrate(i); j++)	{
			if ((!j) || (max < logl[j]))	{
				max = logl[j];
			}
		}
		double cumul = 0;
		double* p = new double[GetNrate(i)];
		for (int j=0; j<GetNrate(i); j++)	{
			cumul += GetRateWeight(i,j) * exp(logl[j] - max);
			p[j] = cumul;
		}
		double u = rnd::GetRandom().Uniform() * cumul;
		int j = 0;
		while ((j<GetNrate(i)) && (p[j] < u)) j++;
		if (j == GetNrate(i))	{
			cerr << "error in SubstitutionProcess::SampleAlloc\n";
			exit(1);
		}
		ratealloc[i] = j;
		delete[] p;
	}
}

//-------------------------------------------------------------------------
//	* sample states for each site, based on vectors of posterior probability stored in double*** t
//	this is assumed to be conditional on rate allocations 
//	(CPU level 2)
//-------------------------------------------------------------------------

int SubstitutionProcess::SiteChooseState(int i, double** t)	{

	int j = ratealloc[i];
	double* tmp = t[j];
	double total = 0;
	for (int k=0; k<GetNstate(i); k++)	{
		total += tmp[k];
	}
	double u = rnd::GetRandom().Uniform() * total;
	double tot = tmp[0];
	int k = 0;
	while ((k<GetNstate(i)) && (tot < u))	{
		k++;
		if (k==GetNstate(i))	{
			cerr << "error in SubstitutionProcess::ChooseState\n";
			exit(1);
		}
		tot += tmp[k];
	}
	for (int l=0; l<GetNstate(i); l++)	{
		tmp[l] = 0;
	}
	tmp[k] =  1;
	return k;
}
*/

