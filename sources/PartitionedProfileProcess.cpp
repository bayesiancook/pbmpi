
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionedProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionedProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PartitionedProfileProcess::Create(int indim, PartitionScheme inscheme)	{
	if (! profile)	{
		ProfileProcess::Create(inscheme.GetNsite(),indim);
		PartitionProcess::Create(inscheme);
		allocprofile = new double[GetNpart() * GetDim()];
		profile = new double*[GetNpart()];
		for (int i=0; i<GetNpart(); i++)	{
			profile[i] = allocprofile + i*GetDim();
			// profile[i] = new double[GetDim()];
		}
		dirweight = new double[GetDim()];
		logstatprior = new double[GetNpart()];
		profilesuffstatlogprob = new double[GetNpart()];
		fixstat = new bool[GetNpart()];

		SampleProfile();

		for(int p = 0; p < inscheme.Npart; p++)
			SetStat(p, inscheme.partType[p]);
	}
}

void PartitionedProfileProcess::Delete()	{
	if (profile)	{
		delete[] profilesuffstatlogprob;
		delete[] logstatprior;
		delete[] allocprofile;
		delete[] profile;
		delete[] fixstat;
		profile = 0;

		PartitionedProfileProcess::Delete();
		ProfileProcess::Delete();
	}
}

double PartitionedProfileProcess::GetMeanDirWeight()	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += dirweight[k];
	}
	return total;
}

double PartitionedProfileProcess::GetStatEnt()	{
	double total = 0;
	for (int k=0; k<GetNpart(); k++)	{
		total += GetPartNsite(k) * GetStatEnt(k);
	}
	return total / GetNsite();
}

double PartitionedProfileProcess::GetStatEnt(int k)	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		if (profile[k][i] <= 0)	{
			cerr << "error: 0 entry in profile\n";
			cerr << profile[k][i] << '\n';
			exit(1);
		}
		total -= profile[k][i] * log(profile[k][i]);
	}
	if (isnan(total))	{
		cerr << "entropy is nan\n";
		exit(1);
	}
	return  total;
}

double PartitionedProfileProcess::GetCenterStatEnt()	{
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

void PartitionedProfileProcess::RenormalizeProfiles()	{
	for (int i=0; i<GetNpart(); i++)	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += profile[i][k];
		}
		for (int k=0; k<GetDim(); k++)	{
			profile[i][k] /= total;
		}
	}
}

void PartitionedProfileProcess::SampleProfile()	{
	SampleHyper();
	SampleStat();
	// UpdateComponents();
}

void PartitionedProfileProcess::SampleStat()	{
	for (int i=0; i<GetNpart(); i++)	{
		SampleStat(profile[i]);
	}
}

void PartitionedProfileProcess::SampleStat(int i)	{
	SampleStat(profile[i]);
}

void PartitionedProfileProcess::SampleStat(double* prof, double statmin)	{
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

double PartitionedProfileProcess::LogProfilePrior()	{
	double total = 0;
	if(nfreestat > 1)
	{
		total += LogHyperPrior();
		total += LogStatPrior();
	}
	return total;
}

double PartitionedProfileProcess::LogStatPrior()	{

	double total = 0;
	for (int i=0; i<GetNpart(); i++)	{
		if(!fixstat[i])
			total += LogStatPrior(i);
	}
	return total;
}

double PartitionedProfileProcess::LogStatPrior(int cat)	{
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

double PartitionedProfileProcess::ProfileSuffStatLogProb()	{
	// simply, sum over all components
	for (int i=0; i<GetNpart(); i++)	{
		ProfileSuffStatLogProb(i);
	}
	double total = 0;
	for (int i=0; i<GetNpart(); i++)	{
		total += profilesuffstatlogprob[i];
	}
	return total;
}

double PartitionedProfileProcess::MoveDirWeights(double tuning, int nrep)	{
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
	return naccepted / nrep / GetDim();
}

void PartitionedProfileProcess::SampleHyper()	{
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = 1.0;
	}
}

double PartitionedProfileProcess::LogHyperPrior()	{
	double total = 0;
	double sum = 0;
	for (int k=0; k<GetDim(); k++)	{
		total -= dirweight[k];
		sum += dirweight[k];
	}
	if (sum < GetMinTotWeight())	{
		// total += InfProb;
		total -= 1.0 / 0;
	}
	return total;
}

double PartitionedProfileProcess::MoveHyper(double tuning, int nrep)	{
	double total = 0;
	total += MoveDirWeights(tuning,nrep);
	return total;
}

void PartitionedProfileProcess::SetStat(int inpart, string type)	{

	int p = inpart;
	if (type != "None")	{
		fixstat[inpart] = true;

		const double * Stat = 0;

		if (type == "dayhoff")	{
			Stat = DAYHOFF_Stat;
		}

		else if (type == "dcmut")	{
			Stat = DCMUT_Stat;
		}

		else if (type == "jtt")	{
			Stat = JTT_Stat;
		}

		else if (type == "mtrev")	{
			Stat = MTREV_Stat;
		}

		else if (type == "wag")	{
			Stat = WAG_Stat;
		}

		else if (type == "rtrev")	{
			Stat = RTREV_Stat;
		}

		else if (type == "cprev")	{
			Stat = CPREV_Stat;
		}

		else if (type == "vt")	{
			Stat = VT_Stat;
		}

		else if (type == "blosum62")	{
			Stat = BLOSUM62_Stat;
		}

		else if (type == "mtmam")	{
			Stat = MTMAM_Stat;
		}

		else if (type == "lg")	{
			Stat = LG_Stat;
		}

		else if (type == "mtart")	{
			Stat = MTART_Stat;
		}

		else if (type == "mtzoa")	{
			Stat = MTZOA_Stat;
		}

		else if (type == "pmb")	{
			Stat = PMB_Stat;
		}

		else if (type == "hivb")	{
			Stat = HIVB_Stat;
		}

		else if (type == "hivw")	{
			Stat = HIVW_Stat;
		}

		else if (type == "jttdcmut")	{
			Stat = JTTDCMUT_Stat;
		}

		else if (type == "flu")	{
			Stat = FLU_Stat;
		}

		if(Stat != 0){
			if (GetDim() != 20)	{
				if (! GetMyid())	{
					cerr << "error : empirical matrices only apply to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<GetDim(); i++)	{
				profile[p][i] = Stat[i];
			}
		}
		else	{
			ifstream is(type.c_str());
			if (!is)	{
				cerr << "unrecognized file for profiles : " << type << '\n';
				exit(1);
			}

			int permut[GetDim()];
			for (int k=0; k<GetDim(); k++)	{
				string c;
				is >> c;
				permut[k] = GetStateSpace()->GetState(c);
			}
			for(int pi = 0; pi < p; pi++)
			{
				for (int l=0; l<GetDim(); l++)	{
					double tmp;
					is >> tmp;
					if (tmp < 0)	{
						if (! GetMyid())	{
							cerr << "error when reading profile from " << type << '\n';
							cerr << tmp << '\n';
							cerr << '\n';
						}
						MPI_Finalize();
						exit(1);
					}
					profile[p][l] = tmp;

				}
			}
		}
	}
	else
	{
		fixstat[inpart] = false;
	}

	nfreestat = GetNpart();
	for(int p = 0; p < GetNpart(); p++)
	{
		nfreestat -= fixstat[p];
	}
}

