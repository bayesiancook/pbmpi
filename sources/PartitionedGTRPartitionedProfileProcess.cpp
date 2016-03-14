
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionedGTRPartitionedProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"
#include <algorithm>


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionedGTRPartitionedProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PartitionedGTRPartitionedProfileProcess::Create(int indim, PartitionScheme rrscheme, PartitionScheme statscheme)	{
	if (! matrixarray)	{
		PartitionedProfileProcess::Create(indim, statscheme);
		PartitionedGTRProfileProcess::Create(indim, rrscheme);
		matrixarray = new SubMatrix**[PartitionedGTRProfileProcess::GetNpart()];
		for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
			matrixarray[p] = new SubMatrix*[PartitionedProfileProcess::GetNpart()];
			for (int k=0; k<PartitionedProfileProcess::GetNpart(); k++)	{
				matrixarray[p][k] = 0;
			}
		}
		// SampleProfile();
	}

	partitionmap.clear();
	partitionmap.resize(rrscheme.Npart);
	for(size_t p = 0; p < rrscheme.Npart; p++)
	{
		std::vector<int> partsites = rrscheme.partSites[p];

		for(size_t i = 0; i < partsites.size(); i++)
		{
			int statpart = statscheme.sitePart[partsites[i]];
			if(std::find(partitionmap[p].begin(), partitionmap[p].end(), statpart) == partitionmap[p].end())
			{
				partitionmap[p].push_back(statpart);
			}
		}
	}
}

void PartitionedGTRPartitionedProfileProcess::Delete() {
	if (matrixarray)	{
		for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
			for (int k=0; k< partitionmap[p].size(); k++)	{
				delete matrixarray[p][partitionmap[p][k]];
			}
			delete [] matrixarray[p];
		}
		delete[] matrixarray;
		matrixarray = 0;

		PartitionedGTRProfileProcess::Delete();
		PartitionedProfileProcess::Delete();
	}
}

void PartitionedGTRPartitionedProfileProcess::CreateMatrix(int p, int k)	{
	if (matrixarray[p][k])	{
		cerr << "error in gen path suff stat gtr dp profile process: matrixarray is not 0\n";
		cerr << matrixarray[p][k]->GetNstate() << '\n';
		exit(1);
	}
	// matrixarray[k] = new LGSubMatrix(profile[k],false);
	matrixarray[p][k] = new GTRSubMatrix(GetDim(),rr[p],profile[k],false);
}

void PartitionedGTRPartitionedProfileProcess::UpdateMatrix(int p, int k)	{
	/*
	if (! matrixarray[k])	{
		cerr << "error in GeneralPathSuffStatPartitionedGTRMixtureProfileProcess::UpdateComponent: null matrix\n";
		exit(1);
	}
	*/
	if (matrixarray[p][k])	{
	GetGTRMatrix(p, k)->CorruptMatrix();
	}
	// GetGTRMatrix(k)->UpdateMatrix();
}

double PartitionedGTRPartitionedProfileProcess::GlobalMoveProfile(double tuning, int n, int nrep)	{

	assert(GetMyid() == 0);

	int Nocc = PartitionedProfileProcess::GetNpart();

	// send PROFILE_MOVE Message with n and nrep and tuning
	MESSAGE signal = PROFILE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int* itmp = new int[2];
	itmp[0] = n;
	itmp[1] = nrep;
	MPI_Bcast(itmp,2,MPI_INT,0,MPI_COMM_WORLD);
	delete[] itmp;

	double* dtmp = new double[1+Nocc*GetDim()];
	dtmp[0] = tuning;
	int k = 1;
	for (int i=0; i< Nocc; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			dtmp[k] = profile[i][j];
			k++;
		}
	}
	MPI_Bcast(dtmp,1+Nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	delete[] dtmp;

	// split Ncomponent items among GetNprocs() - 1 slaves
	int width = Nocc/(GetNprocs()-1);
	int maxwidth = 0;
	int cmin[GetNprocs()-1];
	int cmax[GetNprocs()-1];
	int dmin[GetNprocs()-1];
	int dmax[GetNprocs()-1];

	for(int i=0; i<GetNprocs()-1; ++i) {

		int ddmin = width * i;
		int ddmax = (i == GetNprocs() - 2) ? Nocc : width * (i+1);
		if (maxwidth < (ddmax - ddmin))	{
			maxwidth = ddmax - ddmin;
		}

		int k = -1;
		int ccmin = -1;
		while ((ccmin< Nocc) && (k<ddmin))	{
			ccmin++;
			if (ccmin == Nocc)	{
				cerr << "error in matmixslavemoveprofile: overflow\n";
				exit(1);
			}
			k++;
		}
		int ccmax = ccmin;
		if (ddmax == Nocc)	{
			ccmax = Nocc;
		}
		else	{
			while ((ccmax< Nocc) && (k<ddmax))	{
				ccmax++;
				k++;
			}
		}

		cmin[i] = ccmin;
		cmax[i] = ccmax;
		dmin[i] = ddmin;
		dmax[i] = ddmax;
		if ((ccmax - ccmin) != (dmax[i] - dmin[i]))	{
			cerr << "error: non matching numbers: " << i << '\t' << ccmin << '\t' << ccmax << '\t' << ddmin << '\t' << ddmax << '\n';
			cerr << '\n';
			exit(1);
		}
	}

	// collect final values of profiles (+ total acceptance rate) from slaves
	MPI_Status stat;
	int bigdim = maxwidth * GetDim();
	double* tmp = new double[bigdim+1]; // (+1 for the acceptance rate)
	double total = 0;
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(tmp,(dmax[i-1]-dmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		int l = 0;
		for(int j=cmin[i-1]; j<cmax[i-1]; ++j) {
			for (int k=0; k<GetDim(); k++)	{
				profile[j][k] = tmp[l];
				l++;
			}
		}
		total += tmp[l]; // (sum all acceptance rates)
	}
	delete[] tmp;
	// return average acceptance rate
	return total / Nocc;
}

void PartitionedGTRPartitionedProfileProcess::SlaveMoveProfile()	{

	assert(GetMyid() > 0);

	// parse arguments sent by master

	int* itmp = new int[2];
	MPI_Bcast(itmp,2,MPI_INT,0,MPI_COMM_WORLD);
	int n = itmp[0];
	int nrep = itmp[1];
	delete[] itmp;

	int Nocc = PartitionedProfileProcess::GetNpart();

	double* dtmp = new double[1 + Nocc*GetDim()];
	MPI_Bcast(dtmp,1+Nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	double tuning = dtmp[0];
	int k = 1;
	for (int i=0; i<Nocc; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			profile[i][j] = dtmp[k];
			k++;
		}
	}
	delete[] dtmp;

	// determine the range of components to move
	int width = Nocc/(GetNprocs()-1);
	int dmin = width * (GetMyid() - 1);
	int dmax = (GetMyid() == GetNprocs() - 1) ? Nocc : width * GetMyid();

	k = -1;
	int cmin = -1;
	while ((cmin<Nocc) && (k<dmin))	{
		cmin++;
		if (cmin == Nocc)	{
			cerr << "error in matmixslavemoveprofile: overflow\n";
			exit(1);
		}
		k++;
	}
	int cmax = cmin;
	if (dmax == Nocc)	{
		cmax = Nocc;
	}
	else	{
		while ((cmax<Nocc) && (k<dmax))	{
			cmax++;
			k++;
		}
	}
	int nocc = 0;
	for (int j=cmin; j<cmax; j++)	{
		nocc++;
	}
	if (nocc != (dmax - dmin))	{
		cerr << "error : mismatch in nocc\n";
		exit(1);
	}

	// update sufficient statistics
	UpdateModeProfileSuffStat();
	// move components in the range just computed
	double total = 0;
	for (int i=cmin; i<cmax; i++)	{
		if(!fixstat[i])
			total += MoveProfile(i,tuning,n,nrep);
	}

	// send the new values of the profiles, plus the total success rate (total)
	double* tmp = new double[(dmax - dmin) * GetDim() + 1];
	int l = 0;
	for (int i=cmin; i<cmax; i++)	{
		for (int k=0; k<GetDim(); k++)	{
			tmp[l] = profile[i][k];
			l++;
		}
	}
	tmp[l] = total;
	MPI_Send(tmp,(dmax-dmin)*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	delete[] tmp;
}

double PartitionedGTRPartitionedProfileProcess::MoveProfile(double tuning, int n, int nrep)	{
	double total = 0;
	for (int i=0; i<PartitionedProfileProcess::GetNpart(); i++)	{
		if(!fixstat[i])
			total += MoveProfile(i,tuning,n,nrep);
	}
	return total / PartitionedProfileProcess::GetNpart();
}

double PartitionedGTRPartitionedProfileProcess::MoveProfile(int cat, double tuning, int n, int nrep)	{
// double PartitionedGTRPartitionedProfileProcess::MoveProfile(int cat, double tuning, int n, int nrep, RngStream& g)	{
	int naccepted = 0;
	double* bk = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bk[k] = profile[cat][k];
	} 
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogStatPrior(cat) - ProfileSuffStatLogProb(cat);
		double loghastings = ProfileProposeMove(profile[cat],tuning,n,0);
		// double loghastings = ProfileProposeMove(profile[cat],tuning,n,0,g);
		UpdateComponent(cat);
		deltalogprob += LogStatPrior(cat) + ProfileSuffStatLogProb(cat);
		deltalogprob += loghastings;
		// int accepted = (g.RandU01() < exp(deltalogprob));
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted ++;
			for (int k=0; k<GetDim(); k++)	{
				bk[k] = profile[cat][k];
			} 
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[cat][k] = bk[k];
			} 
			UpdateComponent(cat);
		}
	}
	delete[] bk;
	return naccepted / nrep;
}

