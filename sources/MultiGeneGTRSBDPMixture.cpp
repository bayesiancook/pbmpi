
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <cassert>
#include "Parallel.h"
#include <string.h>

#include "MultiGeneGTRSBDPMixture.h"

MultiGeneGTRSBDPMixture::MultiGeneGTRSBDPMixture(string indatafile, string treefile, string inname, int innratecat, int infixtopo, int inkappaprior, double ingibbsfactor, int indc, int me, int np)	{
	myid = me;
	nprocs = np;
	name = inname;
	fixtopo = infixtopo;
	dc = indc;
	kappaprior = inkappaprior;
	Nstate = 0;
	datafile = indatafile;
	nratecat = innratecat;
	gibbsfactor = ingibbsfactor;
	AllocateAlignments(datafile,treefile,dc);
	
	Create(GetGlobalNsite(0),Nstate,nratecat,indc);
	if (! myid)	{
		GlobalUpdateParameters();
		GlobalSample();
		GlobalUnfold();
	}

	MakeFiles();
}

MultiGeneGTRSBDPMixture::MultiGeneGTRSBDPMixture(istream& is, int me, int np)	{

	myid = me;
	nprocs = np;
	FromStreamHeader(is);

	Nstate = 0;
	AllocateAlignments(datafile,"None",dc);
	
	Create(GetGlobalNsite(0),Nstate,nratecat,dc);
	if (! myid)	{
		FromStream(is);
		GlobalUpdateParameters();
		GlobalUnfold();
	}
}

void MultiGeneGTRSBDPMixture::Create(int inNsite, int Nstate, int nratecat, int dc)	{

	// create Nsite Nstate
	ExpoConjugateGTRSBDPProfileProcess::Create(inNsite,Nstate);
	if (!myid)	{
		SampleProfile();
		// CreateMatrices();
		// ???
	}
	else	{
		// CreateMatrices();
		process = new GenePhyloProcess*[Ngene];
		// process = new ExpoConjugateGTRGeneGammaPhyloProcess*[Ngene];
		int offset = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				process[gene] = new ExpoConjugateGTRGeneGammaPhyloProcess(genename[gene], treename[gene], nratecat, fixtopo, dc, this, offset);
				process[gene]->SetGibbsFactor(gibbsfactor);
			}
			else	{
				process[gene] = 0;
			}
			offset += genesize[gene];
		}
	}
	CreateSuffStat();
}

void MultiGeneGTRSBDPMixture::CreateSuffStat()	{

	allocsiteprofilesuffstatcount = new int[GetGlobalNsite(0) * GetDim()];
	allocsiteprofilesuffstatbeta = new double[GetGlobalNsite(0) * GetDim()];

	siteprofilesuffstatcount = new int*[GetGlobalNsite(0)];
	siteprofilesuffstatbeta  = new double*[GetGlobalNsite(0)];
	for (int i=0; i<GetGlobalNsite(0); i++)	{
		siteprofilesuffstatcount[i] = allocsiteprofilesuffstatcount + i*GetDim();
		siteprofilesuffstatbeta[i] = allocsiteprofilesuffstatbeta + i*GetDim();
	}

	if (! myid)	{
		tmpallocsiteprofilesuffstatcount = new int[GetGlobalNsite(0) * GetDim()];
		tmpallocsiteprofilesuffstatbeta = new double[GetGlobalNsite(0) * GetDim()];

		tmpsiteprofilesuffstatcount = new int*[GetGlobalNsite(0)];
		tmpsiteprofilesuffstatbeta  = new double*[GetGlobalNsite(0)];
		for (int i=0; i<GetGlobalNsite(0); i++)	{
			tmpsiteprofilesuffstatcount[i] = tmpallocsiteprofilesuffstatcount + i*GetDim();
			tmpsiteprofilesuffstatbeta[i] = tmpallocsiteprofilesuffstatbeta + i*GetDim();
		}

		tmprrsuffstatcount = new int[GetNrr()];
		tmprrsuffstatbeta = new double[GetNrr()];
	}
}

void MultiGeneGTRSBDPMixture::DeleteSuffStat()	{
	
	delete[] allocsiteprofilesuffstatcount;
	delete[] allocsiteprofilesuffstatbeta;

	for (int i=0; i<GetGlobalNsite(0); i++)	{
		delete[] siteprofilesuffstatcount[i];
		delete[] siteprofilesuffstatbeta[i];
	}
	delete[] siteprofilesuffstatcount;
	delete[] siteprofilesuffstatbeta;

	if (! myid)	{
		delete[] tmprrsuffstatcount;
		delete[] tmprrsuffstatbeta;
	}
}

void MultiGeneGTRSBDPMixture::GlobalUpdateParameters()	{

	assert(myid == 0);
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int nd = GetNrr() + GetNmodeMax() * (GetDim() + 1) + GetDim() + 1;
	double dvector[nd]; 
	int k = 0;
	for(int i=0; i<GetNrr(); ++i) {
		dvector[k] = rr[i];
		k++;
	}

	for(int i=0; i<GetNmodeMax(); ++i) {
		for(int j=0; j<GetDim(); ++j) {
			dvector[k] = profile[i][j];
			k++;
		}
		dvector[k] = weight[i];
		k++;
	}
	dvector[k] = kappa;
	k++;
	for (int i=0; i<GetDim(); i++)	{
		dvector[k] = dirweight[i];
		k++;
	}
	if (k != nd)	{
		cerr << "error in EMGTRMixture global update mix params\n";
		exit(1);
	}

	MPI_Bcast(&Ncomponent,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(alloc,GetGlobalNsite(0),MPI_INT,0,MPI_COMM_WORLD);
}


void MultiGeneGTRSBDPMixture::SlaveUpdateParameters()	{
	int nd = GetNrr() + GetNmodeMax() * (GetDim() + 1) + GetDim() + 1;
	double* dvector = new double[nd];
	MPI_Bcast(&Ncomponent,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(alloc,GetGlobalNsite(0),MPI_INT,0,MPI_COMM_WORLD);

	int k = 0;
	for(int i=0; i<GetNrr(); ++i) {
		rr[i] = dvector[k];
		k++;
	}

	for(int i=0; i<GetNmodeMax(); ++i) {
		for(int j=0; j<GetDim(); ++j) {
			profile[i][j] = dvector[k];
			k++;
		}
		weight[i] = dvector[k];
		k++;
	}
	kappa = dvector[k];
	k++;
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[k];
		k++;
	}
	if (k != nd)	{
		cerr << "error in EMGTRMixture slave update mix params\n";
		exit(1);
	}

	delete[] dvector;
	CreateMatrices();
	UpdateMatrices();

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SetMixtureParameters();
		}
	}
}

void MultiGeneGTRSBDPMixture::GlobalUpdateSiteProfileSuffStat()	{

	// send signal
	assert(myid == 0);
	MESSAGE signal = UPDATE_SPROFILE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// collect from slaves
	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpallocsiteprofilesuffstatcount, GetGlobalNsite(0) * GetDim(),MPI_INT,j,TAG1,MPI_COMM_WORLD,&stat);
		MPI_Recv(tmpallocsiteprofilesuffstatbeta, GetGlobalNsite(0) * GetDim(),MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		int offset = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				for (int i=0; i<genesize[gene]; i++)	{
					for (int k=0; k<GetDim(); k++)	{
						siteprofilesuffstatcount[i+offset][k] = tmpsiteprofilesuffstatcount[i+offset][k];
						siteprofilesuffstatbeta[i+offset][k] = tmpsiteprofilesuffstatbeta[i+offset][k];
					}
				}
			}
			offset += genesize[gene];
		}
	}

	// MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(allocsiteprofilesuffstatcount,GetGlobalNsite(0)*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocsiteprofilesuffstatbeta,GetGlobalNsite(0)*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	/*
	int index = 0;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(allocsiteprofilesuffstatcount + index * GetDim(), GetGlobalNsite(j) * GetDim(),MPI_INT,j,TAG1,MPI_COMM_WORLD,&stat);
		MPI_Recv(allocsiteprofilesuffstatbeta + index * GetDim(), GetGlobalNsite(j) * GetDim(),MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		index += GetGlobalNsite(j);
	}
	if (index != GetGlobalNsite())	{
		cerr << "error : non matching length\n";
		cerr << index << '\t' << GetGlobalNsite() << '\n';
		exit(1);
	}
	*/
}

void MultiGeneGTRSBDPMixture::SlaveUpdateSiteProfileSuffStat()	{

	// compute gene specific site profile suff stats
	int offset = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ExpoConjugateGTRGeneGammaPhyloProcess* localprocess = GetGTRProcess(gene);
			localprocess->UpdateSiteProfileSuffStat();
			for (int i=0; i<genesize[gene]; i++)	{
				for (int k=0; k<GetDim(); k++)	{
					siteprofilesuffstatcount[i+offset][k] = localprocess->siteprofilesuffstatcount[i][k];
					siteprofilesuffstatbeta[i+offset][k] = localprocess->siteprofilesuffstatbeta[i][k];
				}
			}
		}
		offset += genesize[gene];
	}

	MPI_Send(allocsiteprofilesuffstatcount, GetGlobalNsite(0) * GetDim(), MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(allocsiteprofilesuffstatbeta, GetGlobalNsite(0) * GetDim(), MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	// MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(allocsiteprofilesuffstatcount,GetGlobalNsite(0)*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocsiteprofilesuffstatbeta,GetGlobalNsite(0)*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	/*
	offset = 0;
	for (int j=1; j<myid; j++)	{
		offset += GetGlobalNsite(j);
	}

	// send back to master
	MPI_Send(allocsiteprofilesuffstatcount + offset * GetDim(), GetGlobalNsite() * GetDim(), MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(allocsiteprofilesuffstatbeta + offset * GetDim(), GetGlobalNsite() * GetDim(), MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	*/
}

void MultiGeneGTRSBDPMixture::GlobalUpdateRRSuffStat()	{

	// send signal
	assert(myid == 0);
	MESSAGE signal = UPDATE_RRATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// collect from slaves
	if (! rrsuffstatcount)	{
		cerr << "rr suff stat arrays not created\n";
		exit(1);
	}
	for (int k=0; k<GetNrr(); k++)	{
		rrsuffstatcount[k] = 0;
		rrsuffstatbeta[k] = 0;
	}
	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmprrsuffstatcount, GetNrr(),MPI_INT,j,TAG1,MPI_COMM_WORLD,&stat);
		MPI_Recv(tmprrsuffstatbeta, GetNrr(),MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int k=0; k<GetNrr(); k++)	{
			rrsuffstatcount[k] += tmprrsuffstatcount[k];
			rrsuffstatbeta[k] += tmprrsuffstatbeta[k];
		}
	}
}

void MultiGeneGTRSBDPMixture::SlaveUpdateRRSuffStat()	{

	// compute gene specific site profile suff stats
	if (! rrsuffstatcount)	{
		cerr << "rr suff stat arrays not created\n";
		exit(1);
	}
	for (int k=0; k<GetNrr(); k++)	{
		rrsuffstatcount[k] = 0;
		rrsuffstatbeta[k] = 0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ExpoConjugateGTRGeneGammaPhyloProcess* localprocess = GetGTRProcess(gene);
			localprocess->UpdateRRSuffStat();
			for (int k=0; k<GetNrr(); k++)	{
				rrsuffstatcount[k] += localprocess->rrsuffstatcount[k];
				rrsuffstatbeta[k] += localprocess->rrsuffstatbeta[k];
			}
		}
	}

	// send back to master
	MPI_Send(rrsuffstatcount, GetNrr(), MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(rrsuffstatbeta, GetNrr(), MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneGTRSBDPMixture::LengthRelRateMove(double tuning, int nrep)	{

	assert(myid == 0);

	double naccept = 0;

	for (int rep=0; rep<nrep; rep++)	{
		double deltalogratio = - LogRRPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);

		for (int i=0; i<GetNrr(); i++)	{
			rr[i] /= e;
		}

		deltalogratio += LogRRPrior();
		deltalogratio -= GetNrr() * m;

		MESSAGE signal = LENGTHFACTOR;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&e,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Status stat;
		double tmp = 0;
		for (int j=1; j<nprocs; j++)	{
			MPI_Recv(&tmp, 1, MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
			deltalogratio += tmp;
		}

		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

		if (accepted)	{
			naccept++;
		}
		else	{
			for (int i=0; i<GetNrr(); i++)	{
				rr[i] *= e;
			}
			double inve = 1.0 / e;
			MESSAGE signal = LENGTHFACTOR;
			MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&inve,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Status stat;
			double tmp = 0;
			for (int j=1; j<nprocs; j++)	{
				MPI_Recv(&tmp, 1, MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
			}
		}	
	}
	return naccept / nrep;
}

void MultiGeneGTRSBDPMixture::SlaveLengthFactorMove()	{
	// get factor
	double factor = 0;
	MPI_Bcast(&factor,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	double logratio = 0;
	int degree = 0; 
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			logratio -= process[gene]->LogLengthPrior();
			degree += process[gene]->MoveAllBranches(factor);
			logratio += process[gene]->LogLengthPrior();
		}
	}
	logratio += degree * log(factor);
	MPI_Send(&logratio,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneGTRSBDPMixture::SlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case LENGTHFACTOR:
		SlaveLengthFactorMove();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	
	default:
		MultiGeneMixture::SlaveExecute(signal);
	}
}


/*
double MultiGeneGTRSBDPMixture::GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep)	{

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}

	int NAccepted = 0;

	// define threshold between GIbbs and MH
	int K0 = GetNmodeMax();
	if (epsilon)	{
		double r = kappa / (1 + kappa);
		K0 = (int) (log(epsilon) / log(r));
		if (K0 >= GetNmodeMax())	{
			K0 = GetNmodeMax();
		}
	}
	// K0 = GetNmodeMax();

	// send mixmove signal and tuning parameters
	MESSAGE signal = MIX_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int itmp[4];
	itmp[0] = nrep;
	itmp[1] = nallocrep;
	itmp[2] = K0;
	itmp[3] = nprofilerep;
	MPI_Bcast(itmp,4,MPI_INT,0,MPI_COMM_WORLD);

	double* tmp = new double[Ncomponent * GetDim() + 1];

	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		// mpi send message for realloc move
		// mpi send profiles and weights
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);
		// MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

		// here slaves do realloc moves

		// mpi receive new allocations
		MPI_Status stat;
		int tmpalloc[GetGlobalNsite()+1];
		for(int j=1; j<GetNprocs(); ++j) {
			MPI_Recv(tmpalloc,GetGlobalNsite(),MPI_INT,j,TAG1,MPI_COMM_WORLD,&stat);
			int offset = 0;
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == j)	{
					for (int i=0; i<genesize[gene]; i++)	{
						alloc[i+offset] = tmpalloc[i+offset];
					}
				}
				offset += genesize[gene];
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// broadcast new allocations
		MPI_Bcast(alloc,GetGlobalNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// here slaves do profile moves

		// split Ncomponent items among GetNprocs() - 1 slaves
		UpdateOccupancyNumbers();
		int Nocc = GetNOccupiedComponent();
		int width = Nocc/(GetNprocs()-1);
		int cmin[GetNprocs()-1];
		int cmax[GetNprocs()-1];
		int dmin[GetNprocs()-1];
		int dmax[GetNprocs()-1];

		for(int i=0; i<GetNprocs()-1; ++i) {
			int ddmin = width * i;
			int ddmax = (i == GetNprocs() - 2) ? Nocc : width * (i+1);

			int k = -1;
			int ccmin = -1;
			while ((ccmin<Ncomponent) && (k<ddmin))	{
				ccmin++;
				if (ccmin == Ncomponent)	{
					cerr << "error in matmixslavemoveprofile: overflow\n";
					exit(1);
				}
				if (occupancy[ccmin])	{
					k++;
				}
			}
			int ccmax = ccmin;
			if (ddmax == Nocc)	{
				ccmax = Ncomponent;
			}
			else	{
				while ((ccmax<Ncomponent) && (k<ddmax))	{
					ccmax++;
					if (occupancy[ccmax])	{
						k++;
					}
				}
			}

			cmin[i] = ccmin;
			cmax[i] = ccmax;
			dmin[i] = ddmin;
			dmax[i] = ddmax;
			int nocc = 0;
			for (int j=ccmin; j<ccmax; j++)	{
				if (occupancy[j])	{
					nocc++;
				}
			}
			if (nocc != (dmax[i] - dmin[i]))	{
				cerr << "error: non matching numbers: " << i << '\t' << nocc << '\t' << ccmin << '\t' << ccmax << '\t' << ddmin << '\t' << ddmax << '\n';
				for (int j=ccmin; j<ccmax; j++)	{
					cerr << occupancy[j] << '\t';
				}
				cerr << '\n';
				exit(1);
			}
		}

		// collect final values of profiles (+ total acceptance rate) from slaves
		MPI_Status stat2;
		double total = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,(dmax[i-1]-dmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat2);
			int l = 0;
			for(int j=cmin[i-1]; j<cmax[i-1]; ++j) {
				if (occupancy[j])	{
					double tot = 0;
					for (int k=0; k<GetDim(); k++)	{
						profile[j][k] = tmp[l];
						tot += profile[j][k];
						l++;
					}
					if (fabs(tot - 1) > 1e-6)	{
						cerr << "normalization error : " << tot -1 << '\n';
						cerr << "upon receiving\n";
						exit(1);
					}
				}
			}
			total += tmp[l]; // (sum all acceptance rates)
		}

		// resample empty profiles
		ResampleEmptyProfiles();

		MPI_Barrier(MPI_COMM_WORLD);
		// resend all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	// check that profiles are normalized
	for (int k=0; k<Ncomponent; k++)	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			tot += profile[k][i];
		}
		if (fabs(tot - 1) > 1e-6)	{
			cerr << "normalization error : " << tot -1 << '\n';
			exit(1);
		}
	}
	delete[] tmp;

	// CreateMatrices();
	UpdateMatrices();

	return ((double) NAccepted) / GetGlobalNsite() / nrep;
}

void MultiGeneGTRSBDPMixture::SlaveMixMove()	{

	int itmp[4];
	MPI_Bcast(itmp,4,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
	int K0 = itmp[2];
	int nprofilerep = itmp[3];

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];
	double* tmp = new double[Ncomponent * GetDim() + 1];

	for (int rep=0; rep<nrep; rep++)	{

		// realloc move

		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);
		// MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		double totq = 0;
		for (int mode=K0; mode<GetNmodeMax(); mode++)	{
			totq += weight[mode];
		}

		int NAccepted = 0;

		for (int allocrep=0; allocrep<nallocrep; allocrep++)	{

			int offset = 0;
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == myid)	{
					for (int i=0; i<genesize[gene]; i++)	{
						int site = i+offset;

						int bk = alloc[site];

						double max = 0;
						// double mean = 0;
						for (int mode = 0; mode<K0; mode++)	{
							mLogSamplingArray[mode] = LogStatProb(site,mode);
							if ((!mode) || (max < mLogSamplingArray[mode]))	{
								max = mLogSamplingArray[mode];
							}
							// mean += mLogSamplingArray[mode];
						}
						// mean /= K0;

						double total = 0;
						for (int mode = 0; mode<K0; mode++)	{
							double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
							total += p;
							cumul[mode] = total;
						}
						if (isnan(total))	{
							cerr << "nan\n";
						}

						// double M = exp(mean- max);
						double M = 1;
						total += M * totq;
						double q = total * rnd::GetRandom().Uniform();
						int mode = 0;
						while ( (mode<K0) && (q > cumul[mode])) mode++;
						if (mode == K0)	{
							mode--;
							double r = (q - cumul[mode]) / M;
							while (r > 0)	{
								mode++;
								r -= weight[mode];
							}
						}

						// MH 
						double logratio = 0;
						if (mode >= K0)	{
							logratio += LogStatProb(site,mode) - max - log(M);
						}
						if (bk >= K0)	{
							logratio -= LogStatProb(site,bk) - max - log(M);
						}
						
						if (log(rnd::GetRandom().Uniform()) > logratio)	{
							mode = bk;
						}

						int Accepted = (mode != bk);
						if (Accepted)	{
							NAccepted ++;
						}
						alloc[site] = mode;

					}
				}
				offset += genesize[gene];
			}

		}

		MPI_Send(alloc,GetGlobalNsite(0),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		// profile move

		// receive new allocations
		MPI_Bcast(alloc,GetGlobalNsite(0),MPI_INT,0,MPI_COMM_WORLD);

		// determine the range of components to move
		UpdateOccupancyNumbers();
		int Nocc = GetNOccupiedComponent();
		int width = Nocc/(GetNprocs()-1);
		int dmin = width * (GetMyid() - 1);
		int dmax = (GetMyid() == GetNprocs() - 1) ? Nocc : width * GetMyid();

		int k = -1;
		int cmin = -1;
		while ((cmin<Ncomponent) && (k<dmin))	{
			cmin++;
			if (cmin == Ncomponent)	{
				cerr << "error in matmixslavemoveprofile: overflow\n";
				exit(1);
			}
			if (occupancy[cmin])	{
				k++;
			}
		}
		int cmax = cmin;
		if (dmax == Nocc)	{
			cmax = Ncomponent;
		}
		else	{
			while ((cmax<Ncomponent) && (k<dmax))	{
				cmax++;
				if (occupancy[cmax])	{
					k++;
				}
			}
		}
		int nocc = 0;
		for (int j=cmin; j<cmax; j++)	{
			if (occupancy[j])	{
				nocc++;
			}
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
			if (occupancy[i])	{
				total += MoveProfile(i,1,1,nprofilerep);
				total += MoveProfile(i,1,3,nprofilerep);
				total += MoveProfile(i,0.1,3,nprofilerep);
			}
		}


		// send the new values of the profiles, plus the total success rate (total)
		int l = 0;
		for (int i=cmin; i<cmax; i++)	{
			if (occupancy[i])	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tmp[l] = profile[i][k];
					tot += tmp[l];
					l++;
				}
				if (fabs(tot - 1) > 1e-6)	{
					cerr << "normalization error : " << tot -1 << '\n';
					cerr << "upon sending\n";
					cerr << cmin << '\t' << cmax << '\n';
					exit(1);
				}
			}
		}
		tmp[l] = total;
		
		MPI_Send(tmp,(dmax-dmin)*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// rereceive all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		UpdateMatrices();
	}

	delete[] cumul;
	delete[] mLogSamplingArray;
	delete[] tmp;
}

*/
