
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonSBDPProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"

double PoissonSBDPProfileProcess::GlobalMixMove(int nrep, int nallocrep, double epsilon)	{

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
	int itmp[3];
	itmp[0] = nrep;
	itmp[1] = nallocrep;
	itmp[2] = K0;
	MPI_Bcast(itmp,3,MPI_INT,0,MPI_COMM_WORLD);

	// split Nsite among GetNprocs()-1 slaves
	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
	}

	/*
	ResampleEmptyProfiles();
	MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/

	double* tmp = new double[Ncomponent * GetDim() + 1];

	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// here slaves do realloc moves

		// mpi receive new allocations
		MPI_Status stat;
		int tmpalloc[GetNsite()+1];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=smin[i-1]; j<smax[i-1]; ++j) {
				alloc[j] = tmpalloc[j];
				if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
					cerr << "alloc overflow\n";
					exit(1);
				}
			}
		}

		// MPI_Barrier(MPI_COMM_WORLD);

		// broadcast new allocations
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// here slaves do profile moves

		UpdateOccupancyNumbers();

		// collect final values of profiles (+ total acceptance rate) from slaves

		// split Nmode among GetNprocs()-1 slaves
		int mwidth = GetNcomponent()/(GetNprocs()-1);
		int mmin[GetNprocs()-1];
		int mmax[GetNprocs()-1];
		for(int i=0; i<GetNprocs()-1; ++i) {
			mmin[i] = mwidth*i;
			mmax[i] = mwidth*(1+i);
			if (i == (GetNprocs()-2)) mmax[i] = GetNcomponent();
		}

		MPI_Status stat2;
		double total = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,(mmax[i-1]-mmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat2);
			int l = 0;
			for(int j=mmin[i-1]; j<mmax[i-1]; ++j) {
				for (int k=0; k<GetDim(); k++)	{
					profile[j][k] = tmp[l];
					l++;
				}
			}
			total += tmp[l]; // (sum all acceptance rates)
		}

		// MPI_Barrier(MPI_COMM_WORLD);
		// resend all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	delete[] tmp;

	return ((double) NAccepted) / GetNsite() / nrep;
}

void PoissonSBDPProfileProcess::SlaveMixMove()	{

	int itmp[3];
	MPI_Bcast(itmp,3,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
	int K0 = itmp[2];

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];
	double* tmp = new double[Ncomponent * GetDim() + 1];

	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
	}

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

			for (int site=smin[GetMyid()-1]; site<smax[GetMyid()-1]; site++)	{
			// for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

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
				if (std::isnan(total))	{
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
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		// MPI_Barrier(MPI_COMM_WORLD);
		// profile move

		// receive new allocations
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// determine the range of components to move
		UpdateOccupancyNumbers();

		// update sufficient statistics
		UpdateModeProfileSuffStat();

		// split Nmode among GetNprocs()-1 slaves
		int mwidth = GetNcomponent()/(GetNprocs()-1);
		int mmin[GetNprocs()-1];
		int mmax[GetNprocs()-1];
		for(int i=0; i<GetNprocs()-1; ++i) {
			mmin[i] = mwidth*i;
			mmax[i] = mwidth*(1+i);
			if (i == (GetNprocs()-2)) mmax[i] = GetNcomponent();
		}

		double total = 0;
		for (int mode=mmin[GetMyid()-1]; mode<mmax[GetMyid()-1]; mode++)	{
			total += MoveProfile(mode);
		}
		int l = 0;
		for (int mode=mmin[GetMyid()-1]; mode<mmax[GetMyid()-1]; mode++)	{
			for (int k=0; k<GetDim(); k++)	{
				tmp[l] = profile[mode][k];
				l++;
			}
		}
		tmp[l] = total;

		MPI_Send(tmp,(mmax[GetMyid()-1] - mmin[GetMyid()-1])*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);

		// rereceive all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	delete[] cumul;
	delete[] mLogSamplingArray;
	delete[] tmp;
}
