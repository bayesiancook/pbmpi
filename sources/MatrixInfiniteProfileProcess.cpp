
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixFiniteProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"


void MatrixFiniteProfileProcess::ResetNcomponent()	{

	int kmax = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (kmax < alloc[i])	{
			kmax = alloc[i];
		}
	}
	for (int k=kmax+1; k<Ncomponent; k++)	{
		DeleteComponent(k);
	}
	Ncomponent = kmax+1;

}

// K0: determines the transition between posterior and prior Gibbs
double MatrixFiniteProfileProcess::GlobalIncrementalFiniteMove(int nrep, int K0)	{

	assert(GetMyid() == 0);

	// send command and arguments
	MESSAGE signal = REALLOC_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&K0,1,MPI_INT,0,MPI_COMM_WORLD);

	// split Nsite among GetNprocs()-1 slaves
	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
	}

	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{

		// Ncomponent now equal to last occupied component 
		ResetNcomponent();

		// using the Vs and ps
		ResampleWeights();
		// here should have a cumulprod;

		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// MPI loop here:
		int nreceived = 0;
		while (nreceived < GetNprocs() - 1)	{
			// receive message
			MESSAGE signal;
			MPI_Status stat;

			int sender;
			MPI_Recv(&signal,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			MPI_Recv(&sender,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);

			if (signal == GIVEMEMORE)	{

				int K;
				MPI_Recv(&K,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
				double r;
				MPI_Recv(&r,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);

				int mode = K-1;
				while (r > 0)	{
					mode++;
					if (mode >= Ncomponent)	{
						if (mode >= GetNmodeMax())	{
							cerr << "nmode max overflow\n";
							exit(1);
						}
						CreateComponent(Ncomponent);
						// create weight;
						Ncomponent++;
					}
					r -= weight[mode];
				}

				// send
				MPI_Send(mode,1,MPI_INT,sender,TAG1,MPI_COMM_WORLD);
				MPI_Send(Ncomponent,1,MPI_INT,sender,TAG1,MPI_COMM_WORLD);
				MPI_Send(weight,Ncomponent,MPI_DOUBLE,sender,TAG1,MPI_COMM_WORLD);
			}
			// if message is slave has finished
			// fillup the array where it has finished
			else if (signal == REALLOC_DONE)	{
				int tmpalloc[GetNsite()];
				MPI_Recv(tmpalloc,GetNsite(),MPI_INT,sender,TAG1,MPI_COMM_WORLD,&stat);
				for(int j=smin[i-1]; j<smax[i-1]; ++j) {
					alloc[j] = tmpalloc[j];
					if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
						cerr << "alloc overflow\n";
						exit(1);
					}
				}
				nreceived++;
			}
			else	{
				cerr << "error in master realloc\n";
				exit(1);
			}
		}
	}
	
	// final cleanup
	UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();

	// CHECK that: might be useful depending on the exact submodel
	// ResampleWeights();
	// UpdateMatrices();

	return ((double) NAccepted) / GetNsite() / nrep;
}

double MatrixFiniteProfileProcess::SlaveIncrementalFiniteMove()	{

	// parse argument sent by master
	int nrep;
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);
	int K0;
	MPI_Bcast(&K0,1,MPI_INT,0,MPI_COMM_WORLD);

	int NAccepted = 0;

	vector<double> mLogSamplingArray;
	vector<double> cumul;

	for (int rep=0; rep<nrep; rep++)	{

		// receive weights sent by master
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		// do the incremental reallocation move on my site range
		for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

			int bk = alloc[site];

			double max = 0;
			for (int mode = 0; mode<K0; mode++)	{
				mLogSamplingArray[mode] =  LogStatProb(site,mode);
				if ((!mode) || (max < mLogSamplingArray[mode]))	{
					max = mLogSamplingArray[mode];
				}
			}

			double total = 0;
			for (int mode = 0; mode<K0; mode++)	{
				double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
				total += p;
				cumul[mode] = total;
			}

			total += 1 - totp;

			double q = total * rnd::GetRandom().Uniform();
			int mode = 0;
			while ( (mode<K0) && (q > cumul[mode])) mode++;
			if (mode == K0)	{
				double r = (q - (1 - totp)) / total;
				while (r > 0)	{
					if (mode >= Ncomponent)	{
						MESSAGE gimmemore = GIVEMEMORE;
						MPI_Send(gimmemore,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
						
						MPI_Recv(&mode,1,MPI_INT,0,TAG1,MPI_COMM_WORLD,&stat);
						MPI_Recv(&Ncomponent,1,MPI_INT,0,TAG1,MPI_COMM_WORLD,&stat);
						MPI_Recv(weight,Ncomponent,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);
						r = 0;
					}
					else	{
						r -= weight[mode];
						mode ++;
					}
				}
			}

			int Accepted = (mode != bk);
			if (Accepted)	{
				NAccepted ++;
			}
			alloc[site] = mode;
		}

		// send incmovedone message
		MESSAGE done = REALLOC_DONE;
		MPI_Send(done,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);

		// send back new allocations 
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
	
	return ((double) NAccepted) / GetNsite() / nrep;
}

double MatrixFiniteProfileProcess::IncrementalFiniteMove(int nrep, int K0)	{

	int NAccepted = 0;

	vector<double> mLogSamplingArray;
	vector<double> cumul;

	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		for (int site=0; site<GetNsite(); site++)	{

			int bk = alloc[site];

			double max = 0;
			for (int mode = 0; mode<K0; mode++)	{
				mLogSamplingArray[mode] =  LogStatProb(site,mode);
				if ((!mode) || (max < mLogSamplingArray[mode]))	{
					max = mLogSamplingArray[mode];
				}
			}

			double total = 0;
			for (int mode = 0; mode<K0; mode++)	{
				double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
				total += p;
				cumul[mode] = total;
			}

			total += 1 - totp;

			double q = total * rnd::GetRandom().Uniform();
			int mode = 0;
			while ( (mode<K0) && (q > cumul[mode])) mode++;
			if (mode == K0)	{
				double r = (q - (1 - totp)) / total;
				while (r > 0)	{
					mode++;
					if (mode >= Ncomponent)	{
						if (mode >= GetNmodeMax())	{
							cerr << "nmode max overflow\n";
							exit(1);
						}
						CreateComponent(Ncomponent);
						// create weight;
						Ncomponent++;
					}
					r -= weight[mode];
				}
			}

			int Accepted = (mode != bk);
			if (Accepted)	{
				NAccepted ++;
			}
			alloc[site] = mode;
		}
	}
	
	UpdateModeProfileSuffStat();
	return ((double) NAccepted) / GetNsite() / nrep;
}
