
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonFiniteProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonFiniteProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PoissonFiniteProfileProcess::ToStream(ostream& os)	{

	os << Ncomponent << '\n';
	for (int j=0; j<GetDim(); j++)	{
		os << dirweight[j] << '\t';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			os << profile[i][j] << '\t';
		}
		os << '\n';
	}

	for (int i=0; i<GetNsite(); i++)	{
		os << alloc[i] << '\t';
	}
	os << '\n';
}

void PoissonFiniteProfileProcess::FromStream(istream& is)	{

	is >> Ncomponent;
	
	for (int i=0; i<GetDim(); i++)	{
		is >> dirweight[i];
	}

	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			is >> profile[i][j];
		}
	}

	for (int i=0; i<GetNsite(); i++)	{
		is >> alloc[i];
	}

	ResampleWeights();
	// CHECK some update here ?
}


double PoissonFiniteProfileProcess::GlobalIncrementalFiniteMove(int nrep)	{

	assert(GetMyid() == 0);

	// send command and arguments
	MESSAGE signal = REALLOC_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

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

		// resample component weights based on current site allocations
		// and send them to slaves
		UpdateOccupancyNumbers();
		ResampleWeights();
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// receive new site allocations from slave
		MPI_Status stat;
		int tmpalloc[GetNsite()];
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
	}
	
	// final cleanup
	UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();

	// CHECK that: might be useful depending on the exact submodel
	// ResampleWeights();
	// UpdateMatrices();

	return ((double) NAccepted) / GetNsite() / nrep;
}

double PoissonFiniteProfileProcess::SlaveIncrementalFiniteMove()	{

	// parse argument sent by master
	int nrep;
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int NAccepted = 0;

	double* bigarray = new double[Ncomponent * GetNsite()];
	double* bigcumul = new double[Ncomponent * GetNsite()];

	for (int rep=0; rep<nrep; rep++)	{

		// receive weights sent by master
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// do the incremental reallocation move on my site range
		for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

			double* mLogSamplingArray = bigarray + site * Ncomponent;
			double* cumul = bigcumul + site * Ncomponent;

			int bk = alloc[site];

			double max = 0;
			for (int mode = 0; mode < Ncomponent; mode++)	{
				mLogSamplingArray[mode] =  LogStatProb(site,mode);
				if ((!mode) || (max < mLogSamplingArray[mode]))	{
					max = mLogSamplingArray[mode];
				}
			}

			double total = 0;

			for (int mode = 0; mode < Ncomponent; mode++)	{
				double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
				total += p;
				cumul[mode] = total;
			}

			double q = total * rnd::GetRandom().Uniform();
			int mode = 0;
			while ( (mode<Ncomponent) && (q > cumul[mode])) mode++;
			if (mode == Ncomponent)	{
				cerr << "error in switch mode: gibbs overflow\n";
				exit(1);
			}

			int Accepted = (mode != bk);
			if (Accepted)	{
				NAccepted ++;
			}
			alloc[site] = mode;
		}

		// send new allocations to master
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
	
	delete[] bigarray;
	delete[] bigcumul;
	return ((double) NAccepted) / GetNsite() / nrep;
}

