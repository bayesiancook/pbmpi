
/********************

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"

#include "Random.h"
#include "PhyloProcess.h"
#include <string>

#include <cassert>
#include "Parallel.h"

#include "TexTab.h"

bool PhyloProcess::ActiveSite(int i)  {
    return (minsitecutoff == -1) || ((steppingrank[i] >= minsitecutoff) && (steppingrank[i] < maxsitecutoff));
}


void PhyloProcess::GlobalPrepareStepping()   {
	assert(myid == 0);
    if (bkdata) {
        cerr << "error in SlaveBackupData: bkdata already exists\n";
        exit(1);
    }

	MESSAGE signal = PREPARESTEPPING;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    bkdata = new SequenceAlignment(GetData());

    steppingrank = new int[GetNsite()];
    rnd::GetRandom().DrawFromUrn(steppingrank, GetNsite(), GetNsite());
    for (int i=0; i<GetNsite(); i++)    {
        steppingrank[i] = i;
    }
	MPI_Bcast(steppingrank,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlavePrepareStepping()	{
    if (bkdata) {
        cerr << "error in SlaveBackupData: bkdata already exists\n";
        exit(1);
    }
    bkdata = new SequenceAlignment(GetData());
    steppingrank = new int[GetNsite()];
	MPI_Bcast(steppingrank,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalSetSteppingFraction(int cutoff1, int cutoff2)    {
	MESSAGE signal = SETSTEPPINGFRAC;
    int cutoff[2];
    cutoff[0] = cutoff1;
    cutoff[1] = cutoff2;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(cutoff,2,MPI_INT,0,MPI_COMM_WORLD);
    SetSteppingFraction(cutoff1, cutoff2);
}

void PhyloProcess::SlaveSetSteppingFraction()    {
    int cutoff[2];
	MPI_Bcast(cutoff,2,MPI_INT,0,MPI_COMM_WORLD);
    SetSteppingFraction(cutoff[0], cutoff[1]);
}

void PhyloProcess::SetSteppingFraction(int cutoff1, int cutoff2)  {

    for (int i=0; i<GetNsite(); i++)    {
        if ((steppingrank[i] >= cutoff1) && (steppingrank[i] < cutoff2))  {
            for (int j=0; j<GetNtaxa(); j++)    {
                GetData()->SetState(j, i, bkdata->GetState(j,i));
            }
        }
        else    {
            for (int j=0; j<GetNtaxa(); j++)    {
                GetData()->SetState(j, i, -1);
            }
        }
    }
    minsitecutoff = cutoff1;
    maxsitecutoff = cutoff2;
    // fill missing map
    FillMissingMap();
}

void PhyloProcess::GlobalSetEmpiricalFrac(double infrac)    {
	MESSAGE signal = EMPIRICALFRAC;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&infrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    SetEmpiricalFrac(infrac);
}

void PhyloProcess::SlaveSetEmpiricalFrac()  {
    double frac;
	MPI_Bcast(&frac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    SetEmpiricalFrac(frac);
}
