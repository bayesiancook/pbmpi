
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

void PhyloProcess::GlobalSetSteppingFraction(int cutoff)    {
	MESSAGE signal = SETSTEPPINGFRAC;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&cutoff,1,MPI_INT,0,MPI_COMM_WORLD);
    SetSteppingFraction(cutoff);
}

void PhyloProcess::SlaveSetSteppingFraction()    {
    int cutoff;
	MPI_Bcast(&cutoff,1,MPI_INT,0,MPI_COMM_WORLD);
    SetSteppingFraction(cutoff);
}

void PhyloProcess::SetSteppingFraction(int cutoff)  {

    for (int i=0; i<GetNsite(); i++)    {
        if (steppingrank[i] < cutoff)  {
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
    // fill missing map
    FillMissingMap();
}

