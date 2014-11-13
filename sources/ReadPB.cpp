
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"
#include "Model.h"


int main(int argc, char* argv[])	{

	int myid,nprocs;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int blockcounts[2] = {1,3};
	MPI_Datatype types[2] = {MPI_DOUBLE,MPI_INT};
	MPI_Aint dtex,displacements[2];
	
	displacements[0] = (MPI_Aint) 0;
	MPI_Type_extent(MPI_DOUBLE,&dtex);
	displacements[1] = dtex;
	MPI_Type_struct(2,blockcounts,displacements,types,&Propagate_arg);
	MPI_Type_commit(&Propagate_arg); 

	string name = argv[argc-1];
	
	Model* model = new Model(name,myid,nprocs);
	if (myid == 0) {
		model->ReadPB(argc,argv);
		MESSAGE signal = KILL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		model->WaitLoop();
	}

	MPI_Finalize();
}

