
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MPIMODULE_H
#define MPIMODULE_H

class MPIModule {

	protected:

	MPIModule() {}
	virtual ~MPIModule() {}

	virtual int GetNprocs() = 0;
	virtual int GetMyid() = 0;
	virtual int GetSiteMin(int proc = -1) = 0;
	virtual int GetSiteMax(int proc = -1) = 0;
	virtual void MakeMPIPartition() = 0;

};


#endif

