
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ZIP_H
#define ZIP_H

#include "SequenceAlignment.h"

class ZippedSequenceAlignment : public SequenceAlignment	{

	public:

	ZippedSequenceAlignment(SequenceAlignment* infrom)	{

		from = infrom;

		Ntaxa = from->Ntaxa;
		Nsite = from->Nsite;
		taxset = from->taxset;
		statespace = from->statespace;
		Nstate = statespace->GetNstate();

        CreateZipArrays();
		ComputeZipArrays();
	}

	~ZippedSequenceAlignment()	{
        DeleteZipArrays();
	}

	SequenceAlignment* GetTemplate() {return from;}
	int GetZipSize(int site) {return ZipSize[site];}
	int GetOrbitSize(int site) {return OrbitSize[site];}
	int GetStateFromZip(int site, int state) {
		if (state >= OrbitSize[site])	{
			cerr << "error in getstate from zip\n";
			exit(1);
		}
		return Indices[site][state];
	}
	bool InOrbit(int site, int state)	{
		return Orbit[site][state];
	}

	void ComputeZipArrays();

	private:

    void CreateZipArrays();
    void DeleteZipArrays();

	SequenceAlignment* from;
	int Nstate;

	double SpeedFactor;

	bool** Orbit;
	int* OrbitSize;
	double MeanOrbitSize;
	int* ZipSize;
	int** Indices;
	int** ZipIndices;
	int** ZipData;

	
	

};


#endif

