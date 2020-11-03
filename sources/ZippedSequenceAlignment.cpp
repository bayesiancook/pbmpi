
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ZippedSequenceAlignment.h"
#include <iostream>


void ZippedSequenceAlignment::CreateZipArrays()	{

	Data = new int*[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		Data[i] = new int[Nsite];
	}

	Indices = new int*[Nsite];
	ZipIndices = new int*[Nsite];
	for (int i=0; i<Nsite; i++)	{
		Indices[i] = new int[Nstate];
		ZipIndices[i] = new int[Nstate];
	}

	ZipSize = new int[Nsite];

	OrbitSize = new int[Nsite];
	Orbit = new bool*[Nsite];
	for (int i=0; i<Nsite; i++)	{
		Orbit[i] = new bool[Nstate];
	}
}

void ZippedSequenceAlignment::DeleteZipArrays()	{

	for (int i=0; i<Ntaxa; i++)	{
        delete[] Data[i];
    }
    delete[] Data;

	for (int i=0; i<Nsite; i++)	{
        delete[] ZipIndices[i];
        delete[] Indices[i];
	}
    delete[] ZipIndices;
    delete[] Indices;

    delete[] ZipSize;

    delete[] OrbitSize;

	for (int i=0; i<Nsite; i++)	{
        delete[] Orbit[i];
	}
    delete[] Orbit;
}

void ZippedSequenceAlignment::ComputeZipArrays()	{

	for (int i=0; i<Nsite; i++)	{

		OrbitSize[i] = 0;
		for (int k=0 ; k< Nstate; k++)	{
			Orbit[i][k]= false;
		}

		for (int j=0; j<Ntaxa; j++)	{
			int d = GetTemplate()->GetState(j,i);
			if (d != unknown)	{
				if (! Orbit[i][d])	{
					Orbit[i][d] = true;
					Indices[i][OrbitSize[i]] = d;
					OrbitSize[i] ++;
				}
			}
		}

		// sort Indices[i]
		for (int j=0; j<OrbitSize[i]; j++)	{
			for (int k=OrbitSize[i]-1; k>j; k--)	{
				if (Indices[i][j] > Indices[i][k])	{
					int tmp = Indices[i][j];
					Indices[i][j] = Indices[i][k];
					Indices[i][k] = tmp;
				}
			}
		}

		// reverse translation table
		for (int j=0; j<OrbitSize[i]; j++)	{
			ZipIndices[i][Indices[i][j]] = j;
		}
		for (int j=0; j<Nstate; j++)	{
			if (! Orbit[i][j])	{
				ZipIndices[i][j] = OrbitSize[i];
			}
		}	

		if (OrbitSize[i] < Nstate)	{
			ZipSize[i] = OrbitSize[i] + 1;
		}
		else	{
			ZipSize[i] = OrbitSize[i];
		}

		for (int j=0; j<Ntaxa; j++)	{
			Data[j][i] = -2;
			int d = GetTemplate()->GetState(j,i);
			if (d == unknown)	{
				Data[j][i] = unknown;
			}
			else	{
				for (int k=0; k<OrbitSize[i]; k++)	{
					if (Indices[i][k] == d)	{
						Data[j][i] = k;
					}
				}

				// Zip identity:
				// Indices[i][Data[j][i]] == template->Data[j][i] , for every i and j
				// within their respective range

			}

			// here, may be check that ZipData != -2

			if (Data[j][i] == -2)	{
				cerr << "error in zip data making\n";
			}
		}
		
		int observed[Nstate];
		int orbitsize= OrbitSize[i];
		for (int k=0; k<Nstate; k++)	{
			observed[k] = 0;
		}
		for (int k=0; k<orbitsize; k++)	{
			observed[Indices[i][k]] = 1;
		}
		for (int k=0; k<Nstate; k++)	{
			if (! observed[k])	{
				Indices[i][orbitsize++] = k;
			}
		}
		if (orbitsize != Nstate)	{
			cerr << "error in RegisterWithData\n";
			cerr << "site : " << i << '\n';
			cerr << '\n';
			exit(1);
		}
	}

	double temp = 0;
	for (int i=0; i<Nsite; i++)	{
		temp += ((double) ZipSize[i] * ZipSize[i]) / Nstate / Nstate;
	}
	SpeedFactor =  temp / Nsite ;
	int nconst = 0;
	for (int i=0; i<Nsite; i++)	{
		if (OrbitSize[i] == 1)	{
			nconst++;
		}
	}
	MeanOrbitSize = 0;
	for (int i=0; i<Nsite; i++)	{
		MeanOrbitSize += OrbitSize[i];
	}
	MeanOrbitSize /= Nsite;

}

