
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRPARTPROFILE_H
#define PARTGTRPARTPROFILE_H

#include "PartitionedGTRProfileProcess.h"
#include "PartitionedProfileProcess.h"

class PartitionedGTRPartitionedProfileProcess : public virtual PartitionedGTRProfileProcess, public virtual PartitionedProfileProcess {

	public:

	PartitionedGTRPartitionedProfileProcess() : matrixarray(0) {}
	virtual ~PartitionedGTRPartitionedProfileProcess() {}

	SubMatrix* GetMatrix(int site)	{
		if (! matrixarray[PartitionedGTRProfileProcess::GetSitePart(site)][PartitionedProfileProcess::GetSitePart(site)]) 	{
			cerr << "error in get matrix : null partitioned matrix\n";
			exit(1);
		}
		return matrixarray[PartitionedGTRProfileProcess::GetSitePart(site)][PartitionedProfileProcess::GetSitePart(site)];
	}

	protected:

	// called at the beginning and the end of the run
	virtual void Create(int indim, PartitionScheme rrscheme, PartitionScheme statscheme);
	virtual void Delete();

	// simply creates/deletes GTR matrices for all currently existing components
	void CreateMatrix(int p, int k);
	virtual void UpdateMatrix(int p, int k);

	GTRSubMatrix* GetGTRMatrix(int p, int k)	{
		GTRSubMatrix* tmp = dynamic_cast<GTRSubMatrix*>(matrixarray[p][k]);
		if (!tmp)	{
			cerr << "error in GetGTRMatrix: null matrix \n";
			exit(1);
		}
		return tmp;
	}

	double GetNormRate(int p, int k)	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i+1; j<GetDim(); j++)	{
				tot += rr[p][rrindex(i,j,GetDim())] * profile[k][i] * profile[k][j];
			}
		}
		return 2*tot;
	}


	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	double GlobalMoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	void SlaveMoveProfile();
	double MoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	double MoveProfile(int cat, double tuning, int n, int nrep);


	virtual void UpdateModeProfileSuffStat() = 0;

	// should be called each time global parameters are modified
	virtual void UpdateMatrices()	{
		for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
			if(IsPartitionUnmasked(p))
			{
				for (int k=0; k<PartitionedProfileProcess::GetNpart(); k++)	{
					UpdateMatrix(p, k);
				}
			}
		}
	}

	virtual void CreateMatrices()	{
		for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
			if(IsPartitionUnmasked(p))
			{
				for (int k=0; k<PartitionedProfileProcess::GetNpart(); k++)	{
					if (! matrixarray[p][k])	{
						CreateMatrix(p, k);
					}
				}
			}
		}
		// check to delete unoccupied combinations
		/*
		for (int k=GetNcomponent(); k<GetNmodeMax(); k++)	{
			for (int p=0; p<PartitionedProfileProcess::GetNpart(); p++)	{
				if (matrixarray[k])	{
					DeleteMatrix(k, p);
				}
				matrixarray[k] = 0;
			}
		}*/
	}

	virtual void DeleteMatrices()	{
		for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
			if(IsPartitionUnmasked(p))
			{
				for (int k=0; k<PartitionedProfileProcess::GetNpart(); k++)	{
					DeleteMatrix(p, k);
				}
			}
		}
	}

	virtual void DeleteMatrix(int p, int k)	{
		delete matrixarray[p][k];
		matrixarray[p][k] = 0;
	}

	SubMatrix*** matrixarray;
};

#endif
