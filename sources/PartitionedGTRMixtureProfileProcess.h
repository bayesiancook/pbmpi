
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef PARTGTRMIXTUREPROFILE_H
#define PARTGTRMIXTUREPROFILE_H


#include "PartitionedGTRProfileProcess.h"
#include "MixtureProfileProcess.h"

// general superclass for GTR-like Dirichlet-process mixture on profiles
class PartitionedGTRMixtureProfileProcess : public virtual PartitionedGTRProfileProcess, public virtual MixtureProfileProcess {

	public:

	PartitionedGTRMixtureProfileProcess() : matrixarray(0) {}
	virtual ~PartitionedGTRMixtureProfileProcess() {}

	SubMatrix* GetMatrix(int site)	{
		if (! matrixarray[GetSitePart(site)][alloc[site]]) 	{
			cerr << "error in get matrix : null matrix , site " << site << " and alloc : " << alloc[site] << '\n';
			exit(1);
		}
		return matrixarray[GetSitePart(site)][alloc[site]];
	}

	protected:

	virtual void Create(int indim, PartitionScheme rrscheme);
	virtual void Delete();

	void SwapComponents(int cat1, int cat2);

	// simply creates/deletes GTR matrices for all currently existing components
	void CreateMatrix(int p, int k);
	virtual void UpdateMatrix(int p, int k);

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
		for (int p=0; p<GetNpart(); p++)	{
			if(!IsPartitionMasked(p))
			{
				for (int k=0; k<GetNcomponent(); k++)	{
					UpdateMatrix(p, k);
				}
			}
		}
	}

	virtual void CreateMatrices()	{
		for (int p=0; p<GetNpart(); p++)	{
			if(!IsPartitionMasked(p))
			{
				for (int k=0; k<GetNcomponent(); k++)
				{
					if (! matrixarray[p][k])
					{
						CreateMatrix(p, k);
					}
				}

				for (int k=GetNcomponent(); k<GetNmodeMax(); k++)	{
					if (matrixarray[p][k])	{
						DeleteMatrix(p, k);
					}
					matrixarray[p][k] = 0;
				}
			}
		}
	}

	virtual void DeleteMatrices()	{
		for (int p=0; p<GetNpart(); p++)	{
			if(!IsPartitionMasked(p))
			{
				for (int k=0; k<GetNcomponent(); k++)	{
					DeleteMatrix(p, k);
				}
			}
		}
	}

	virtual void DeleteMatrix(int p, int k)	{
		delete matrixarray[p][k];
		matrixarray[p][k] = 0;
	}

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

	SubMatrix*** matrixarray;

};

#endif

