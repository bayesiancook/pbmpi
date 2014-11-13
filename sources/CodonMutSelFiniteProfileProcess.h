
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSELFINITEPROFILE_H
#define CODONMUTSELFINITEPROFILE_H

#include "MatrixFiniteProfileProcess.h"
#include "CodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class CodonMutSelFiniteProfileProcess : public virtual MatrixFiniteProfileProcess, public virtual CodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	CodonMutSelFiniteProfileProcess() {}
	virtual ~CodonMutSelFiniteProfileProcess() {}

	protected:

	void Create(int innsite, int indim)	{
		cerr << "In two-argument Create of AACodonMutSelFiniteProfileProcess. Should not be here.\n";
		exit(1);
	}

	void Create(int innsite, int indim, int ncat, int infixncomp, int inempmix, string inmixtype, CodonStateSpace* instatespace)	{
		MatrixFiniteProfileProcess::Create(innsite,indim,ncat,infixncomp,inempmix,inmixtype);
		//MatrixFiniteProfileProcess::Create(innsite,indim,ncat);
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create(innsite,indim);
		CodonMutSelProfileProcess::Create(innsite,indim,instatespace);
	}
	
	void Delete()	{
		CodonMutSelProfileProcess::Delete();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		MatrixFiniteProfileProcess::Delete();
	}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		GlobalUpdateSiteProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{

			// mutation rates
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			MoveNucRR(tuning,2);
			MoveNucStat(tuning,2);
			MoveNucRR(tuning*0.9,2);
			MoveNucStat(tuning*0.9,2);
			MoveNucRR(tuning*0.8,2);
			MoveNucStat(tuning*0.8,2);
			
			MoveNucRR(tuning*0.7,2);
			MoveNucStat(tuning*0.7,2);
			MoveNucRR(tuning*0.6,2);
			MoveNucStat(tuning*0.6,2);
			MoveNucRR(tuning*0.5,2);
			MoveNucStat(tuning*0.5,2);
			
			MoveNucRR(tuning*0.1,2);
			MoveNucStat(tuning*0.1,2);
			MoveNucRR(tuning*0.09,2);
			MoveNucStat(tuning*0.09,2);
			MoveNucRR(tuning*0.08,2);
			MoveNucStat(tuning*0.08,2);

			MoveNucRR(tuning*0.07,2);
			MoveNucStat(tuning*0.07,2);
			MoveNucRR(tuning*0.06,2);
			MoveNucStat(tuning*0.06,2);
			MoveNucRR(tuning*0.05,2);
			MoveNucStat(tuning*0.05,2);

			// allocations
			if (Ncomponent != GetNsite())	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				GlobalIncrementalFiniteMove(3);
			}
			if (!fixncomp)	{
				MoveNcomponent(10);
			}

			if (!empmix)	{
				// profiles
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				UpdateModeProfileSuffStat();
				GlobalMoveProfile(1,1,100);
				GlobalMoveProfile(1,3,100);
				GlobalMoveProfile(0.1,3,100);

				// hyperparameters
				if (dirweightprior == 0)	{
					MoveHyper(tuning,10);
				}
			}		
		}
		return 1;
	}

	void ToStream(ostream& os) {
		for (int i=0; i<Nnuc; i++)	{
			os << GetNucStat(i) << '\t';
		}
		os << '\n';
		os << '\n';
		for (int i=0; i<GetNnucrr(); i++)	{
			os << GetNucRR(i) << '\t';
		}
		os << '\n';
		os << '\n';		

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

	void FromStream(istream& is) {
		for (int i=0; i<Nnuc; i++)	{
			is >> nucstat[i];
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			is >> nucrr[i];
		}

		is >> Ncomponent;
		for (int j=0; j<GetDim(); j++)	{
			is >> dirweight[j];
		}
		for (int i=0; i<Ncomponent; i++)	{
			for (int j=0; j<GetDim(); j++)	{
				is >> profile[i][j];
			}
		}
		for (int i=0; i<GetNsite(); i++)	{
			is >> alloc[i];
		}
	}


	//void ProfilesToStream(ostream& os) {
	//	for (int i=0; i<GetNsite(); i++)	{
	//		for (int j=0; j<GetDim(); j++)	{
	//			os << profile[alloc[i]][j] << '\t';
	//		}
	//		os << '\n';
	//	}	
	//}

	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelFiniteProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new CodonMutSelProfileSubMatrix(statespace,nucrr,nucstat,profile[k],true);
		//matrixarray[k] = new CodonMutSelProfileSubMatrix(statespace,nucrr,nucstat,profile[k],false);
	}

	void UpdateMatrix(int k)	{
		matrixarray[k]->CorruptMatrix();
	}
	/*
	void CreateMatrices()	{
		cerr << "in codon mut sel create matrices\n";
		exit(1);
	}

	void DeleteMatrices()	{
		cerr << "in codon mut sel delete matrices\n";
		exit(1);
	}
	*/
	
};

#endif

