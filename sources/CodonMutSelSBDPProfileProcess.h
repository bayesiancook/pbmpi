
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSELSBDPPROFILE_H
#define CODONMUTSELSBDPPROFILE_H

#include "MatrixSBDPProfileProcess.h"
#include "CodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class CodonMutSelSBDPProfileProcess : public virtual MatrixSBDPProfileProcess, public virtual CodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	CodonMutSelSBDPProfileProcess() {}
	virtual ~CodonMutSelSBDPProfileProcess() {}

	protected:

	void Create(int innsite, int indim)	{
		cerr << "In two-argument Create of CodonMutSelSBDPProfileProcess. Should not be here.\n";
		exit(1);
	}

	void Create(int innsite, int indim, CodonStateSpace* instatespace)	{
		MatrixSBDPProfileProcess::Create(innsite,indim);
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create(innsite,indim);
		CodonMutSelProfileProcess::Create(innsite,indim,instatespace);
	}
	
	void Delete()	{
		CodonMutSelProfileProcess::Delete();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		MatrixSBDPProfileProcess::Delete();
	}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		//GlobalUpdateSiteProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{

			// mutation rates
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			MoveNucRR(tuning*0.5, 2);
			MoveNucRR(tuning*0.4, 2);
			MoveNucRR(tuning*0.3, 2);
			MoveNucRR(tuning*0.2, 2);
			MoveNucRR(tuning*0.1, 2);
			MoveNucRR(tuning*0.05, 2);
			MoveNucRR(tuning*0.04, 2);
			MoveNucRR(tuning*0.03, 2);
			MoveNucRR(tuning*0.02, 2);
			MoveNucRR(tuning*0.01, 2);
			MoveNucRR(tuning*0.05, n);
			MoveNucRR(tuning*0.04, n);
			MoveNucRR(tuning*0.03, n);
			MoveNucRR(tuning*0.02, n);
			MoveNucRR(tuning*0.01, n);
			MoveNucStat(tuning*0.5,n);
			MoveNucStat(tuning*0.4,n);
			MoveNucStat(tuning*0.3,n);
			MoveNucStat(tuning*0.2,n);
			MoveNucStat(tuning*0.1,n);
			MoveNucStat(tuning*0.05,n);
			MoveNucStat(tuning*0.04,n);
			MoveNucStat(tuning*0.03,n);
			MoveNucStat(tuning*0.02,n);
			MoveNucStat(tuning*0.01,n);

			// allocations
			//GlobalUpdateParameters();
			//GlobalUpdateSiteProfileSuffStat();//*/
			//IncrementalDPMove(3,GetNsite() / 2);
			// GlobalIncrementalFiniteMove(3);
			//MoveNcomponent(10);

			// label switch moves
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			GlobalMixMove(5,1,0.001,40);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);

			// profiles
			/*
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			GlobalMoveProfile(1,1,100);
			GlobalMoveProfile(1,3,100);
			GlobalMoveProfile(0.1,3,100);
			*/

			// hyperparameters
			//GlobalUpdateParameters();
			//GlobalUpdateSiteProfileSuffStat();
			MoveHyper(tuning,10);		
			MoveHyper(tuning*0.5,10);		
			MoveHyper(tuning*0.4,10);		
			MoveHyper(tuning*0.3,10);		
			MoveHyper(tuning*0.2,10);		
			MoveHyper(tuning*0.1,10);		
		}
		return 1;
	}

	void ToStream(ostream& os)	{
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

		os << kappa << '\n';
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
	void FromStream(istream& is)	{
		for (int i=0; i<Nnuc; i++)	{
			is >> nucstat[i];
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			is >> nucrr[i];
		}

		is >> kappa;
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


	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in CodonMutSelSBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		//matrixarray[k] = new AAMutSelProfileSubMatrix(statespace,nucrr,nucstat,profile[k],false);
		matrixarray[k] = new CodonMutSelProfileSubMatrix(statespace,nucrr,nucstat,profile[k],true);
	}
	
	virtual void SwapComponents(int cat1, int cat2)	{
		MatrixSBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	void UpdateMatrix(int k)	{
		matrixarray[k]->CorruptMatrix();
	}
};

#endif

