
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSBDPPROFILE_H
#define AACODONMUTSELSBDPPROFILE_H

#include "MatrixSBDPProfileProcess.h"
#include "AACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class AACodonMutSelSBDPProfileProcess : public virtual MatrixSBDPProfileProcess, public virtual AACodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	AACodonMutSelSBDPProfileProcess() {}
	virtual ~AACodonMutSelSBDPProfileProcess() {}

	protected:

	void Create(int innsite, int indim)	{
		cerr << "In two-argument Create of AACodonMutSelSBDPProfileProcess. Should not be here.\n";
		exit(1);
	}

	void Create(int innsite, int indim, CodonStateSpace* instatespace, int infixcodonprofile, int infixomega)	{
		MatrixSBDPProfileProcess::Create(innsite,indim);
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create(innsite,indim);
		AACodonMutSelProfileProcess::Create(innsite,indim,instatespace);
		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
	}
	
	void Delete()	{
		AACodonMutSelProfileProcess::Delete();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		MatrixSBDPProfileProcess::Delete();
	}

	/*
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		//GlobalUpdateSiteProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{

			// mutation rates
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			if (! fixcodonprofile)	{
				//MoveCodonProfile(tuning,30,10);
				MoveNucStatCodonProfile(tuning,30,10);
				MoveNucRR(tuning,2);
				//MoveCodonProfile(tuning*0.5,30,10);
				MoveNucStatCodonProfile(tuning*0.5,30,10);
				MoveNucRR(tuning*0.5,2);
				//MoveCodonProfile(tuning*0.1,30,10);
				MoveNucStatCodonProfile(tuning*0.1,30,10);
				MoveNucRR(tuning*0.1,2);
				//MoveCodonProfile(tuning*0.01,30,10);
				MoveNucStatCodonProfile(tuning*0.01,30,10);
				MoveNucRR(tuning*0.01,2);
				//MoveCodonProfile(tuning*0.001,30,10);
				MoveNucStatCodonProfile(tuning*0.001,30,10);
				MoveNucRR(tuning*0.001,2);
			}
			else	{
				MoveNucRR(tuning,2);
				MoveNucStat(tuning,2);
				//MoveNucRR(tuning*0.8,2);
				//MoveNucStat(tuning*0.8,2);
				//MoveNucRR(tuning*0.6,2);
				//MoveNucStat(tuning*0.6,2);
				//MoveNucRR(tuning*0.4,2);
				//MoveNucStat(tuning*0.4,2);
				//MoveNucRR(tuning*0.2,2);
				//MoveNucStat(tuning*0.2,2);
			
				MoveNucRR(tuning*0.1,2);
				MoveNucStat(tuning*0.1,2);
				//MoveNucRR(tuning*0.08,2);
				//MoveNucStat(tuning*0.08,2);
				//MoveNucRR(tuning*0.06,2);
				//MoveNucStat(tuning*0.06,2);
				//MoveNucRR(tuning*0.04,2);
				//MoveNucStat(tuning*0.04,2);
				//MoveNucRR(tuning*0.02,2);
				//MoveNucStat(tuning*0.02,2);

				MoveNucRR(tuning*0.01,2);
				MoveNucStat(tuning*0.01,2);
				//MoveNucRR(tuning*0.008,2);
				//MoveNucStat(tuning*0.008,2);
				//MoveNucRR(tuning*0.006,2);
				//MoveNucStat(tuning*0.006,2);
				//MoveNucRR(tuning*0.004,2);
				//MoveNucStat(tuning*0.004,2);
				//MoveNucRR(tuning*0.002,2);
				//MoveNucStat(tuning*0.002,2);

				MoveNucRR(tuning*0.001,2);
				MoveNucStat(tuning*0.001,2);
				//MoveNucRR(tuning*0.0008,2);
				//MoveNucStat(tuning*0.0008,2);
				//MoveNucRR(tuning*0.0006,2);
				//MoveNucStat(tuning*0.0006,2);
				//MoveNucRR(tuning*0.0004,2);
				//MoveNucStat(tuning*0.0004,2);
				//MoveNucRR(tuning*0.0002,2);
				//MoveNucStat(tuning*0.0002,2);
			}

			if (! fixomega)	{
				MoveOmega(tuning*2);
				MoveOmega(tuning);
				MoveOmega(tuning*0.1);
				MoveOmega(tuning*0.01);
			}

			// allocations
			//GlobalUpdateParameters();
			//GlobalUpdateSiteProfileSuffStat();
			//IncrementalDPMove(3,GetNsite() / 2);
			// GlobalIncrementalFiniteMove(3);
			//MoveNcomponent(10);


			// label switch moves
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			//cerr << "GlobalMixMove\n";
			//cerr.flush();
			GlobalMixMove(5,1,0.001,40);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);


			// hyperparameters
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			if (dirweightprior == 0)	{
				MoveHyper(tuning,10);		
				MoveHyper(tuning*0.5,10);		
				MoveHyper(tuning*0.4,10);		
				MoveHyper(tuning*0.3,10);		
				MoveHyper(tuning*0.2,10);		
				MoveHyper(tuning*0.1,10);
			}
			else	{
				MoveKappa(tuning,10);
				MoveKappa(tuning*0.5,10);		
				MoveKappa(tuning*0.4,10);		
				MoveKappa(tuning*0.3,10);		
				MoveKappa(tuning*0.2,10);		
				MoveKappa(tuning*0.1,10);
			}
		}
		return 1;
	}
	*/

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		for (int rep=0; rep<nrep; rep++)	{

			// mutation rates
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			if (! fixcodonprofile)	{
				MoveCodonProfile(tuning,30,10);
				MoveNucStatCodonProfile(tuning,30,10);
				MoveNucRR(tuning,2);
				MoveCodonProfile(tuning*0.5,30,10);
				MoveNucStatCodonProfile(tuning*0.5,30,10);
				MoveNucRR(tuning*0.5,2);
				MoveCodonProfile(tuning*0.1,30,10);
				MoveNucStatCodonProfile(tuning*0.1,30,10);
				MoveNucRR(tuning*0.1,2);
				MoveCodonProfile(tuning*0.01,30,10);
				MoveNucStatCodonProfile(tuning*0.01,30,10);
				MoveNucRR(tuning*0.01,2);
				MoveCodonProfile(tuning*0.001,30,10);
				MoveNucStatCodonProfile(tuning*0.001,30,10);
				MoveNucRR(tuning*0.001,2);
			}
			else	{
				MoveNucRR(tuning,2);
				MoveNucStat(tuning,2);
			
				MoveNucRR(tuning*0.1,2);
				MoveNucStat(tuning*0.1,2);
			}

			if (! fixomega)	{
				MoveOmega(tuning);
				MoveOmega(tuning*0.3);
			}

			// label switch moves
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			GlobalMixMove(1,1,0.001,40);
			MoveOccupiedCompAlloc(2);
			MoveAdjacentCompAlloc(2);

			// hyperparameters
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			if (dirweightprior == 0)	{
				MoveHyper(tuning,10);		
				MoveHyper(tuning*0.5,10);		
				MoveHyper(tuning*0.4,10);		
				MoveHyper(tuning*0.3,10);		
				MoveHyper(tuning*0.2,10);		
				MoveHyper(tuning*0.1,10);
			}
			else	{
				MoveKappa(tuning,10);
				MoveKappa(tuning*0.5,10);		
				MoveKappa(tuning*0.4,10);		
				MoveKappa(tuning*0.3,10);		
				MoveKappa(tuning*0.2,10);		
				MoveKappa(tuning*0.1,10);
			}
		}
		return 1;
	}

	virtual double MoveMixtureOnly(double tuning = 1, int n = 1, int nrep = 1)	{
		//GlobalUpdateSiteProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{

			// label switch moves
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();

			//cerr << "GlobalMixMove\n";
			//cerr.flush();
			GlobalMixMove(5,1,0.001,40);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);


			// hyperparameters
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
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

		for (int i=0; i<AACodonMutSelProfileProcess::statespace->GetNstate(); i++)	{
			os << codonprofile[i] << '\t';
		}
		os << '\n';
		os << '\n';
		os << *omega << '\n';

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
		for (int i=0; i<AACodonMutSelProfileProcess::statespace->GetNstate(); i++)	{
			is >> codonprofile[i];
		}

		is >> *omega;
		is >> kappa;
		is >> Ncomponent;
		for (int j=0; j<GetDim(); j++)	{
			is >> dirweight[j];
		}
		for (int i=0; i<Ncomponent; i++)	{
            double tot = 0;
			for (int j=0; j<GetDim(); j++)	{
				is >> profile[i][j];
                tot += profile[i][j];
			}
            /*
            if (fabs(tot-1) > 1e-5) {
                cerr << "normalization error when reading profiles\n";
                cerr << tot-1 << '\n';
                exit(1);
            }
            */
			for (int j=0; j<GetDim(); j++)	{
                profile[i][j] /= tot;
            }
		}
		for (int i=0; i<GetNsite(); i++)	{
			is >> alloc[i];
		}
		ResampleWeights();
	}


	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelSBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AACodonMutSelProfileSubMatrix(statespace,nucrr,nucstat,codonprofile,profile[k],omega,true);
		//matrixarray[k] = new AACodonMutSelProfileSubMatrix(statespace,nucrr,nucstat,codonprofile,profile[k],omega,false);
	}

	virtual void SwapComponents(int cat1, int cat2)	{
		MatrixSBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	void UpdateMatrix(int k)	{
		matrixarray[k]->CorruptMatrix();
		// matrixarray[k]->UpdateMatrix();
	}

	
	GeneticCodeType codetype;
	int fixcodonprofile;
	int fixomega;

};

#endif

