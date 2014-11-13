
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AAMUTSELMIXTUREPROFILE_H
#define AAMUTSELMIXTUREPROFILE_H

#include "AAMutSelProfileProcess.h"
#include "MixtureProfileProcess.h"

class AAMutSelMixtureProfileProcess : public virtual AAMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRMixtureProfileProcess

	public:

	AAMutSelMixtureProfileProcess() {}
	virtual ~AAMutSelMixtureProfileProcess() {}

	protected:

	void Create(int innsite, int indim)	{
		cerr << "In two-argument Create of AAMutSelMixtureProfileProcess. Should not be here.\n";
		exit(1);
	}

	void Create(int innsite, int indim, CodonStateSpace* instatespace)	{
		AAMutSelProfileProcess::Create(innsite,indim,instatespace);
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create(innsite,indim);
	}
	
	void Delete()	{
		AAMutSelProfileProcess::Delete();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
	}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		GlobalUpdateSiteProfileSuffStat();
		for (int rep=0; rep<nrep; rep++)	{
			UpdateModeProfileSuffStat();
			//MoveNucRR(tuning);
			//MoveNucRR(tuning*0.5);
			//MoveNucRR(tuning*0.1);
			MoveNucRR(tuning, 2);
			MoveNucRR(tuning, 3);
			MoveNucRR(tuning*0.5, 2);
			MoveNucRR(tuning*0.5, 3);
			MoveNucRR(tuning*0.1, 2);
			MoveNucRR(tuning*0.1, 3);
			MoveNucStat(tuning,n);
			MoveNucStat(tuning*0.5,n);
			MoveNucStat(tuning*0.1,n);
			MoveProfile(1,1,100);
			MoveProfile(1,3,100);
			MoveProfile(0.1,3,100);
			// MoveHyper(tuning,10);		
			IncrementalMixtureMove(5);
			// MoveNcomponent(100);
		}
		return 1;
	}


	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AAMutSelMixtureProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AAMutSelProfileSubMatrix(statespace,nucrr,nucstat,profile[k],false);
		//matrixarray[k] = new AAMutSelProfileSubMatrix(statespace,nucrr,nucstat,profile[k],true);
		// matrix is normalized for testing...
	}

	void UpdateMatrix(int k)	{
		matrixarray[k]->CorruptMatrix();
	}

	void CreateMatrices()	{}
	void DeleteMatrices()	{}
	
	
};

#endif

