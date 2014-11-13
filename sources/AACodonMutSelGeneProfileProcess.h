
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELGENEPROFILEPROCESS_H
#define AACODONMUTSELGENEPROFILEPROCESS_H

#include "GeneProfileProcess.h"
#include "AACodonMutSelProfileProcess.h"

class AACodonMutSelGeneProfileProcess : public virtual AACodonMutSelProfileProcess, public virtual GeneProfileProcess	{

	public:

	AACodonMutSelGeneProfileProcess() : aacmutselhub(0) {}
	virtual ~AACodonMutSelGeneProfileProcess() {}

	void SetMixtureParameters()	{
		for (int i=0; i<Nnuc; i++)	{
			nucstat[i] = aacmutselhub->GetNucStat(i);	
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			nucrr[i] = aacmutselhub->GetNucRR(i);
		}
		for (int i=0; i<statespace->GetNstate(); i++)	{
			codonprofile[i] = aacmutselhub->GetCodonProfileEntry(i);
		}
		*omega = aacmutselhub->GetOmega();
		/*
		for (int i=0; i<GetNsite(); i++)	{
			alloc[i] = expohub->alloc[i+offset];
		}
		*/
	}

	protected:

	virtual void Create(int innsite, int innstate, CodonStateSpace* instatespace, int infixbl, int infixcodonprofile, int infixomega, AACodonMutSelProfileProcess* inhub, int inoffset)	{
		AACodonMutSelProfileProcess::Create(innsite,innstate,instatespace);
		GeneProfileProcess::Create(innsite,innstate,inhub,inoffset);
		aacmutselhub = inhub;
		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		fixbl = infixbl;
	}

	virtual void Delete()	{
		GeneProfileProcess::Delete();
		AACodonMutSelProfileProcess::Delete();
	}

	virtual SubMatrix* GetMatrix(int site)	{
		return aacmutselhub->GetMatrix(site+offset);
	}

	virtual void UpdateMatrices() {
		cerr << "error: in gene update matrices\n";
		exit(1);
	}

	virtual void CreateMatrices()	{
		// cerr << "errror : in gene create matrices\n";
		// exit(1);
	}

	virtual void DeleteMatrices()	{
		// cerr << "errror : in gene delete matrices\n";
		// exit(1);
	}

	/*
	virtual const int* GetSiteProfileSuffStatCount(int site)	{
		return expohub->GetSiteProfileSuffStatCount(site+offset);
	}

	virtual const double* GetSiteProfileSuffStatBeta(int site)	{
		return expohub->GetSiteProfileSuffStatBeta(site+offset);
	}
	*/

	AACodonMutSelProfileProcess* aacmutselhub;
	int fixcodonprofile;
	int fixomega;
	int fixbl;

};


#endif

