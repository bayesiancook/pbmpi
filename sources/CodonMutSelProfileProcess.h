
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSEL_H
#define CODONMUTSEL_H

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"

class CodonMutSelProfileProcess : public virtual GeneralPathSuffStatMatrixProfileProcess	{

	// y mettre les variables globales (taux de mutation essentiellement)

	// s'inspirer de GTRProfileProcess et GeneralPathSuffStatGTRProfileProcess

	public:

	CodonMutSelProfileProcess() : nucrr(0), nucstat(0), statespace(0) {}
	virtual ~CodonMutSelProfileProcess() {}

	int GetNnucrr()	{
		return Nnucrr;
	}	

	const double* GetNucRR()	{
		if (! nucrr)	{
			cerr << "error : getnucrr\n";
			exit(1);
		}
		return nucrr;
	}

	const double GetNucRR(int i)	{
		if ( (! nucrr) || (i>GetNnucrr()) || (i<0) )	{
			cerr << "error : getnucrr i\n";
			exit(1);
		}
		return nucrr[i];
	}

	const double* GetNucStat()	{
		if (! nucstat)	{
			cerr << "error : getnucstat\n";
			exit(1);
		}
		return nucstat;
	}

	const double GetNucStat(int i)	{
		if ( (! nucstat) || (i>Nnuc) || (i<0) )	{
			cerr << "error : getnucstat i\n";
			exit(1);
		}
		return nucstat[i];
	}


	protected:

	virtual void Create(int innsite, int indim, CodonStateSpace* instatespace);
	virtual void Delete();
	virtual double GetNormalizationFactor()	{return 1.0;}

	// nuc relative rates
	virtual double LogNucRRPrior();
	virtual void SampleNucRR();

	// nuc stationaries
	virtual double LogNucStatPrior();
	virtual void SampleNucStat();

	double MoveNucRR(double tuning); 
	double MoveNucRR(double tuning, int n); 
	double MoveNucStat(double tuning, int n);
	
	int Nnucrr;
	double* nucrr;
	double* nucstat;
	CodonStateSpace* statespace;
};

#endif
