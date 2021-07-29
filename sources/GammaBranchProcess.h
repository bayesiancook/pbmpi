
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GAMMABRANCH_H
#define GAMMABRANCH_H

#include "TaxonSet.h"
#include "BranchProcess.h"

class GammaBranchProcess : public virtual BranchProcess	{

	public:

	GammaBranchProcess() : branchempalpha(0), branchempbeta(0), betaprior(0) {}
	virtual ~GammaBranchProcess() {}

	double LogBranchLengthPrior(const Branch* branch);

	double LogHyperPrior();

	double GetBranchAlpha() {return branchalpha;}
	double GetBranchBeta() {return branchbeta;}

	void SampleLength();
	void SampleLength(const Branch* branch);
    void PriorSampleLength();

	void ToStreamWithLengths(ostream& os, const Link* from);

	// conjugate sampling
	double Move(double tuning = 1, int nrep=1);
	double MoveBranchBeta(double tuning, int nrep);
	double MoveLength();

	double NonMPIMove(double tuning = 1, int nrep=1);
	double NonMPIMoveLength();

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual void Create(Tree* intree, double inalpha = 1, double inbeta = 10);

	virtual void Delete();

	double branchalpha;
	double branchbeta;

    double* branchempalpha;
    double* branchempbeta;

	int betaprior;
};

#endif

