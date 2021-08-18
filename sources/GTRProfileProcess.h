
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GTRPROFILE_H
#define GTRPROFILE_H

#include "GTRSubMatrix.h"
#include "MatrixProfileProcess.h"

// superclass for all GTR-like models
class GTRProfileProcess : public virtual MatrixProfileProcess {

	public:

	GTRProfileProcess() : rr(0), fixrr(false), rrtype("None"), emprralpha(0), emprrbeta(0) {}
	virtual ~GTRProfileProcess() {}

	int GetNrr()	{
		return Nrr;
	}

	const double* GetRR() {
		if (! rr)	{
			cerr << "error : getrr\n";
			exit(1);
		}
		return rr;
	}

	double GetRRMean()	{
		double mean = 0;
		for (int i=0; i<Nrr; i++)	{
			mean += rr[i];
		}
		mean /= Nrr;
		return mean;
	}

	void SetRR(string type);

	double GetRRVarCoeff()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<Nrr; i++)	{
			mean += rr[i];
			var += rr[i] * rr[i];
		}
		mean /= Nrr;
		var /= Nrr;
		var -= mean*mean;
		var /= mean * mean;
		return var;
	}

	double GetRREntropy()	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			total += rr[i];
		}
		double ent = 0;
		for (int i=0; i<Nrr; i++)	{
			double tmp = rr[i] / total;
			if (tmp > 1e-6)	{
				ent -= tmp * log(tmp);
			}
		}
		return ent;
	}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	protected:

	virtual void Create(int innsite, int indim);
	virtual void Delete();

	// relative rates
	virtual double LogRRPrior();
	void SampleRR();
    void PriorSampleRR();

	// assumes that site-specific sufficient statistics are already updated
	// collect them into more compact sufficient statistics
	// (depending on what is permitted by the type of sufficient satistics used)
	virtual void UpdateRRSuffStat() = 0;
	virtual void GlobalUpdateRRSuffStat() = 0;
	virtual void SlaveUpdateRRSuffStat() = 0;

	// Metropolis or Gibbs Sampling algorithm,
	// (depending on what is permitted by the type of sufficient satistics used)
	virtual void MoveRR() = 0;

	int Nrr;
	double* rr;

	bool fixrr;
	string rrtype;
    double* emprralpha;
    double* emprrbeta;

};

#endif

