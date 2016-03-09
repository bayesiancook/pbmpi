
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRPROFILE_H
#define PARTGTRPROFILE_H

#include "GTRSubMatrix.h"
#include "MatrixProfileProcess.h"
#include "Partition.h"

// superclass for all GTR-like models
class PartitionedGTRProfileProcess : public virtual MatrixProfileProcess, public virtual PartitionProcess {

	public:

	PartitionedGTRProfileProcess() : rr(0) {}
	virtual ~PartitionedGTRProfileProcess() {}

	int GetNrr()	{
		return Nrr;
	}

	const double* GetRR(int inpart) {
		if (! rr)	{
			cerr << "error : getrr\n";
			exit(1);
		}
		return rr[inpart];
	}

	double GetRRMean()	{
		double grandMean = 0.0;
		for(int p = 0; p < GetNpart(); p++)
		{
			double mean = 0;
			for (int i=0; i<Nrr; i++)	{
				mean += rr[p][i];
			}
			mean /= Nrr;

			grandMean += mean;
		}
		return grandMean/GetNpart();
	}

	void SetRR(int inpart, string type);

	double GetRRVarCoeff(int inpart)	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<Nrr; i++)	{
			mean += rr[inpart][i];
			var += rr[inpart][i] * rr[inpart][i];
		}
		mean /= Nrr;
		var /= Nrr;
		var -= mean*mean;
		var /= mean * mean;
		return var;
	}

	double GetRREntropy()	{
		double mean_ent = 0.0;
		for(int p = 0; p < GetNpart(); p++)
		{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				total += rr[p][i];
			}
			double ent = 0;
			for (int i=0; i<Nrr; i++)	{
				double tmp = rr[p][i] / total;
				if (tmp > 1e-6)	{
					ent -= tmp * log(tmp);
				}
			}

			mean_ent += ent;
		}
		return mean_ent/GetNpart();
	}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	protected:

	virtual void Create(int indim, PartitionScheme inscheme);
	virtual void Delete();

	// relative rates
	virtual double LogRRPrior();
	virtual void SampleRR();

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
	double** rr;

	bool* fixrr;

	int nfreerr;

};

#endif

