
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef PARTEXPCONGTRPARTPROFILE_H
#define PARTEXPCONGTRPARTPROFILE_H

#include "PartitionedGTRPartitionedProfileProcess.h"
#include "PartitionedExpoConjugateGTRProfileProcess.h"

// Exponential conjugate GTR models
class PartitionedExpoConjugateGTRPartitionedProfileProcess : public virtual PartitionedGTRPartitionedProfileProcess, public virtual PartitionedExpoConjugateGTRProfileProcess {

	public:

	PartitionedExpoConjugateGTRPartitionedProfileProcess() : profilesuffstatcount(0), profilesuffstatbeta(0) {}
	virtual ~PartitionedExpoConjugateGTRPartitionedProfileProcess() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		totchrono.Start();
		// GlobalUpdateParameters();
		// GlobalUpdateSiteProfileSuffStat();
		// UpdateModeProfileSuffStat();
		if(nfreestat + nfreerr > 0)
		{
			for (int rep=0; rep<nrep; rep++)	{

				// relative rates
				if(nfreerr > 0)
				{
					GlobalUpdateParameters();
					GlobalUpdateSiteProfileSuffStat();
					// useless in an expo suff stat context
					// UpdateModeProfileSuffStat();
					MoveRR();
				}

				// profiles
				if(nfreestat > 0)
				{
					GlobalUpdateParameters();
					GlobalUpdateSiteProfileSuffStat();
					UpdateModeProfileSuffStat();
					profilechrono.Start();
					GlobalMoveProfile(1,1,100);
					GlobalMoveProfile(1,3,100);
					GlobalMoveProfile(0.1,3,100);
					profilechrono.Stop();

					if(nfreestat > 1)
						MoveHyper(tuning,10);
				}
			}
		}
		totchrono.Stop();
		return 1;
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	// matrices are not used during Sufficient-statistics based updates
	// all matrices are deleted upon 'Collapsing' (pruning->suffstat)
	// and re-created upon 'Unfolding' (suffstat->pruning) state-vectors
	// using CreateMatrices() and DeleteMatrices()
	virtual void CreateComponent(int k) {
		if (activesuffstat)	{
			for (int l=0; l<GetDim(); l++)	{
				profilesuffstatcount[k][l] = 0;
				profilesuffstatbeta[k][l] = 0;
			}
		}
		else	{
		}
		SampleStat(k);
	}
	virtual void DeleteComponent(int k) {
	}
	virtual void UpdateComponent(int k) {
		// update suff stats ?
	};

	double ProfileSuffStatLogProb(int cat);
	virtual double LogStatProb(int site, int cat);

	virtual void Create(int innsite, int indim)	{
		cerr << "in create 2 arguments\n";
		exit(1);
	}

	virtual void Create(int indim, PartitionScheme rrscheme, PartitionScheme statscheme);

	virtual void Delete();

	// collects site-specific suffstats and pools them componentwise
	void UpdateModeProfileSuffStat();

	// component-specific sufficient statistics
	int** profilesuffstatcount;
	double** profilesuffstatbeta;

	Chrono totchrono;
	Chrono profilechrono;

};

#endif

