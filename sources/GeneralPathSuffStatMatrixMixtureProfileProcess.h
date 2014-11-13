
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSMATMIXTUREPROFILE_H
#define GENPATHSSMATMIXTUREPROFILE_H

#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "MatrixMixtureProfileProcess.h"

// general subclass for all Matrix Mixture mixtures using generic sufficient statistics
class GeneralPathSuffStatMatrixMixtureProfileProcess : public virtual MatrixMixtureProfileProcess, public virtual GeneralPathSuffStatMatrixProfileProcess  {

	public:

	GeneralPathSuffStatMatrixMixtureProfileProcess() : profilerootcount(0), profilepaircount(0) , profilewaitingtime(0) {}
	virtual ~GeneralPathSuffStatMatrixMixtureProfileProcess() {}

	protected:

	virtual void Create(int innsite, int indim);
	virtual void Delete();

	virtual void CreateComponent(int k)	{
		occupancy[k] = 0;
		SampleStat(k);
		// useful?
		if (activesuffstat)	{
			profilepaircount[k].clear();
			profilerootcount[k].clear();
			profilewaitingtime[k].clear();
		}
		CreateMatrix(k);
		UpdateMatrix(k);
	}

	virtual void DeleteComponent(int k)	{
		DeleteMatrix(k);
	}

	// update component k (in particular, will be used for updating the matrix
	// typically, after the profile has changed)
	virtual void UpdateComponent(int k)	{
		UpdateMatrix(k);
	}

	// necessary to keep track of componentwise sufficient statistics 
	// when updating the structure of the mixture
	void AddSite(int site, int cat);
	void RemoveSite(int site, int cat);

	// collects site-specific suffstats and pools them componentwise
	void UpdateModeProfileSuffStat();

	double ProfileSuffStatLogProb(int cat);
	void SwapComponents(int cat1, int cat2);

	virtual double LogStatProb(int site, int cat);

	// implemented in phyloprocess
	// virtual double logSiteProbPath(int site, SubMatrix* mat) = 0;

	// componentwise
	map<int,int>* profilerootcount;
	map< pair<int,int>, int>* profilepaircount;
	map<int,double>* profilewaitingtime;

};

#endif

