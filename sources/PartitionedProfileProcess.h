
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTPROFILE_H
#define PARTPROFILE_H

#include <cmath>
#include "Partition.h"
#include "ProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class PartitionedProfileProcess: public virtual ProfileProcess, public PartitionProcess	{

	public:

	PartitionedProfileProcess() : profile(0) {}
	virtual ~PartitionedProfileProcess(){}

	double* GetProfile(int site)	{
		return profile[GetSitePart(site)];
	}

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt(GetSitePart(site));}

	// summary statistic: mean entropy over all profiles
	double GetStatEnt();
	double GetStatEnt(int k);
	double GetCenterStatEnt();

	double GetMeanDirWeight();

	void RenormalizeProfiles();

	// generic Move function
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1) = 0;
	virtual double MoveHyper(double tuning, int nrep); // added virtual
	double MoveDirWeights(double tuning, int nrep);

	protected:

	void SetStat(int inpart, string type);

	virtual void UpdateModeProfileSuffStat() = 0;

	// implements a pure virtual defined in ProfileProcess
	double ProfileSuffStatLogProb();

	// suffstat lnL of site <site> when allocated to component <cat>
	virtual double LogStatProb(int site, int cat) = 0;

	// the component suff stat log prob is yet to be implemented in subclasses
	virtual double ProfileSuffStatLogProb(int cat) = 0;

	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create(int indim, PartitionScheme statscheme);
	virtual void Delete();

	// in certain models,
	// the matrix associated to each component should be created on the fly
	// all other component-specific variables are static (see above, NmodeMax)
	virtual void CreateComponent(int k) = 0;
	// virtual void CreateComponent(int k, double* instat) = 0;
	virtual void DeleteComponent(int k) = 0;
	virtual void UpdateComponent(int k) = 0;

	void UpdateComponents()	{
		for (int k=0; k<GetNpart(); k++)	{
			UpdateComponent(k);
		}
	}

	// sample all aspects of the mixture (number of components, composition) from the prior
	void SampleProfile();

	virtual void SampleHyper();
	virtual void SampleStat();
	void SampleStat(int cat);
	void SampleStat(double* stat, double statmin = 0);

	double LogProfilePrior();

	virtual double LogHyperPrior();

	// dirichlet prior ~ Dirichlet (dirweight[0]... dirweight[GetDim()-1])
	virtual double LogStatPrior();

	virtual double LogStatPrior(int cat);

	double** profile;
	double* allocprofile;
	double* dirweight;
	double* logstatprior;
	double* profilesuffstatlogprob;

	bool* fixstat;

	int nfreestat;
};

#endif

