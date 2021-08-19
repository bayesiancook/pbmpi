
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MIXTUREPROFILE_H
#define MIXTUREPROFILE_H

#include <cmath>
#include "ProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class MixtureProfileProcess: public virtual ProfileProcess	{

	public:

	MixtureProfileProcess() : profile(0) {}
	virtual ~MixtureProfileProcess(){}

	double* GetProfile(int site)	{
		return profile[alloc[site]];
	}

	int GetNcomponent() { return Ncomponent;}
	virtual int GetNDisplayedComponent()	{
		return Ncomponent;
	}
	int GetNOccupiedComponent() {
		int n = 0;
		for (int k=0; k<Ncomponent; k++)	{
			if (occupancy[k])	{
				n++;
			}
		}
		return n;
	}

    double GetEffectiveComponentNumber()    {
        double m1 = 0;
        double m2 = 0;
        for (int k=0; k<Ncomponent; k++)    {
            double tmp = ((double) occupancy[k]) / GetNsite();
            m1 += tmp;
            m2 += tmp*tmp;
        }
        if (fabs(m1 - 1) > 1e-6)    {
            cerr << "error: total relative occupancies do not sum to 1\n";
            cerr << m1 << '\n';
            exit(1);
        }
        return 1.0 / m2;
    }

	virtual int GetNmodeMax() {return GetNsite();}
	// virtual int GetNmodeMax() {return 100;}

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt(alloc[site]);}

	// summary statistic: mean entropy over all profiles
	double GetStatEnt();
	double GetStatEnt(int k);
	double GetCenterStatEnt();

	double GetMeanDirWeight();

	void RenormalizeProfiles();

	// generic Move function
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1) = 0;
	double MoveDirWeights(double tuning, int nrep);

	protected:

	virtual void UpdateModeProfileSuffStat() = 0;

	// implements a pure virtual defined in ProfileProcess
	double ProfileSuffStatLogProb();

	// suffstat lnL of site <site> when allocated to component <cat>
	virtual double LogStatProb(int site, int cat) = 0;

	// the component suff stat log prob is yet to be implemented in subclasses
	virtual double ProfileSuffStatLogProb(int cat) = 0;

	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create(int innsite, int indim);
	virtual void Delete();

	// in certain models,
	// the matrix associated to each component should be created on the fly
	// all other component-specific variables are static (see above, NmodeMax)
	virtual void CreateComponent(int k) = 0;
	// virtual void CreateComponent(int k, double* instat) = 0;
	virtual void DeleteComponent(int k) = 0;
	virtual void UpdateComponent(int k) = 0;

	void UpdateComponents()	{
		for (int k=0; k<GetNcomponent(); k++)	{
			UpdateComponent(k);
		}
	}

	virtual double GetAllocEntropy()	{
		double total = 0;
		UpdateOccupancyNumbers();
		for (int k=0; k<GetNcomponent(); k++)	{
			double tmp = ((double) occupancy[k]) / GetNsite();
			if (tmp)	{
				total -= tmp * log(tmp);
			}
		}
		return total;
	}

	virtual void AddSite(int site, int cat)	{
		alloc[site] = cat;
		occupancy[cat]++;
	}
	virtual void RemoveSite(int site, int cat)	{
		occupancy[cat]--;
	}

	// sample all aspects of the mixture (number of components, composition) from the prior
	virtual void SampleProfile();

    virtual void PriorSampleProfile()   {
        PriorSampleHyper();
        PriorSampleAlloc();
        PriorSampleStat();
    }

	virtual void SampleHyper() = 0;
	virtual void SampleAlloc() = 0;
	virtual void SampleStat();
	void SampleStat(int cat);

    virtual void PriorSampleHyper() = 0;
    virtual void PriorSampleAlloc() {
        SampleAlloc();
    }
    virtual void PriorSampleStat()  {
       SampleStat();
    }

	void SampleStat(double* stat, double statmin = 0);
    void SampleEmpiricalStat(double* stat, const double* count, double statsmin = 0);

	void UpdateOccupancyNumbers();
	double ResampleEmptyProfiles();

	virtual double LogProfilePrior();

	virtual double LogHyperPrior() = 0;
	virtual double LogAllocPrior()	{
		return 0;
	}

	// dirichlet prior ~ Dirichlet (dirweight[0]... dirweight[GetDim()-1])
	virtual double LogStatPrior();

	virtual double LogStatPrior(int cat);
    double LogStatPrior(const double* profile);
    double EmpiricalLogStatPrior(const double* profile, const double* count);

	virtual void SwapComponents(int cat1, int cat2);

    virtual void SetEmpiricalDirWeightPrior(double* inalpha, double* inbeta);

	double** profile;
	double* allocprofile;
	double* dirweight;
	int* alloc;
	int* occupancy;
	int Ncomponent;
	double* logstatprior;
	double* profilesuffstatlogprob;
    double* empdirweightalpha;
    double* empdirweightbeta;
};

#endif

