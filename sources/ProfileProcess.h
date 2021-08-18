
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PROFILE_H
#define PROFILE_H

#include "Chrono.h"
#include "StateSpace.h"

#include <iostream>
#include <map>

using namespace std;

// this class takes care of all aspects related to substitution processes (not including branch lengths and site-specific relative rates)
// this includes global parameters (such as global relative exchangeabilities for CATGTR, or mutational parameters for Codon Mutation Selection models)
// and hyperparameters of the mixtures of the mixture models
//
// this code allows for variations across sites only of the profiles
// all other parameters should be homogeneous along the sequence

const double stateps = 1e-100;
const int refnmodemax = 1000;

class ProfileProcess {

	public:

	ProfileProcess() : nsite(0), dim(0), activesuffstat(false), burnin(false), statinfcount(0), totstatcount(0), mintotweight(0), profilefrac(1.0) {}
	virtual ~ProfileProcess() {}

	virtual string GetVersion() = 0;

    virtual void SetProfileFrac(double infrac)  {
        profilefrac = infrac;
    }

	double GetStatInfCount() {
		double tmp = ((double) statinfcount) / totstatcount;
		statinfcount = 0;
		totstatcount = 0;
		return tmp;
	}

	// dimension of the profile(s)
	int GetDim() {return dim;}
	int GetNsite() {return nsite;}
	virtual double* GetProfile(int site) = 0;

	// virtual int GetNOccupiedComponent()  = 0;
	virtual StateSpace* GetStateSpace() = 0;

	// total prior associated to all aspects of the substitution processes across sites
	virtual double LogProfilePrior() = 0;

	// sample all parameters 
	virtual void SampleProfile() = 0;
	virtual void PriorSampleProfile() = 0;

	virtual double GetMeanStationaryEntropy() = 0;
	virtual double GetSiteStationaryEntropy(int site) = 0;

	// implemented in specialized phyloprocess classes
	// will collect the sufficient statistics necessary for updating substitution parameters
	// separately for each site
	// these statistics will be accessed through specific functions declared and implemented in subclasses
	virtual void UpdateSiteProfileSuffStat() = 0;
	virtual void GlobalUpdateSiteProfileSuffStat() = 0;
	virtual void SlaveUpdateSiteProfileSuffStat() = 0;

	virtual void GlobalUpdateParameters() = 0;
	virtual void SlaveUpdateParameters() = 0;

	/*
	virtual void SendCurrentProfileConfig();
	virtual void ReceiveCurrentProfileConfig();
	*/

	// returns that part of the total likelihood that depends on substitution process parameters (not including branch lengths and site rates)
	virtual double ProfileSuffStatLogProb() = 0;

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;

	virtual void SetBurnin(bool in)	{
		burnin = in;
	}

	virtual int isObserved(int site, int k) = 0;

	double GetMinStat(int site)	{

		double* profile = GetProfile(site);
		double min = 1;
		for (int k=0; k<GetDim(); k++)	{
			if (isObserved(site,k))	{
				if (min > profile[k])	{
					min = profile[k];
				}
			}
		}
		return min;
	}

	double GetMinMinStat()	{
		double min = 1;
		for (int i=0; i<GetNsite(); i++)	{
			double tmp = GetMinStat(i);
			if (min > tmp)	{
				min = tmp;
			}
		}
		return min;
	}

	protected:

	virtual void DrawProfileFromPrior() {
		cerr << "error: in ProfileProcess::DrawProfileFromPrior\n";
		exit(1);
	}

	virtual double GetNormalizationFactor()	{cerr << "should not be here\n"; exit(1); return 1;}

	// called at the beginning of the run only
	// but potentially called several times due to multiple inheritance
	// thus, should internally control that no object is created twice
	// (by checking that pointers are null before creating)
	virtual void Create(int innsite, int indim);

	// called at the end of the run only
	// but potentially called several times due to multiple inheritance
	// should check that pointers are not zero before deleting
	virtual void Delete() {}

	virtual int GetNprocs() = 0;
	virtual int GetMyid() = 0;
	virtual int GetSiteMin() = 0;
	virtual int GetSiteMax() = 0;

	virtual double* GetEmpiricalFreq() = 0;

	int nsite;
	int dim;
	bool activesuffstat;

	double GetMinTotWeight() {return mintotweight;}
	void SetMinTotWeight(double in) {
		mintotweight = in;
		/*
		cerr << "limit : " << in << '\n';
		if (in < 0)	{
			cerr << "dim : " << GetDim() << '\n';
			mintotweight = ((double) GetDim()) / 4;
		}
		else	{
			mintotweight = in;
		}
		cerr << "min tot weight : " << mintotweight << '\n';
		*/
	}

	Chrono chronodp;
	Chrono chronostat;
	Chrono chronorr;

	bool burnin;

	int statinfcount;
	int totstatcount;

	double mintotweight;

    double profilefrac;
};


#endif
