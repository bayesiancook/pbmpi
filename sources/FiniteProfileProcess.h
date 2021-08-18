
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef FinitePROFILE_H
#define FinitePROFILE_H

#include <cmath>
#include "MixtureProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class FiniteProfileProcess: public virtual MixtureProfileProcess	{

	public:

	FiniteProfileProcess(int K = 1) : weight(0), fixncomp(false), empmix(false), Ncat(0), nmodemax(refnmodemax), statfix(0), empweight(0), dirweightprior(0), empdirweight(0)  {
		Ncomponent = K;
		weightalpha = 1;
	}
	virtual ~FiniteProfileProcess(){}

	// uniform prior on component number
	double GetLogNPrior() {return 0;}

	void SetFixedNcomponent(bool in = true)	{
		fixncomp = in;
	}

    void SetProfileFrac(double infrac)    {
        profilefrac = infrac;
        if (fixncomp && (Ncomponent == 1))  {
            for (int k=0; k<GetDim(); k++)  {
                dirweight[k] = profilefrac + (1-profilefrac)*empdirweight[k];
            }
        }
    }

	virtual int GetNmodeMax() {return fixncomp ? Ncomponent : nmodemax;} 
	virtual void SetNmodeMax(int n) {nmodemax = n;}

	protected:

	virtual void DrawProfileFromPrior();

	void ReadNcomponent(string name);
	void ReadStatFix(string name);
	void SetStatFix();

	virtual double IncrementalFiniteMove(int nrep) = 0;
	virtual double MoveHyper(double tuning, int nrep); // added virtual
	virtual double MoveWeightAlpha(double tuning, int nrep); // added virtual
	double MoveNcomponent(int nrep);
	virtual void ResampleWeights(); // added virtual

	// static allocation of many component-specific variables
	// such as: profiles, occupancy number
	// basically everything except substitution matrices


	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create(int innsite, int indim)	{
		cerr << "in create 2 arguments\n";
		exit(1);
	}

	virtual void Create(int innsite, int indim, int ncat, int infixncomp = 0, int inempmix = 0, string inmixtype = "None");
	virtual void Delete();

	// multinomial 
	void SampleAlloc();
	void SampleStat();
	virtual void SampleHyper(); // added virtual
	virtual void SampleWeights(); // added virtual

    virtual void PriorSampleHyper();

	virtual double LogHyperPrior();

	double LogWeightPrior();

	double* weight;
	double weightalpha;

	bool fixncomp;
	bool empmix;
	
	int Ncat;

    int nmodemax;

	string mixtype;
	double** statfix;
	double* empweight;

	int dirweightprior;
	// 0 : flexible
	// 1 : rigid
    
    double* empdirweight;

};

#endif

