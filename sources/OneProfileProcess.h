
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ONEPROFILE_H
#define ONEPROFILE_H

#include <cmath>
#include "ProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class OneProfileProcess: public virtual ProfileProcess	{

	public:

	OneProfileProcess() : profile(0) {}
	virtual ~OneProfileProcess(){}

	double* GetProfile(int site)	{
		return profile;
	}

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt();}

	int GetNOccupiedComponent()  {return 1;}

	double GetStatEnt();

	// generic Move function
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1) = 0;

	protected:

	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create(int innsite, int indim);
	virtual void Delete();

	void UpdateProfile()	{
	}


	// sample all aspects of the mixture (number of components, composition) from the prior
	void SampleProfile();
	void SampleStat();

	double LogProfilePrior();

	double* profile;
};

#endif

