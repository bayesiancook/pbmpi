
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GTRMIXTUREPROFILE_H
#define GTRMIXTUREPROFILE_H


#include "GTRProfileProcess.h"
#include "MatrixMixtureProfileProcess.h"

// general superclass for GTR-like Dirichlet-process mixture on profiles
class GTRMixtureProfileProcess : public virtual GTRProfileProcess, public virtual MatrixMixtureProfileProcess {

	public:

	GTRMixtureProfileProcess() {}
	virtual ~GTRMixtureProfileProcess() {}

	protected:

	virtual void Create(int innsite, int indim);
	virtual void Delete();

	virtual void SampleProfile()    {
        if (! fixrr)    {
            SampleRR();
        }
        MixtureProfileProcess::SampleProfile();
    }

    virtual void PriorSampleProfile()   {
        if (! fixrr)    {
            PriorSampleRR();
        }
        MixtureProfileProcess::PriorSampleProfile();
    }

    virtual double LogProfilePrior()    {
        double ret = 0;
        if (! fixrr)    {
            ret += LogRRPrior();
        }
        ret += MixtureProfileProcess::LogProfilePrior();
        return ret;
    }

	// simply creates/deletes GTR matrices for all currently existing components
	void CreateMatrix(int k);
	virtual void UpdateMatrix(int k);

	GTRSubMatrix* GetGTRMatrix(int k)	{
		GTRSubMatrix* tmp = dynamic_cast<GTRSubMatrix*>(matrixarray[k]);
		if (!tmp)	{
			cerr << "error in GetGTRMatrix: null matrix \n";
			exit(1);
		}
		return tmp;
	}

	double GetNormRate(int k)	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i+1; j<GetDim(); j++)	{
				tot += rr[rrindex(i,j,GetDim())] * profile[k][i] * profile[k][j];
			}
		}
		return 2*tot;
	}

	virtual double GetNormalizationFactor()	{
		UpdateOccupancyNumbers();
		double norm = 0;
		int tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (occupancy[k])	{
				norm += (occupancy[k] + 1) * GetNormRate(k);
				tot += occupancy[k] + 1;
			}
		}
		/*
		if (tot != GetNsite() + GetNcomponent())	{
			cerr << "error in norm factor\n";
			exit(1);
		}
		*/
		norm /= tot;
		return norm;
	}

};

#endif

