
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "OneProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* OneProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void OneProfileProcess::Create(int innsite, int indim)	{
	if (! profile)	{
		ProfileProcess::Create(innsite,indim);
		profile = new double[GetDim()];
	}
}

void OneProfileProcess::Delete()	{
	if (profile)	{
		delete[] profile;
		profile = 0;
		ProfileProcess::Delete();
	}
}

double OneProfileProcess::GetStatEnt()	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		total -= profile[i] * log(profile[i]);
	}
	return  total;
}

void OneProfileProcess::SampleProfile()	{
	SampleStat();
}

void OneProfileProcess::SampleStat()	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[k] = rnd::GetRandom().sGamma(1.0);
		total += profile[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[k] /= total;
	}
	total = 0;
	for (int k=0; k<GetDim(); k++)	{
		if (profile[k] < stateps)	{
			profile[k] = stateps;
		}
		total += profile[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[k] /= total;
	}
	// UpdateComponent(i);
}

double OneProfileProcess::LogProfilePrior()	{
	return 0;
}


