
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMatrixOneProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GeneralPathSuffStatMatrixOneProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixOneProfileProcess::UpdateProfileSuffStat()	{

	profilepaircount.clear();
	profilerootcount.clear();
	profilewaitingtime.clear();

	for (int i=0; i<GetNsite(); i++)	{
		
		map<pair<int,int>, int>& paircount = GetSitePairCount(i);
		map<int,double>& waitingtime = GetSiteWaitingTime(i);
		int rootstate = GetSiteRootState(i);

		profilerootcount[rootstate]++;
		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			profilewaitingtime[i->first] += i->second;
		}
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			profilepaircount[i->first] += i->second;
		}
	}
}

double GeneralPathSuffStatMatrixOneProfileProcess::ProfileSuffStatLogProb()	{
	double total = 0;
	SubMatrix* mat = matrix;
	if (! mat)	{
		cerr << "error : null matrix\n";
		exit(1);
	}
	const double* stat = matrix->GetStationary();
	for (map<int,int>::iterator i = profilerootcount.begin(); i!= profilerootcount.end(); i++)	{
		total += i->second * log(stat[i->first]);
	}
	for (map<int,double>::iterator i = profilewaitingtime.begin(); i!= profilewaitingtime.end(); i++)	{
		total += i->second * (*mat)(i->first,i->first);
	}
	for (map<pair<int,int>, int>::iterator i = profilepaircount.begin(); i!= profilepaircount.end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	return total;
}

