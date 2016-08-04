
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GeneralPathSuffStatMatrixMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixMixtureProfileProcess::Create(int innsite, int indim)	{
	if (! profilepaircount)	{
		MatrixMixtureProfileProcess::Create(innsite,indim);
		profilepaircount = new map< pair<int,int>, int>[GetNmodeMax()];
		profilewaitingtime = new map<int,double>[GetNmodeMax()];
		profilerootcount = new map<int,int>[GetNmodeMax()];
	}
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::Delete() {
	if (profilepaircount)	{
		delete[] profilepaircount;
		delete[] profilerootcount;
		delete[] profilewaitingtime;
		profilepaircount = 0;
		profilerootcount = 0;
		profilewaitingtime = 0;
		MatrixMixtureProfileProcess::Delete();
	}
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::UpdateModeProfileSuffStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		profilepaircount[i].clear();
		profilerootcount[i].clear();
		profilewaitingtime[i].clear();
	}
	for (int i=0; i<GetNsite(); i++)	{
		
		map<pair<int,int>, int>& paircount = GetSitePairCount(i);
		map<int,double>& waitingtime = GetSiteWaitingTime(i);
		int rootstate = GetSiteRootState(i);
		/*
		if (rootstate < 0)	{
			cerr << "error : negative root state \n";
			cerr << rootstate << '\n';
			exit(1);
		}
		*/
		int cat = alloc[i];
		if (rootstate != -1)	{
			profilerootcount[cat][rootstate]++;
			for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
				profilewaitingtime[cat][i->first] += i->second;
			}
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				profilepaircount[cat][i->first] += i->second;
			}
		}
	}
}

double GeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	SubMatrix* mat = matrixarray[cat];
	if (! mat)	{
		cerr << "error : null matrix\n";
		cerr << cat << '\t' << Ncomponent << '\n';
		cerr << occupancy[cat] << '\n';
		exit(1);
	}
	const double* stat = matrixarray[cat]->GetStationary();
	for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
		total += i->second * log(stat[i->first]);
	}
	for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
		total += i->second * (*mat)(i->first,i->first);
	}
	for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MatrixMixtureProfileProcess::SwapComponents(cat1,cat2);

	// booh! inefficient copies!
	//cerr << "swap in general path suff stat matrix mixture profile process\n";
	//exit(1);

	map<int,int> roottmp = profilerootcount[cat1];
	profilerootcount[cat1] = profilerootcount[cat2];
	profilerootcount[cat2] = roottmp;

	map<int,double> timetmp = profilewaitingtime[cat1];
	profilewaitingtime[cat1] = profilewaitingtime[cat2];
	profilewaitingtime[cat2] = timetmp;

	map<pair<int,int>,int> pairtmp = profilepaircount[cat1];
	profilepaircount[cat1] = profilepaircount[cat2];
	profilepaircount[cat2] = pairtmp;
}


double GeneralPathSuffStatMatrixMixtureProfileProcess::LogStatProb(int site, int cat)	{
	double total = 0;
	SubMatrix* mat = matrixarray[cat];
	const double* stat = matrixarray[cat]->GetStationary();
	if (GetSiteRootState(site) != -1)	{
		total += log(stat[GetSiteRootState(site)]);

		map<int,double>& waitingtime = GetSiteWaitingTime(site);
		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			total += i->second * (*mat)(i->first,i->first);
		}

		map<pair<int,int>, int>& paircount = GetSitePairCount(site);
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			total += i->second * log((*mat)(i->first.first, i->first.second));
		}
	}
	return total;
}


void GeneralPathSuffStatMatrixMixtureProfileProcess::AddSite(int site, int cat)	{
	alloc[site] = cat;
	occupancy[cat] ++;
		
	/*
	if (activesuffstat)	{
		profilerootcount[cat][GetSiteRootState(site)]++;

		map<int,double>& waitingtime = GetSiteWaitingTime(site);
		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			profilewaitingtime[cat][i->first] += i->second;
		}

		map<pair<int,int>, int>& paircount = GetSitePairCount(site);
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			profilepaircount[cat][i->first] += i->second;
		}
	}
	*/
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::RemoveSite(int site, int cat)	{
	occupancy[cat] --;
		
	/*
	if (activesuffstat)	{
		profilerootcount[cat][GetSiteRootState(site)]--;

		map<int,double>& waitingtime = GetSiteWaitingTime(site);
		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			profilewaitingtime[cat][i->first] -= i->second;
		}

		map<pair<int,int>, int>& paircount = GetSitePairCount(site);
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			profilepaircount[cat][i->first] -= i->second;
		}
	}
	*/
}
