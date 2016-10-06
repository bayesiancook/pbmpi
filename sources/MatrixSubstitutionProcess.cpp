
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixSubstitutionProcess.h"
#include <vector>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Matrix Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//	* Sample substitution mappings for all site, on a given branch
//	and  conditional on states at both ends of the branch
//	(CPU level 3)
//-------------------------------------------------------------------------

// root case (trivial)
BranchSitePath** MatrixSubstitutionProcess::SampleRootPaths(int* state)	{
	// BranchSitePath** patharray = new BranchSitePath*[sitemax - sitemin];
	BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		patharray[i] = new BranchSitePath(state[i]);
	}

	return patharray;
}

// general case
BranchSitePath** MatrixSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time) 	{
	// BranchSitePath** patharray = new BranchSitePath*[sitemax - sitemin];
	BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		double rate = GetRate(i);
		SubMatrix* matrix = GetMatrix(i);
		BranchSitePath* path = ResampleAcceptReject(1000,stateup[i],statedown[i],rate,time,matrix);
		if (! path)	{
			path = ResampleUniformized(stateup[i],statedown[i],rate,time,matrix);
		}
		patharray[i] = path;
	}
	return patharray;
}


// this is an implementation of Nielsen 2002, Syst Biol 51:729
// accept-reject sampling method for drawing a substitution mapping along a branch
// conditional on the states at both ends

BranchSitePath* MatrixSubstitutionProcess::ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix)	{

	int ntrial = 0;
	BranchSitePath* path = 0;
	if (rate * totaltime < 1e-10)	{
	// if (rate * totaltime == 0)	{
		if (stateup != statedown)	{
			cerr << "error in MatrixSubstitutionProcess::ResampleAcceptReject: stateup != statedown, efflength == 0\n";
			exit(1);
		}
		delete path;
		path = new BranchSitePath();
		ntrial++;
		path->Reset(stateup);
	}
	else	{
	do	{
		delete path;
		path = new BranchSitePath();
		ntrial++;
		path->Reset(stateup);
		double t = 0;
		int state = stateup;

		if (state != statedown)	{

			// draw waiting time conditional on at least one substitution
			double q = - rate * (*matrix)(state,state);
			double u = -log(1 - rnd::GetRandom().Uniform() * (1 - exp(-q * totaltime))) / q;

			t += u;
			int newstate = matrix->DrawOneStep(state);
			path->Append(newstate,u/totaltime);
			state = newstate;
		}
		while (t < totaltime)	{

			// draw waiting time
			double q = - rate * (*matrix)(state,state);
			double u = -log (1 - rnd::GetRandom().Uniform()) / q;

			t += u;
			if (t < totaltime)	{
				int newstate = matrix->DrawOneStep(state);
				path->Append(newstate,u/totaltime);
				state = newstate;
			}
			else	{
				t -= u;
				u = totaltime - t;
				path->last->SetRelativeTime(u/totaltime);
				t = totaltime;
			}
		}
	} while ((ntrial < maxtrial) && (path->last->GetState() != statedown));
	}

	// if endstate does not match state at the corresponding end of the branch
	// just force it to match
	// however, this is really dirty !
	// normally, in that case, one should give up with accept-reject
	// and use a uniformized method instead (but not yet adapted to the present code, see below)
	if (path->last->GetState() != statedown)	{
		// fossil
		// path->last->SetState(statedown);
		delete path;
		path = 0;
	}

	return path;
}

BranchSitePath* MatrixSubstitutionProcess::ResampleUniformized(int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix)	{

	double length = rate * totaltime;
	int m = matrix->DrawUniformizedSubstitutionNumber(stateup, statedown, length);

	vector<double> y(m+1);
	for (int r=0; r<m; r++)	{
		y[r] = rnd::GetRandom().Uniform();
	}
	y[m] = 1;
	sort(y.begin(),y.end());

	int state = stateup;

	BranchSitePath* path = new BranchSitePath();
	path->Reset(stateup);

	double t = y[0];
	for (int r=0; r<m; r++)	{
		int k = (r== m-1) ? statedown : matrix->DrawUniformizedTransition(state,statedown,m-r-1);
		if (k != state)	{
			path->Append(k,t);
			t = 0;
		}
		state = k;
		t += y[r+1] - y[r];
	}
	path->last->SetRelativeTime(t);
	return path;
}
