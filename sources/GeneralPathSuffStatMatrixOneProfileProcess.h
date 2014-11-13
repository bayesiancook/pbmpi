
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSMATONEPROFILE_H
#define GENPATHSSMATONEPROFILE_H

#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "MatrixOneProfileProcess.h"

// general subclass for all Matrix One mixtures using generic sufficient statistics
class GeneralPathSuffStatMatrixOneProfileProcess : public virtual MatrixOneProfileProcess, public virtual GeneralPathSuffStatMatrixProfileProcess  {

	public:

	GeneralPathSuffStatMatrixOneProfileProcess() {}
	virtual ~GeneralPathSuffStatMatrixOneProfileProcess() {}

	protected:

	// collects site-specific suffstats and pools them componentwise
	void UpdateProfileSuffStat();

	double ProfileSuffStatLogProb();

	map<int,int> profilerootcount;
	map< pair<int,int>, int> profilepaircount;
	map<int,double> profilewaitingtime;

};

#endif

