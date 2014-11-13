
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSMATPROFILE_H
#define GENPATHSSMATPROFILE_H

#include "MatrixProfileProcess.h"

// superclass for all matrix implementations using generic sufficient statistics
// generic sufficient statistics are: total time in each state, number of transitions between each pair of states, number of times in each state at the root
class GeneralPathSuffStatMatrixProfileProcess : public virtual MatrixProfileProcess	{

	public:

	GeneralPathSuffStatMatrixProfileProcess() {}
	virtual ~GeneralPathSuffStatMatrixProfileProcess() {}

	// will be implemented in phyloprocess
	// return the sufficient statistics for a given site
	virtual map<pair<int,int>,int>& GetSitePairCount(int site) = 0;
	virtual map<int,double>& GetSiteWaitingTime(int site) = 0;
	virtual int GetSiteRootState(int site) = 0;

};

#endif

