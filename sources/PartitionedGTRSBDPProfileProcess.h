
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRSBDPPROFILE_H
#define PARTGTRSBDPPROFILE_H

#include "PartitionedGTRMixtureProfileProcess.h"
#include "SBDPProfileProcess.h"

class PartitionedGTRSBDPProfileProcess : public virtual PartitionedGTRMixtureProfileProcess, public virtual SBDPProfileProcess {

	public:

	PartitionedGTRSBDPProfileProcess() {}
	virtual ~PartitionedGTRSBDPProfileProcess() {}

	protected:

	virtual void SwapComponents(int cat1, int cat2);

	virtual void Create(int indim, PartitionScheme rrscheme)	{
		PartitionedGTRMixtureProfileProcess::Create(indim, rrscheme);
		SBDPProfileProcess::Create(rrscheme.GetNsite(),indim);
	}

	virtual void Delete()	{
		SBDPProfileProcess::Delete();
		PartitionedGTRMixtureProfileProcess::Delete();
	}

	double MixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual double GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual void SlaveMixMove();

	double GlobalIncrementalDPMove(int nrep, double epsilon);
	void SlaveIncrementalDPMove();
	double IncrementalDPMove(int nrep, double epsilon);

};

#endif
