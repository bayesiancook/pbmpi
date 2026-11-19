
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXSUB_H
#define MATRIXSUB_H

#include "SubstitutionProcess.h"
#include "MatrixProfileProcess.h"
#include <vector>
#include <unordered_map>
#include <functional>

class MatrixSubstitutionProcess : public virtual SubstitutionProcess, public virtual MatrixProfileProcess	{

	public:

	MatrixSubstitutionProcess() {}
	virtual ~MatrixSubstitutionProcess() {}

	virtual int GetNstate(int site) {return GetMatrix(site)->GetNstate();}
	// virtual int GetNstate() {return GetMatrix(0)->GetNstate();}

	virtual const double* GetStationary(int site)	{
		return GetMatrix(site)->GetStationary();
	}

	protected:

	// CPU Level 3: implementations of likelihood propagation and substitution mapping methods
	void Propagate(double*** from, double*** to, double time, bool condalloc = false);
	void PropagateFromLeaf(const int* state, double*** from, double*** to, double time, bool condalloc = false);

	// common implementation; leafstate == 0 : general case
	void PropagateImpl(double*** from, double*** to, double time, bool condalloc, const int* leafstate);

	private:

	// scratch space for PropagateImpl, kept between calls to avoid reallocations:
	// values of exp(time * rate * eigenvalue), computed once per (matrix, rate) pair and per call
	struct ExpKey	{
		const SubMatrix* matrix;
		double rate;
		bool operator==(const ExpKey& o) const {return (matrix == o.matrix) && (rate == o.rate);}
	};
	struct ExpKeyHash	{
		size_t operator()(const ExpKey& k) const {
			return std::hash<const void*>()(k.matrix) ^ (std::hash<double>()(k.rate) * 0x9e3779b97f4a7c15ULL);
		}
	};
	std::unordered_map<ExpKey,size_t,ExpKeyHash> expindex;
	std::vector<double> expbuffer;
	std::vector<double> propaux;

	protected:

	void SitePropagate(int site, double** from, double** to, double time, bool condalloc = false);

	BranchSitePath** SamplePaths(int* stateup, int* statedown, double time);
	BranchSitePath** SampleRootPaths(int* rootstate);
	BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix);
	BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix);

	void SimuPropagate(int* stateup, int* statedown, double time);
};

#endif

