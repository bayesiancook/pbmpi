
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef BRANCH_H
#define BRANCH_H

#include "Tree.h"
#include "Chrono.h"

class BranchProcess : public NewickTree {

	public:

	BranchProcess() : tree(0), blarray(0), lengthfrac(1.0) {}
	virtual ~BranchProcess() {}

	virtual string GetVersion() = 0;

	// 
	NewickTree* GetLengthTree() {return this;}
	Tree* GetTree() {return tree;}
	Link* GetRoot() {return tree->GetRoot();}
	const Link* GetRoot() const {return tree->GetRoot();}

	string GetBranchName(const Link* link) const;
	string GetNodeName(const Link* link) const;

	double GetLength(const Branch* branch) {return branch ? blarray[branch->GetIndex()] : 0;}
	double GetLength(const Branch* branch) const {return branch ? blarray[branch->GetIndex()] : 0;}
	void SetLength(const Branch* branch, double inlength) {
		if (! branch)	{
			cerr << "error in branch process: null branch\n";
			exit(1);
		}
		if (! blarray)	{
			cerr << "array not created\n";
			exit(1);
		}
		blarray[branch->GetIndex()] = inlength;
	}

	void Backup();
	void Restore();

	double ProposeMove(const Branch* branch, double tuning);
	void MoveBranch(const Branch* branch, double factor);
	void Restore(const Branch* branch);

	virtual double GetNormFactor()  = 0;

	// return number of branches
	int MoveAllBranches(double factor);
	int RecursiveMoveAllBranches(const Link* from, double e);

	virtual void SampleLength(const Branch* branch) = 0;
	virtual double LogBranchLengthPrior(const Branch* branch) = 0;

	double GetTotalLength()	{
		return RecursiveTotalLength(GetRoot());
	}

	double GetRenormTotalLength()	{
		return RecursiveTotalLength(GetRoot()) * GetNormFactor();
	}

	void RenormalizeBranchLengths()	{
		double tmp = GetNormFactor();
		RecursiveNormalizeBranchLengths(GetRoot(),tmp);
	}

	void DenormalizeBranchLengths()	{
		double tmp = GetNormFactor();
		RecursiveNormalizeBranchLengths(GetRoot(),1.0 / tmp);
	}

	void RecursiveNormalizeBranchLengths(const Link* from, double factor);
	
	double RecursiveTotalLength(const Link* from);

	void SetLengthsFromNames()	{
		if (blarray)	{
			RecursiveSetLengthsFromNames(GetRoot());
		}
		else	{
			cerr << "set lengths from names called without blarray\n";
		}
	}
	
	void SetNamesFromLengths()	{
		if (blarray)	{
			RecursiveSetNamesFromLengths(GetRoot());
		}
		else	{
			cerr << "set names from length called without blarray\n";
		}
	}

	void RecursiveSetLengthsFromNames(const Link* from);
	void RecursiveSetNamesFromLengths(const Link* from);

	void SampleLength()	{
		RecursiveSampleLength(GetRoot());
	}

    virtual void PriorSampleLength() = 0;

    double LogBranchPrior() {
        return LogHyperPrior() + LogLengthPrior();
    }

    virtual double LogHyperPrior() = 0;

	double LogLengthPrior()	{
		return RecursiveLogLengthPrior(GetRoot());
	}

	virtual void GlobalUpdateBranchLengthSuffStat() = 0;
	virtual void SlaveUpdateBranchLengthSuffStat() = 0;
	virtual void UpdateBranchLengthSuffStat() = 0;

	/*
	virtual double GetBranchLengthSuffStatBeta(const Branch* branch) = 0;
	virtual int GetBranchLengthSuffStatCount(const Branch* branch) = 0;
	*/

	virtual double GetBranchLengthSuffStatBeta(int index) = 0;
	virtual int GetBranchLengthSuffStatCount(int index) = 0;

	// implements a map<const Branch*, double>

	// Move function ?
	// how about the tuning parameters ?

	double LengthSuffStatLogProb();

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;

	virtual const TaxonSet* GetTaxonSet() const = 0;

	protected:

	int GetNbranch()	{
		return tree->GetNbranch();
	}

	int GetNnode()	{
		return tree->GetNnode();
	}

	int GetNlink()	{
		return tree->GetNlink();
	}

	virtual void Create(Tree* intree) {
		tree = intree;
		blarray = new double[GetNbranch()];
		bkarray = new double[GetNbranch()];
		blarray[0] = 0;
		SetLengthsFromNames();
	}
	virtual void Delete() {
		delete[] blarray;
		delete[] bkarray;
	}

	double RecursiveLogLengthPrior(const Link* from);
	void RecursiveSampleLength(const Link* from);

    void SetLengthFrac(double infrac)   {
        lengthfrac = infrac;
    }

	Tree* tree;
	double* blarray;
	double* bkarray;

    double lengthfrac;

	Chrono chronolength;

};

#endif

