
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef CONJUGATEPATH_H
#define CONJUGATEPATH_H

#include "RandomBranchSitePath.h"
#include "Conjugate.h"
#include "RandomSubMatrix.h"
#include "PhyloProcess.h"

#include <utility>
#include <map>

#include "CodonSubMatrix.h"

class PathConjugate : public DSemiConjugatePrior<void>	{

	public:

	PathConjugate(Var<PosReal>* inlength, Var<PosReal>* inrate, RandomSubMatrix* inmatrix, Var<Profile>* instationary) {
		// myprocess = inprocess;
		length = inlength;
		rate = inrate;
		matrix = inmatrix;
		stationary = instationary;
		if ((! matrix) && (! stationary))	{
			cerr << "error in RandomBranchSitePath: should specify a matrix or a stationary\n";
			exit(1);
		}

		Register(length);
		Register(rate);
		Register(matrix);
		Register(stationary);

		CreateSuffStat();
	}

	~PathConjugate() {
		DeleteSuffStat();
	}

	void specialUpdate()	{
	}

	// accessors?

	// PhyloProcess* GetPhyloProcess() {return myprocess;}
	Var<PosReal>* GetRandomRate() {return rate;}
	Var<PosReal>* GetRandomTime() {return length;}
	RandomSubMatrix* GetRandomSubMatrix() {return matrix;}
	Var<Profile>* GetRandomStationary() {return stationary;}

	virtual SubMatrix* 		GetMatrix()		{return matrix;}
	virtual double 			GetRate() 		{return rate ?((double) rate->val()) : 1;}	
	virtual const double* 		GetStationary()		{return stationary ? stationary->GetArray() : matrix->GetStationary();}
	virtual double 			GetTime()		{return length ? ((double) length->val()) : 0;}

	int GetNstate() {return matrix->GetNstate();}
	// StateSpace* GetStateSpace() {return matrix->GetStateSpace();}

	// suff stats:
	// for each state: total number of root occurrences
	// for each state: total waiting time
	// for each pair of state : total number of transitions

	void CreateSuffStat()	{
	}

	void DeleteSuffStat()	{
	}

	void ResetSufficientStatistic()	{
		rootcount.clear();
		paircount.clear();
		waitingtime.clear();
	}

	void SaveSufficientStatistic()	{
		cerr << "save sufficient\n";
		exit(1);
	}

	void RestoreSufficientStatistic()	{
		cerr << "restore sufficient\n";
		exit(1);
	}
	
	void IncrementRootCount(int state)	{
		rootcount[state]++;
	}

	void AddRootCount(int state, int count)	{
		rootcount[state] += count;
	}

	void AddWaitingTime(int state, double time)	{
		waitingtime[state] += time;
	}

	void IncrementPairCount(int state1, int state2)	{
		paircount[ pair<int,int>(state1,state2)]++;
	}

	void AddPairCount(int state1, int state2, int count)	{
		paircount[ pair<int,int>(state1,state2)] += count;
	}

	int GetPairCount(int state1, int state2)	{
		return  paircount[ pair<int,int>(state1,state2)];
	}

	void PrintSuffStat(ostream& os)	{

		// if is root
		if (! GetTime())	{
			os << rootcount.size() << '\t';
			for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
				os << i->first << '\t' << i->second << '\t';
			}
			os << '\t';
		}
		else	{
			os << waitingtime.size() << '\t';
			for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
				os << i->first << '\t' << i->second << '\t';
			}
			os << '\t';
			os << paircount.size() << '\t';
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				os << i->first.first << '\t' << i->first.second << '\t' << i->second << '\t';
			}
			os << '\t';
		}
	}
				
	void AddSuffStatFromStream(istream& is)	{

		// if is root
		if (! GetTime())	{
			int n;
			is >> n;
			for (int i=0; i<n; i++)	{
				int f, s;
				is >> f >> s;
				AddRootCount(f,s);
			}
		}
		else	{
			int n;
			is >> n;
			for (int i=0; i<n; i++)	{
				int f;
				double s;
				is >> f >> s;
				AddWaitingTime(f,s);
			}
			is >> n;
			for (int i=0; i<n; i++)	{
				int f,s,t;
				is >> f >> s >> t;
				AddPairCount(f,s,t);
			}
		}
	}
				
	/*
	void PrintSuffStat(ostream& os)	{
		for (int state=0; state<GetMatrix()->GetNstate(); state++)	{
			os << rootcount[state] << '\t';
		}
		for (int state=0; state<GetMatrix()->GetNstate(); state++)	{
			os << waitingtime[state] << '\t';
		}
		for (int state1=0; state1<GetMatrix()->GetNstate(); state1++)	{
			for (int state2=0; state2<GetMatrix()->GetNstate(); state2++)	{
				if (state1 != state2)	{
					os << GetPairCount(state1,state2) << '\t';
				}
			}
		}
		os << '\t';
	}
	*/

	int GetSynTotal(int nuc1, int nuc2)	{
		CodonSubMatrix* codmat = dynamic_cast<CodonSubMatrix*>(GetMatrix());
		if (! codmat)	{
			cerr << "error in conjugate path: not codon\n";
			cerr << GetMatrix() << '\n';
			cerr << GetMatrix()->GetNstate() << '\n';
			exit(1);
		}
		CodonStateSpace* statespace = codmat->GetCodonStateSpace();
		int total = 0;
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			int state1 = i->first.first;
			int state2 = i->first.second;
			int n = i->second;
			bool incl = true;
			if ((nuc1 != -1) && (nuc2 != -1))	{
				int pos = statespace->GetDifferingPosition(state1,state2);
				int n1 = statespace->GetCodonPosition(pos,state1);
				int n2 = statespace->GetCodonPosition(pos,state2);
				if (((nuc1 == n1) && (nuc2 == n2)) || ((nuc1 == n2) && (nuc2 == n1)))	{
					incl = true;
				}
				else {
					incl = false;
				}
			}
			if (incl && statespace->Synonymous(state1,state2))	{
				total += n;
			}
		}
		return total;
	}

	int GetNonSynTotal(int nuc1, int nuc2)	{
		CodonSubMatrix* codmat = dynamic_cast<CodonSubMatrix*>(GetMatrix());
		if (! codmat)	{
			cerr << "error in conjugate path: not codon\n";
			cerr << GetMatrix() << '\n';
			cerr << GetMatrix()->GetNstate() << '\n';
			exit(1);
		}
		CodonStateSpace* statespace = codmat->GetCodonStateSpace();
		int total = 0;
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			int state1 = i->first.first;
			int state2 = i->first.second;
			int n = i->second;
			bool incl = true;
			if ((nuc1 != -1) && (nuc2 != -1))	{
				int pos = statespace->GetDifferingPosition(state1,state2);
				int n1 = statespace->GetCodonPosition(pos,state1);
				int n2 = statespace->GetCodonPosition(pos,state2);
				if (((nuc1 == n1) && (nuc2 == n2)) || ((nuc1 == n2) && (nuc2 == n1)))	{
					incl = true;
				}
				else {
					incl = false;
				}
			}
			if (incl && (! statespace->Synonymous(state1,state2)))	{
				total += n;
			}
		}
		return total;
	}

	double SuffStatLogProb()	{
		double total = 0;

		// root part
		if (! GetTime())	{
			const double* rootstat = GetStationary();
			for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
				total += i->second * log(rootstat[i->first]);
			}
		}

		// non root part
		if (GetTime())	{
			int totnsub = 0;
			double totscalestat = 0;
			SubMatrix& mat = *GetMatrix();
			for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
				totscalestat += i->second * mat(i->first,i->first);
			}
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				total += i->second * log(mat(i->first.first, i->first.second));
				totnsub += i->second;
			}
			total += GetRate() * GetTime() * totscalestat;
			if (totnsub)	{
				total += totnsub * log(GetRate() * GetTime());
			}
		}
		// return 0;
		return total;
	}

	double logProb()	{
		if (isActive())	{
			return SuffStatLogProb();
		}
		return 0;
	}

	protected:

	map<int,int> rootcount;
	map< pair<int,int>, int> paircount;
	map<int,double> waitingtime;
	
	Var<PosReal>* rate;
	Var<PosReal>* length;
	RandomSubMatrix* matrix;
	Var<Profile>* stationary;
};


class ConjugateRandomBranchSitePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath	{

	public:

	ConjugateRandomBranchSitePath(PhyloProcess* inprocess, PathConjugate* inpathconj) : RandomBranchSitePath(inprocess)	{

		pathconj = inpathconj;

		if (!pathconj)	{
			cerr << "error in ConjugateRandomBranchSitePath\n";
			exit(1);
		}

		length = pathconj->GetRandomTime();
		rate = pathconj->GetRandomRate();
		matrix = pathconj->GetRandomSubMatrix();
		stationary = pathconj->GetRandomStationary();
		active_flag = true;

		if ((! matrix) && (! stationary))	{
			cerr << "error in RandomBranchSitePath: should specify a matrix or a stationary\n";
			exit(1);
		}

		conjugate_up.insert(pathconj);
		Register(pathconj); 
	}

	protected:

	void AddSufficientStatistic(SemiConjPrior* parent) {
		if (parent != pathconj)	{
			cerr << "error in ConjugateRandomBranchSitePath::AddSufficientStatistic\n";
			exit(1);
		}

		if (isRoot())	{
			pathconj->IncrementRootCount(Init()->GetState());
		}
		else	{
			Plink* link = Init();
			while (link)	{
				int state = link->GetState();
				pathconj->AddWaitingTime(state,GetRelativeTime(link));
				if (link != last)	{
				int newstate = link->Next()->GetState();
					pathconj->IncrementPairCount(state,newstate);
				}
				link = link->Next();
			}
		}
	}

	private:

	PathConjugate* pathconj;
};


class PathConjugateTree : public BranchValPtrTree<PathConjugate>	{

	public:

	PathConjugateTree(LengthTree* inlengthtree, SequenceAlignment* indata)	{
		lengthtree = inlengthtree;
		data = indata;
	}

	~PathConjugateTree()	{
	}

	LengthTree* GetLengthTree() {return lengthtree;}
	SequenceAlignment* GetData() {return data;}

	Tree* GetTree() {return lengthtree->GetTree();}

	void ActivateSufficientStatistic()	{
		ActivateSufficientStatistic(GetRoot());
	}

	void InactivateSufficientStatistic()	{
		InactivateSufficientStatistic(GetRoot());
	}

	void PrintSynNonSyn(ostream& os, int nuc1 = -1, int nuc2 = -1)	{
		RecursivePrintSynNonSyn(os,GetRoot(), nuc1, nuc2);
		os << '\n';
	}

	void PrintSuffStat(ostream& os)	{
		RecursivePrintSuffStat(os,GetRoot());
		os << '\n';
	}

	double GetLogProb()	{
		return RecursiveGetLogProb(GetRoot());
	}

	void ReadFromFile(string infile)	{
		ifstream is(infile.c_str());
		ActivateSufficientStatistic();
		int N;
		is >> N;
		for (int i=0; i<N; i++)	{
			RecursiveAddSuffStatFromStream(is,GetRoot());
		}
		/*
		ofstream os("check");
		PrintSuffStat(os);
		*/
	}

	protected:

	double RecursiveGetLogProb(const Link* from)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLogProb(link->Out());
		}
		total += GetBranchVal(from->GetBranch())->logProb();
		return total;
	}

	void RecursiveAddSuffStatFromStream(istream& is, const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddSuffStatFromStream(is,link->Out());
		}
		GetBranchVal(from->GetBranch())->AddSuffStatFromStream(is);
	}

	void RecursivePrintSuffStat(ostream& os, const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivePrintSuffStat(os,link->Out());
		}
		GetBranchVal(from->GetBranch())->PrintSuffStat(os);
	}

	void RecursivePrintSynNonSyn(ostream& os, const Link* from, int nuc1, int nuc2)	{

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivePrintSynNonSyn(os,link->Out(), nuc1, nuc2);
		}
		if (! from->isRoot())	{
			os << GetBranchVal(from->GetBranch())->GetNonSynTotal(nuc1, nuc2) << '\t' << GetBranchVal(from->GetBranch())->GetSynTotal(nuc1, nuc2) << '\t';
		}
	}

	PathConjugate* GetPathConjugate(const Branch* branch)	{
		return dynamic_cast<PathConjugate*>(GetBranchVal(branch));
	}

	void ActivateSufficientStatistic(const Link* from)	{
		GetPathConjugate(from->GetBranch())->ActivateSufficientStatistic();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ActivateSufficientStatistic(link->Out());
		}
	}

	void InactivateSufficientStatistic(const Link* from)	{
		GetPathConjugate(from->GetBranch())->InactivateSufficientStatistic();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			InactivateSufficientStatistic(link->Out());
		}
	}

	LengthTree* lengthtree;
	SequenceAlignment* data;

};

class BranchMatrixPathConjugateTree : public PathConjugateTree	{

	public:

	BranchMatrixPathConjugateTree(LengthTree* inlengthtree, BranchValPtrTree<RandomSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PathConjugateTree(inlengthtree, indata) {
		matrixtree = inmatrixtree;
		SetWithRoot(true);
		if (! matrixtree->WithRoot())	{
			cerr << "error in PathConjugateTree: matrixtree does not have root value\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~BranchMatrixPathConjugateTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	PathConjugate* CreateBranchVal(const Link* link)	{
		return new PathConjugate(lengthtree->GetBranchVal(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}
		
	BranchValPtrTree<RandomSubMatrix>* matrixtree;

};

class OneMatrixPathConjugateTree : public PathConjugateTree	{

	public:

	OneMatrixPathConjugateTree(LengthTree* inlengthtree, RandomSubMatrix* inmatrix,  SequenceAlignment* indata) : PathConjugateTree(inlengthtree, indata) {
		matrix = inmatrix;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~OneMatrixPathConjugateTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	PathConjugate* CreateBranchVal(const Link* link)	{
		return new PathConjugate(lengthtree->GetBranchVal(link->GetBranch()), 0, matrix, 0);
	}
		
	RandomSubMatrix* matrix;

};


class PathConjugatePhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	PathConjugatePhyloProcess(PathConjugateTree* inpathconjtree) : PhyloProcess(inpathconjtree->GetLengthTree(), inpathconjtree->GetData())	{
		pathconjtree = inpathconjtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		if (! pathconjtree->GetBranchVal(link->GetBranch()))	{
			cerr << "error in create branch val\n";
			exit(1);
		}
		return  new ConjugateRandomBranchSitePath(this,pathconjtree->GetBranchVal(link->GetBranch()));
	}

	protected:
	PathConjugateTree* pathconjtree;

};

#endif

