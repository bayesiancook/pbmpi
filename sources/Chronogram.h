
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef CHRONOGRAM_H
#define CHRONOGRAM_H

#include "ValTree.h"
#include "Tree.h"
#include <math.h>
#include <algorithm>
#include "RandomTypes.h"

template<class U, class V> class NodeBranchProcess : public MCMC , public NodeBranchValPtrTree<Rvar <U>, Dvar <V> > {

	public:
	U val(const Node* node)	{
		return this->GetNodeVal(node)->val();
	}

	void setval(const Node* node, U inval)	{
		this->GetNodeVal(node)->setval(inval);
	}

	V val(const Branch* branch)	{
		return this->GetBranchVal(branch)->val();
	}

	void setval(const Branch* branch, V inval)	{
		this->GetBranchVal(branch)->setval(inval);
	}

	virtual void drawSample()	{
		// drawSample(this->GetRoot());
	}

	virtual double GetLogProb(const Link* from)	{
		double total = this->GetNodeVal(from->GetNode())->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetLogProb(link->Out());
		}
		return total;
	}

	virtual double GetLogProb()	{
		return GetLogProb(this->GetRoot());
	}

};

class NodeDate : public Rvar<PosReal>{
	
	protected:
	
	double min, max;

	public:

	NodeDate(Var<PosReal>* inRate){
		Register(inRate);
		SetName("date");
	}

	void drawSample(){}
	

	virtual void SetMinMax(double inMin, double inMax){
		min = inMin;
		max = inMax;
		if (max < min)	{
			cerr << "error : max < min : " << min << '\t' << max << '\n';
			exit(1);
		}
	}

	double ProposeMove(double tuning)	{
		double m = (max - min) * tuning * (rnd::GetRandom().Uniform() - 0.5);
		value += m;
		while ((value < min) || (value > max)){
			if (value < min){
				value = 2 * min - value;
			}
			if (value > max){
				value = 2 * max - value;
			}
		}
		return 0;
	}

	void EmptyMove()	{
		value = (max -  min) * rnd::GetRandom().Uniform() + min;
	}

	virtual double logProb(){return 0;}

};


class BranchLength : public Dvar<PosReal>{

	private:
	Var<PosReal>* Before;
	Var<PosReal>* After;
	Var<PosReal>* Rate;
	

	public:
	BranchLength(Var<PosReal>* inBefore, Var<PosReal>* inAfter, Var<PosReal>* inRate){
		Before = inBefore;
		After = inAfter;
		Rate = inRate;
		Register(inRate);
		Register(inBefore);
		Register(inAfter);
		SetName("time");
	}


	~BranchLength(){};


	void specialUpdate(){
		double newval = (Before->val() - After->val()) * Rate->val();
		setval(newval);
	}
};



class Chronogram : public NodeBranchProcess<PosReal, PosReal> {

	public:

	Chronogram() {}

	Chronogram(Tree* intree, Var<PosReal>* inRate) {
		SetWithRoot(false);
		tree = intree;
		rate = inRate;
		// Create objects
		totnode = 0;
		RecursiveCreateNode(GetRoot());
		cerr << "total number of nodes : " << totnode << '\n';
		RecursiveCreateBranch(GetRoot());
		// Set values and make the tree ultra-metric
		double maxage = RecursiveSetNodesValues(GetRoot());
		RecursiveEqualizeLeafNodes(GetRoot(),maxage);
		RecursiveNormalizeTree(GetRoot(),maxage,true);
		RecursiveUpdateBranches(GetRoot());
		/*
		cerr << GetMinDeltaTime(GetRoot()) << '\n';
		while (GetMinDeltaTime(GetRoot()) < 2e-3)	{
			EmptyMove();
			cerr << GetMinDeltaTime(GetRoot()) << '\n';
		}
		*/
	}

	Chronogram(Tree* intree) {
		rate = new Const<PosReal>(1);
		SetWithRoot(false);
		tree = intree;
		// Create objects
		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());
		// Set values and make the tree ultra-metric
		double maxage = RecursiveSetNodesValues(GetRoot());
		RecursiveEqualizeLeafNodes(GetRoot(),maxage);
		RecursiveNormalizeTree(GetRoot(),maxage,true);
		RecursiveUpdateBranches(GetRoot());
		PrintAges(GetRoot(),maxage);
	}

	void PrintAges(const Link* from,double rootage)	{
		if (! from->isLeaf())	{
			cout << *GetNodeDate(from->GetNode()) * rootage << '\n';
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				PrintAges(link->Out(),rootage);
			}
		}
	}

	static bool compposreal(Var<PosReal>* a, Var<PosReal>* b) {return (((double) a->val()) > ((double) b->val()));}

	vector<Var<PosReal>*> GetDateList()	{
		vector<Var<PosReal>*> leaflist = GetLeafList();
		vector<Var<PosReal>*> datelist;
		datelist.push_back(leaflist[0]);
		PushAges(GetRoot(),datelist);
		sort(datelist.begin(), datelist.end(), Chronogram::compposreal);
		if (datelist.size() != GetTree()->GetSize() - 1)	{
			cerr << "error : size of date list : " << datelist.size() << '\t' << GetTree()->GetSize() -1 << '\n';
		}
		return datelist;
	}

	void PushAges(const Link* from, vector<Var<PosReal>*>& datelist)	{
		if (! from->isLeaf())	{
			datelist.push_back(GetNodeDate(from->GetNode()));
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				PushAges(link->Out(), datelist);
			}
		}
	}

	vector<Var<PosReal>*> GetLeafList()	{
		vector<Var<PosReal>*> leaflist;
		PushLeaves(GetRoot(),leaflist);
		if (leaflist.size() != GetTree()->GetSize())	{
			cerr << "error : size of leaf list : " << leaflist.size() << '\t' << GetTree()->GetSize() << '\n';
		}
		return leaflist;
	}

	void PushLeaves(const Link* from, vector<Var<PosReal>*>& leaflist)	{
		if (from->isLeaf())	{
			leaflist.push_back(GetNodeDate(from->GetNode()));
		}
		else	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				PushLeaves(link->Out(), leaflist);
			}
		}
	}

	~Chronogram(){
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree(){
		return tree;
	}

	BranchLength* GetBranchLength(const Branch* branch)	{
		BranchLength* tmp = dynamic_cast<BranchLength*>(GetBranchVal(branch));
		if (!tmp)	{
			cerr << "error in Chronogram::GetBranchLength : dynamic cast \n";
			cerr << "branch : " << branch << '\n';
			exit(1);
		}
		return tmp;
	}


	NodeDate* GetNodeDate(const Node* node)	{
		NodeDate* tmp = dynamic_cast<NodeDate*>(GetNodeVal(node));
		if (! tmp)	{
			cerr << "error in Chronogram::GetNodeDate : dynamic cast\n";
			exit(1);
		}
		return tmp;
	}

	Var<PosReal>* GetRootAge()	{return GetNodeVal(GetRoot()->GetNode());}

	Var<PosReal>* GetRate() {return rate;}

	void MultiplyTimes(const Link* from, double d)	{
		RecursiveMultiplyTimes(from,d);
	}

	double GetAbsoluteTime(const Link* from){
		return GetNodeVal(from->GetNode())->val();
	}

	// called on a node: returns the time span corresponding to the branch arriving at that node
	double GetBranchTimeLength(const Link* from){
		double d = GetAbsoluteTime(from->Out()) - GetAbsoluteTime(from);
		if (d < 0)	{
			d = -d;
			cerr << "error in chronogram: negative time\n";
			cerr << GetAbsoluteTime(from->Out()) << '\t' << GetAbsoluteTime(from) << '\n';
			exit(1);
		}
		return d;
	}

	void Sample()	{
		drawSample(this->GetRoot());
	}
	
	void drawSample(const Link* from){
	}

	void Clamp()	{
		RecursiveClamp(GetRoot());
	}

	void specialUpdate()	{
		RecursiveUpdate(GetRoot());
	}

	double GetTotalTime()	{
		return RecursiveGetTotalTime(GetRoot());
	}

	void RecursiveClamp(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveClamp(link->Out());
		}
		GetNodeVal(from->GetNode())->Clamp();
	}

	double RecursiveGetTotalTime(const Link* from)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetTotalTime(link->Out());
		}
		if (! from->isRoot())	{
			total += GetBranchVal(from->GetBranch())->val();
		}
		return total;
	}

	void RecursiveUpdate(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveUpdate(link->Out());
		}
		if (! from->isRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
	}


	double Move(double tuning){
		int n = 0;
		double tmp = GeneralMoveTime(tuning, GetRoot(),n);
		return tmp / n;
	}

	void EmptyMove(){
		EmptyMove(GetRoot());
	}

	void RecursiveRegister(DAGnode* node, const Link* from)	{
		GetNodeVal(from->GetNode())->Register(node);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(node,link->Out());
		}
	}

	double GetMinDeltaTime(const Link* from)	{
		double min = -1;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = GetBranchTimeLength(link->Out());
			if ((min==-1) || (min > tmp))	{
				min = tmp;
			}
			double tmp2 = GetMinDeltaTime(link->Out());
			if ((tmp2 != -1) && (min > tmp2))	{
				min = tmp2;
			}
		}
		return min;
	}

	protected:

	Rvar<PosReal>* CreateNodeVal (const Link* link){
		totnode++;
		return new NodeDate(rate);
	}

	void RecursiveCreateNode(Link* from) {
		nodeval[from->GetNode()] = CreateNodeVal(from);
		for(Link* link=from->Next(); link!=from; link=link->Next()) {
			RecursiveCreateNode(link->Out());
		}	
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link){
		BranchLength* bl = new BranchLength(GetNodeVal(link->GetNode()),GetNodeVal(link->Out()->GetNode()),rate);
		return bl;
	}

	void RecursiveCreateBranch(Link* from) {
		for(Link* link=from->Next(); link!=from; link=link->Next()) {
			branchval[link->GetBranch()] = CreateBranchVal(link);
			RecursiveCreateBranch(link->Out());
		}
	}

	void RecursiveMultiplyTimes(const Link* from, double d)	{
		if (! from->isLeaf())	{
			GetNodeVal(from->GetNode())->setval(d * GetNodeVal(from->GetNode())->val());
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next()) {
			RecursiveMultiplyTimes(link->Out(),d);
		}
	}

	// Set nodes's absolute times based on observed branch lengths*/
	// root is 0 
	// returns max depth
	double RecursiveSetNodesValues(const Link* from){
		double maxage = 0;
		if (from->isRoot()){
			GetNodeVal(from->GetNode())->setval(0);
		}
		else{
			double l = atof(from->GetBranch()->GetName().c_str());
			if (l <= 0)	{
				cerr << "warning : non strictly positive branch length : " << l << '\n';
				cerr << "correcting and setting to 0.001\n";
				l = 0.001;
			}
			double val = l + GetAbsoluteTime(from->Out());
			GetNodeVal(from->GetNode())->setval(val);
			maxage = val;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next()) {
			double tmp = RecursiveSetNodesValues(link->Out());
			if (maxage < tmp)	{
				maxage = tmp;
			}
		}
		return maxage;
	}

	void RecursiveEqualizeLeafNodes(const Link* from, double maxage)	{
		if (from->isLeaf())	{
			if (GetNodeVal(from->GetNode())->val() > maxage)	{
				cerr << "error in leaf age equalization : " << GetNodeVal(from->GetNode())->val() << '\t' << maxage << '\n';
				exit(1);
			}
			GetNodeVal(from->GetNode())->setval(maxage);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next()) {
			RecursiveEqualizeLeafNodes(link->Out(),maxage);
		}
	}

	void RecursiveNormalizeTree(const Link* from, double maxage, bool invert)	{
		if (invert)	{
			GetNodeVal(from->GetNode())->setval(1 - GetNodeVal(from->GetNode())->val() / maxage);
		}
		else	{
			GetNodeVal(from->GetNode())->setval(GetNodeVal(from->GetNode())->val() / maxage);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next()) {
			RecursiveNormalizeTree(link->Out(),maxage, invert);
		}
	}

	void RecursiveUpdateBranches(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next()){
			GetBranchLength(link->GetBranch())->specialUpdate();
			RecursiveUpdateBranches(link->Out());
		}
	}

	void SetMinMax(Link* from){
		double max = GetAbsoluteTime(from->Out());
		double min = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next()){
			if(GetAbsoluteTime(link->Out()) > min){
				min = GetAbsoluteTime(link->Out());
			}
		}
		GetNodeDate(from->GetNode())->SetMinMax(min, max);
	}

	double MoveTime(double tuning, Link* from){
		SetMinMax(from);
		return GetNodeVal(from->GetNode())->Move(tuning);
	}


	double GeneralMoveTime(double tuning, Link* from, int& n){
		if(from->isLeaf()){
			return 0;
		}
		else{
			double retour = 0;
			for(Link* link=from->Next(); link!=from; link=link->Next()){
				retour += GeneralMoveTime(tuning,link->Out(),n);
			}
			if(!from->isRoot()){
				retour += MoveTime(tuning, from);
				n++;
			}
			return retour;
		}
	}

	void EmptyMove(Link* from){
		if(! from->isLeaf()){
			if(!from->isRoot()){
				SetMinMax(from);
				GetNodeDate(from->GetNode())->EmptyMove();
			}
			for(Link* link=from->Next(); link!=from; link=link->Next()){
				EmptyMove(link->Out());
			}
			if(!from->isRoot()){
				SetMinMax(from);
				GetNodeDate(from->GetNode())->EmptyMove();
			}
		}
	}


	protected:

	Tree* tree;
	Var<PosReal>* rate;

	int totnode;
	 
};

#endif
