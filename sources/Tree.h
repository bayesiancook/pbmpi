
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef TREE_H
#define TREE_H

#include <cstdlib>
#include <string>
#include <iostream>
#include "TaxonSet.h"

class Node {

	private:
	int index;
	string name;

	public:

	Node() : index(0) , name("") {}
	Node(string s) : index(0) , name(s) {}
	Node(const Node* from) : index(from->index) , name(from->name) {}

	virtual ~Node() {}
	
	virtual string GetName() const {return name;}
	virtual void SetName(string inname) {name = inname;}
	int GetIndex() const {return index;}
	void SetIndex(int i) {index = i;}

};

class Branch {

	private:
	int index;
	string name;

	public:

	Branch() : index(0) , name("") {}
	Branch(string s) : index(0) , name(s) {}
	Branch(const Branch* from) : index(from->index), name(from->name) {}

	virtual ~Branch() {}
	
	virtual string GetName() const {return name;}
	virtual void SetName(string inname) {name = inname;}
	int GetIndex() const {return index;}
	void SetIndex(int i) {index = i;}
};

class Link	{

	private:
	Link* next;
	Link* out;
	Branch* branch;
	Node* node;
	int index;

	public:

	double* tbl;

	Link()	{
		tbl = 0;
		next = out = this;
		branch = 0;
		node = 0;
	}

	Link(const Link* from)	{
		tbl = 0;
		next = out = this;
		node = 0;
		branch = 0;
	}

	Link* Next() const {return next;}
	Link* Out() const {return out;}
	Branch* GetBranch() const {return branch;}
	Node* GetNode() const {return node;}

	void SetBranch(Branch* inbranch)	{
		/*
		if (branch == inbranch)	{
			cout << "error in Link::SetBranch: branch and new branch are the same\n";
			exit(1);
		}
		delete branch;
		*/
		branch = inbranch;
	}
	void SetNode(Node* innode)	{
		/*
		if (node == innode)	{
			cout << "error in Link::SetNode: node and new node are the same\n";
			exit(1);
		}
		delete node;
		*/
		node = innode;
	}

	string GetNodeLeafSet() const {
		
		string s;
		if (isLeaf()){
			s = GetNode()->GetName()+ '|';
		}
		else{
			for (const Link* link=Next(); link!=this; link=link->Next())	{
				s = s + link->Out()->GetNodeLeafSet();
			}
		}
		return s;
	}

	void SetIndex(int i) {index = i;}

	int GetIndex() const {return index;}

	void SetOut(Link* inout)	{
		out = inout;
	}

	void SetNext(Link* innext)	{
		next = innext;
	}

	void AppendTo(Link* link)	{
		if (link)	{
			link->next = this;
		}
	}

	void Insert(Link* link)	{ // insert link after this
		link->next = next;
		next = link;
	}

	void InsertOut(Link* link)	{ // insert link as out
		link->out = this;
		out = link;
	}


	// Move the next subtree after this link to be the previous of this link.
	// Should be call with 0;
	void Knit()	{
		Link* previous = 0;	
		for (previous=this; previous->next != this; previous=previous->Next());
		previous->next = this->next;
		this->next = next->Next();
		previous->next->next = this;
		
		/*if(!link){
			link = next;
			next=next->next;
			link->next = this;
		}
		if(next == link->next){
			next = link;
		}
		else{
			next->Knit(link);
		}*/
	}

	bool isLeaf() const	{
		return (next == this);
	}

	bool isUnary() const	{
		return (next->Next() == this && !isLeaf());
	}

	bool isRoot() const {
		return (out == this);
	}

	// degree : number of branches connecting to the node associated to this link
	int GetDegree()	const {
		int d = 1;
		const Link* link = next;
		while (link!=this)	{
			d++;
			link=link->next;
		}
		return d;
	}
};


class NewickTree {

	public:

	virtual ~NewickTree() {}
	virtual Link* GetRoot() = 0;
	virtual const Link* GetRoot() const = 0;

	void ToStream(ostream& os) const;
	void ToStream(ostream& os, const Link* from) const;
	double ToStreamSimplified(ostream& os, const Link* from) const;

	string GetLeftMost(const Link* from) const 	{
		if (from->isLeaf())	{
			return GetNodeName(from);
		}
		return GetLeftMost(from->Next()->Out());
	}

	string GetRightMost(const Link* from) const {
		if (from->isLeaf())	{
			return GetNodeName(from);
		}
		const Link* link = from->Next();
		while (link->Next() != from)	{
			link = link->Next();
		}
		return GetRightMost(link->Out());
	}

	static void Simplify()	{
		simplify = true;
	}


	virtual string GetNodeName(const Link* link) const = 0;
	virtual string GetBranchName(const Link* link) const = 0;

	protected:

	static bool simplify;

};


class Tree : public NewickTree {

	public:

	Tree();
	// default constructor: set member pointers to 0

	// void tree;
	Tree(const TaxonSet* intaxset);



	void MakeRandomTree();

	Tree(const Tree* from);
	// clones the entire Link structure
	// but does NOT clone the Nodes and Branches
	// calls RecursiveClone

	Tree(string filename);
	// create a tree by reading into a file (netwick format expected)
	// calls ReadFromStream

	Tree(istream& is);
	// create a tree by reading into a stream (netwick format expected)
	// calls ReadFromStream

	void ReadFromStream(istream& is);
	// reading a tree from a stream:
	// recursively invokes the two following functions

	virtual ~Tree();
	// calls RecursiveDelete
	// but does NOT delete the Nodes and Branches
 
	// Delete the leaf pointing by the next link and set everithing right.
	void DeleteNextLeaf(Link* previous);
	
	Link* Detach(Link* down, Link* up);

	void Attach(Link* down, Link* up, Link* todown, Link* toup);

	void NNIturn(Link* from);

	void RootAtRandom();
	void RootAt(Link* newroot);
	Link* ChooseLinkAtRandom();

	int CountInternalNodes(const Link* from);
	Link* ChooseInternalNode(Link* from, Link*& chosenup, int& n);
	int CountNodes(const Link* from);
	Link* ChooseNode(Link* from, Link*& chosenup, int& n);

	int DrawSubTree(Link*& down, Link*& up);
	void GrepNode(Link* from, Link*& down, Link*& up, int choose);

	bool RecursiveCheckDegree(const Link* from, int test = 3);
	bool CheckRootDegree(int test = 3);

	// Delete the unary Node wich from is paart of and set everithing right.
	void DeleteUnaryNode(Link* from);

	Link* GetRoot() {return root;}
	const Link* GetRoot() const {return root;}
	const TaxonSet* GetTaxonSet() const {return taxset;}

	void MakeTaxonSet()	{
		taxset = new TaxonSet(this);
	}
	
	void RegisterWith(const TaxonSet* taxset,int id = 0);
	// Registers all leaves of the tree with an external TaxonSet
	// the taxon set defines a map between taxon names and indices (between 0 and P-1)
	// the tree is recursively traversed
	// each leaf's name is looked for in the map of the taxon set
	// if not found : an error is emitted
	// otherwise, the leaf's index is set equal to the index of the corresponding taxon
	
	bool RegisterWith(const TaxonSet* taxset, Link* from, int& tot);
	// recursive function called by RegisterWith

	string GetBranchName(const Link* link) const {return link->GetBranch()->GetName();}

	string GetNodeName(const Link* link) const {
		return link->GetNode()->GetName();
		/*
		if (! link->isLeaf())	{
			return link->GetNode()->GetName();
		}
		string s = link->GetNode()->GetName();
		unsigned int l = s.length();
		unsigned int i = 0;
		while ((i < l) && (s[i] != '_')) i++;
		if (i == l)	{
			cerr << "error in get name\n";
			exit(1);
		}
		i++;
		return s.substr(i,l-i);
		*/
	}
	// trivial accessors
	// they can be useful to override, so as to bypass Branch::GetName() and Node::GetName()





	void EraseInternalNodeName();
	void EraseInternalNodeName(Link* from);

	// void Print(ostream& os,const Link* from) const ;
	// void Print(ostream& os) const;
	// printing int netwick format

	unsigned int GetSize()	{
		return GetSize(GetRoot());
	}

	int GetSize(const Link* from)	const {
		if (from->isLeaf())	{
			return 1;
		}
		else	{
			int total = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				total += GetSize(link->Out());
			}
			return total;
		}
		return 0;
	}
				
	int GetFullSize(const Link* from)	const {
		if (from->isLeaf())	{
			return 1;
		}
		else	{
			int total = 1;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				total += GetFullSize(link->Out());
			}
			return total;
		}
		return 0;
	}
				
	const Link* GetLCA(string tax1, string tax2)	{
		bool found1 = false;
		bool found2 = false;
		const Link* link= RecursiveGetLCA(GetRoot(),tax1,tax2,found1,found2);
		// cerr << tax1 << '\t' << tax2 << '\n';
		// Print(cerr,link);
		// cerr << '\n' << '\n';
		return link;
	}

	void Subdivide(Link* from, int Ninterpol);

	string Reduce(Link* from = 0)	{

		if (! from)	{
			from = GetRoot();
		}
		if (from->isLeaf())	{
			cerr << from->GetNode()->GetName() << '\n';;
			return from->GetNode()->GetName();
		}
		else	{
			string name = "None";
			for (Link* link = from->Next(); link!=from; link=link->Next())	{
				string tmp = Reduce(link->Out());
				if (tmp == "diff")	{
					name = "diff";
				}
				else if (name == "None")	{
					name = tmp;
				}
				else if (name != tmp)	{
					name = "diff";
				}
			}
			cerr << '\t' << name << '\n';
			from->GetNode()->SetName(name);
			return name;
		}
		return "";
	}

	void Print(ostream& os, const Link* from = 0)	{
		if (!from)	{
			from = GetRoot();
		}
		if (from->isLeaf())	{
			os << from->GetNode()->GetName();
		}
		else	{
			os << '(';
			for (const Link* link = from->Next(); link!=from; link=link->Next())	{
				Print(os,link->Out());
				if (link->Next() != from)	{
					os << ',';
				}
			}
			os << ')';
		}
		if (from->isRoot())	{
			os << ";\n";
		}
		else{
			os << from->GetBranch()->GetName();
		}
	}


	void PrintReduced(ostream& os, const Link* from = 0)	{
		if (!from)	{
			from = GetRoot();
		}
		if (from->GetNode()->GetName() != "diff")	{
			os << from->GetNode()->GetName();
		}
		else	{
			os << '(';
			for (const Link* link = from->Next(); link!=from; link=link->Next())	{
				PrintReduced(os,link->Out());
				if (link->Next() != from)	{
					os << ',';
				}
			}
			os << ')';
		}
		if (from->isRoot())	{
			os << ";\n";
		}
	}

	void SetIndices()	{
		Nlink = 0;
		Nnode = GetSize();
		Nbranch = 1;
		linkmap.clear();
		nodemap.clear();
		branchmap.clear();
		SetIndices(GetRoot(),Nlink,Nnode,Nbranch);
	}

	int GetNlink()	{
		return Nlink;
	}

	int GetNbranch()	{
		return Nbranch;
	}

	int GetNnode()	{
		return Nnode;
	}

	// maybe

	const Node* GetNode(int index)	{
		return nodemap[index];
	}
	const Branch* GetBranch(int index)	{
		return branchmap[index];
	}

	Link* GetLink(int index)	{
		return linkmap[index];
	}

	protected:

	map<int,const Node*> nodemap;
	map<int,const Branch*> branchmap;
	map<int,Link*> linkmap;

	void CheckIndices(Link* from)	{

		if (! from->isRoot())	{
			if (from->GetBranch() != branchmap[from->GetBranch()->GetIndex()])	{
				cerr << "branch index : " << from->GetBranch()->GetIndex() << '\n';
				exit(1);
			}
		}
		else	{
			if (branchmap[0] != 0)	{
				cerr << "root branch index\n";
				exit(1);
			}
		}

		if (from->GetNode() != nodemap[from->GetNode()->GetIndex()])	{
			cerr << "node index: " << from->GetNode()->GetIndex() << '\n';
			exit(1);
		}

		if (! from->isRoot())	{
			if (from->Out() != linkmap[from->Out()->GetIndex()])	{
				cerr << "link index : " << from->Out()->GetIndex() << '\n';
			}
		}
		if (from != linkmap[from->GetIndex()])	{
			cerr << "link index : " << from->GetIndex() << '\n';
		}


		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			CheckIndices(link->Out());
		}
	}

	void SetIndices(Link* from, int& linkindex, int& nodeindex, int& branchindex)	{

		if (! from->isRoot())	{
			from->GetBranch()->SetIndex(branchindex);
			branchmap[branchindex] = from->GetBranch();
			branchindex++;
		}

		if (! from->isLeaf())	{
			from->GetNode()->SetIndex(nodeindex);
			nodemap[nodeindex] = from->GetNode();
			nodeindex++;
		}
		else	{
			nodemap[from->GetNode()->GetIndex()] = from->GetNode();
		}

		if (! from->isRoot())	{
			from->Out()->SetIndex(linkindex);
			linkmap[linkindex] = from->Out();
			linkindex++;
		}
		from->SetIndex(linkindex);
		linkmap[linkindex] = from;
		linkindex++;

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			SetIndices(link->Out(),linkindex,nodeindex,branchindex);
		}
	}

	// returns 0 if not found
	// returns link if found (then found1 and found2 must 
	const Link* RecursiveGetLCA(const Link* from, string tax1, string tax2, bool& found1, bool& found2)	{
		const Link* ret= 0;
		if (from->isLeaf())	{
			found1 |= (from->GetNode()->GetName() == tax1);
			found2 |= (from->GetNode()->GetName() == tax2);
			if (! ret)	{
				if (found1 && found2)	{
					ret = from;
				}
			}
		}
		else	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				bool tmp1 = false;
				bool tmp2 = false;
				const Link* ret2 = RecursiveGetLCA(link->Out(),tax1,tax2,tmp1,tmp2);
				found1 |= tmp1;
				found2 |= tmp2;
				if (ret2)	{
					if (ret)	{
						cerr << "error : found node twice\n";
						cerr << tax1 << '\t' << tax2 << '\n';
						ToStream(cerr,ret2->Out());
						cerr << '\n';
						ToStream(cerr,ret->Out());
						cerr << '\n';
						exit(1);
					}
					ret = ret2;
				}
			}
			if (! ret)	{
				if (found1 && found2)	{
					ret = from;
				}
			}
		}
		return ret;
	}

	Link* ParseGroup(string input, Link* from);
	// a group is an expression of one of the two following forms:
	//
	// 	(Body)Node_name
	// 	(Body)Node_name:Branch_name
	//
	// where Body is a list, and Node_name and Branch_name are 2 strings
	// Node_name may be an empty string
	//
	// Node_name cannot contain the ':' character, but Branch_name can
	// thus, if the group reads "(BODY)A:B:C"
	// then Node_name = "A" and Branch_name = "B:C"

	Link* ParseList(string input, Node* node);
	// a list is an expression of the form X1,X2,...Xn
	// where Xi is a group 

	void RecursiveClone(const Link* from, Link* to);
	// Used by Tree(const Tree* from)
	// only clone the Links, and their mutual relations
	// does not copy the Node or Branch objects

	// deletes the link structure
	// does not delete the Node or Branch objects
	void RecursiveDelete(Link* from);

	void SetRoot(Link* link) {root = link;}

	// data fields
	// just 2 pointers, to the root and to a list of taxa
	Link* root;
	const TaxonSet* taxset;
	int Nlink;
	int Nnode;
	int Nbranch;
};



#endif // TREE_H
