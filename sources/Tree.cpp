
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <list>
#include <sstream>
#include <fstream>
#include "Tree.h"
#include "Random.h"

bool NewickTree::simplify = false;

void NewickTree::ToStream(ostream& os) const {
	if (simplify)	{
		ToStreamSimplified(os,GetRoot());
	}
	else	{
		ToStream(os,GetRoot());
	}
	os << ";\n";
}

double NewickTree::ToStreamSimplified(ostream& os, const Link* from) const {

	if (!from->isLeaf())	{
		if (from->Next()->Next() == from)	{
			double tot = ToStreamSimplified(os,from->Next()->Out());		
			tot += atof(GetBranchName(from).c_str());
			return tot;
		}
		else	{
			os << '(';
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				double tmp = ToStreamSimplified(os,link->Out());
				os << ':' << tmp;
				if (link->Next() != from)	{
					os << ',';
				}
			}
			os << ')';
		}
	}
	else	{
	}
	os << GetNodeName(from);
	/*
	if (!from->isRoot())	{
		string brval = GetBranchName(from);
		if (brval != "")	{
			os << ':' << brval;
		}
	}
	*/
	if (from->isRoot())	{
		return 0;
	}
	return atof(GetBranchName(from).c_str());
}

void NewickTree::ToStream(ostream& os, const Link* from) const {

	if (!from->isLeaf())	{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStream(os,link->Out());
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	os << GetNodeName(from);
	if (!from->isRoot())	{
		string brval = GetBranchName(from);
		if (brval != "")	{
			os << ':' << brval;
		}
	}
}
		
Tree::Tree()	{
	root = 0;
	taxset = 0;
}

Tree::Tree(const TaxonSet* intaxset)	{

	taxset = intaxset;
	root = new Link();
	root->InsertOut(root);
	Node* node = new Node();
	root->SetNode(node);

}

void Tree::MakeRandomTree()	{

	int Ntaxa = taxset->GetNtaxa();
	int* included = new int[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		included[i] = 0;
	}

	int* triplet = new int[3];
	rnd::GetRandom().DrawFromUrn(triplet,3,Ntaxa);

	included[triplet[0]] = 1;
	included[triplet[1]] = 1;
	included[triplet[2]] = 1;

	Link* link1 = new Link();
	Link* link2 = new Link();
	Link* link3 = new Link();

	Link* linko1 = new Link();
	Link* linko2 = new Link();
	Link* linko3 = new Link();

	root->SetNext(link1);

	link1->SetNext(link2);
	link2->SetNext(link3);
	link3->SetNext(root);

	linko1->SetOut(link1);
	link1->SetOut(linko1);
	linko2->SetOut(link2);
	link2->SetOut(linko2);
	linko3->SetOut(link3);
	link3->SetOut(linko3);

	linko1->SetIndex(triplet[0]);
	linko2->SetIndex(triplet[1]);
	linko3->SetIndex(triplet[2]);

	link1->SetNode(root->GetNode());
	link2->SetNode(root->GetNode());
	link3->SetNode(root->GetNode());

	Branch* b1 = new Branch();
	Branch* b2 = new Branch();
	Branch* b3 = new Branch();

	link1->SetBranch(b1);
	link2->SetBranch(b2);
	link3->SetBranch(b3);

	linko1->SetBranch(b1);
	linko2->SetBranch(b2);
	linko3->SetBranch(b3);

	Node* n1 = new Node(taxset->GetTaxon(triplet[0]));
	Node* n2 = new Node(taxset->GetTaxon(triplet[1]));
	Node* n3 = new Node(taxset->GetTaxon(triplet[2]));

	linko1->SetNode(n1);
	linko2->SetNode(n2);
	linko3->SetNode(n3);

	for (int i=3; i<Ntaxa; i++)	{
		// ToStream(cerr);

		int choose = ((int) ((Ntaxa - i) * rnd::GetRandom().Uniform()));

		for (int j=0; j<=choose; j++)	{
			if (included[j])	{
				choose++;
			}
		}
		if (included[choose])	{
			cerr << "error in random tree\n";
			exit(1);
		}

		included[choose] = 1;

		Link* linkup = new Link();
		Link* linkupo = new Link();
		Link* linkdown = new Link();
		Link* linkdowno = new Link();
		linkdowno->SetNext(linkup);

		Branch* bdown = new Branch();
		Branch* bup = new Branch();
		Node* ndown = new Node(taxset->GetTaxon(choose));
		Node* nup = new Node();

		linkup->SetOut(linkupo);
		linkupo->SetOut(linkup);
		linkdown->SetOut(linkdowno);
		linkdowno->SetOut(linkdown);

		linkdown->SetNode(ndown);
		linkdowno->SetNode(nup);
		linkup->SetNode(nup);

		linkdown->SetBranch(bdown);
		linkdowno->SetBranch(bdown);
		linkup->SetBranch(bup);
		linkupo->SetBranch(bup);

		linkdown->SetIndex(choose);
		
		// choose a place in the tree: any node except root
		int n = ((int) ((2*i-3) * rnd::GetRandom().Uniform())) + 1;
		Link* tmp = 0;
		Link* fromdown = ChooseNode(root,tmp,n);
		Link* fromup = fromdown->Out()->Next();

		// attach (as in gibbs)
		Attach(linkdown,linkup,fromdown,fromup);
	}

	delete[] included;
	delete[] triplet;
}

Tree::Tree(const Tree* from)	{
	taxset = from->GetTaxonSet();
	root = new Link(from->root);
	root->InsertOut(root);
	RecursiveClone(from->root,root);
}

void Tree::RecursiveClone(const Link* from, Link* to)	{
	Node* node = new Node(from->GetNode());
	to->SetNode(node);
	const Link* linkfrom = from->Next();	
	Link* linkto = to;	
	while (linkfrom != from)	{
		Link* newnext = new Link(linkfrom); // newnext points to same node and branch as linkfrom
		newnext->SetNode(node);
		linkto->Insert(newnext);
		Link* newout = new Link(linkfrom->Out()); // idem, same node and branch as linkfrom->Out()
		newout->InsertOut(newnext);
		Branch* branch = new Branch(linkfrom->GetBranch());
		newnext->SetBranch(branch);
		newout->SetBranch(branch);
		RecursiveClone(linkfrom->Out(),newout);
		linkfrom = linkfrom->Next();
		linkto = linkto->Next();
	}	
}


void Tree::RecursiveDelete(Link* from)	{
	if (from)	{
		Link* link = from->Next();
		while (link != from)	{
			delete link->Out()->GetNode();
			delete link->GetBranch();
			RecursiveDelete(link->Out());
			Link* keep = link->Next();
			delete link;
			link = keep;
		}
		delete link;
	}
}

Tree::~Tree()	{
	if (root)	{
		RecursiveDelete(root);	
		root = 0;
	}
}


void Tree::DeleteNextLeaf(Link* previous){
	Link* link = previous->Next();
	if(!link->Out()->isLeaf()){
		cout << "Bad call of DeleteNextLeaf, it must be call on a link pointing on a leaf\n";
		exit(1);
	}
	previous->Next()->Next()->AppendTo(previous);
	delete link->Out()->GetNode();
	delete link->Out();
	delete link->GetBranch();
	delete link;
}
	
void Tree::DeleteUnaryNode(Link* from){
	if(!from->isUnary()){
		cout << "Bad call of DeleteUnaryNode, node is not unary\n";
		exit(1);
	}
	if(from->isRoot()){
		Link* newroot = from->Next()->Out();
		newroot->SetBranch(from->GetBranch());		
		delete from->Next()->GetBranch();
		delete from->GetNode();
		delete from->Next();
		delete from;
		root = newroot;
		root->InsertOut(root);
	}
	else{
		ostringstream sum;
		sum << atof((from->GetBranch()->GetName()).c_str()) + atof((from->Next()->GetBranch()->GetName()).c_str());
		from->GetBranch()->SetName(sum.str());
		from->Next()->Out()->SetBranch(from->GetBranch());
		from->Out()->InsertOut(from->Next()->Out());
		delete from->Next()->GetBranch();
		delete from->GetNode();
		delete from->Next();
		delete from;
	}
}


void Tree::EraseInternalNodeName()	{
	EraseInternalNodeName(GetRoot());
}

void Tree::EraseInternalNodeName(Link * from)	{
	if (! from->isLeaf())	{
		from->GetNode()->SetName("");
	}
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		EraseInternalNodeName(link->Out());
	}
}

void Tree::RegisterWith(const TaxonSet* intaxset,int myid)	{
	taxset = intaxset;
	int tot = 0;
	if(!RegisterWith(taxset,GetRoot(),tot)){
		// if (myid == 0)
		cout << "There is no match between the tree and the sequences.\n";
		cerr << "problem with : " << taxset->GetNtaxa () << '\n';
		ToStream(cerr);
		exit(1);
	}
	if (tot != taxset->GetNtaxa())	{
		// if (myid == 0) {
			cerr << "error : non matching number of taxa : " << tot << '\t' << taxset->GetNtaxa() << '\n';
			cerr << "some taxa in the dataset are not present in the tree\n";
		// }
		exit(1);
	}
	SetIndices();
	CheckIndices(GetRoot());
}

bool Tree::RegisterWith(const TaxonSet* taxset, Link* from, int& tot)	{
	if (from->isLeaf())	{
		int i = taxset->GetTaxonIndex(from->GetNode()->GetName());
		if (i != -1)	{
			from->GetNode()->SetIndex(i);
			tot++;
		}
		return(i != -1);
	}
	else{
		Link* previous = from;
		while(previous->Next() != from)	{
			if(RegisterWith(taxset,previous->Next()->Out(),tot)){
				previous = previous->Next();
			}
			else{
				// cout << "delete !!\n";
				DeleteNextLeaf(previous);
			}
		}
		if(from->isUnary()){
			// cerr << "DELETE UNARY NODE\n";
			DeleteUnaryNode(from);
		}
		return !(from->isLeaf());
	}
}


/*
void Tree::ToStreamWithIndices(ostream& os, const Link* from) {

	if (!from->isLeaf())	{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStreamwithIndices(os,link->Out());
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	os << from->GetNode()->GetIndex();
	if (!from->isRoot())	{
		os << from->GetIndex() << ':' << from->GetBranch()->GetIndex() << from->Out()->GetIndex();
	}
}
		
void Tree::RegisterWithIndicesFromNames(const Link* from)	{

	string tmp = from->GetNode()->GetName();
	if (! IsInt(tmp))	{
		cerr << "error in register with indices from names\n";
		exit(1);
	}
	int temp = atoi(tmp.c_str());
	from->GetNode()->SetIndex(temp);

	if (!from->isRoot())	{
		os << from->GetIndex() << ':' << from->GetBranch()->GetIndex() << from->Out()->GetIndex();
	}
	if (!from->isLeaf())	{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStreamwithIndices(os,link->Out());
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	os << from->GetNode()->GetIndex();
}
*/

typedef list<string>::const_iterator csit;

Tree::Tree(string filename)	{

	root = 0;
	taxset = 0;
	ifstream is(filename.c_str());
	if (!is)	{
		cout << "cannot find file : " << filename << '\n';
	}
	ReadFromStream(is);

	if (! CheckRootDegree())	{
		cerr << "error: root should be of degree three\n";
		cerr << "i.e. tree should have following format: (A,B,C) and not (A,(B,C));\n";
		exit(1);
	}

	if (! RecursiveCheckDegree(GetRoot()))	{
		cerr << "error: input tree is not bifurcating\n";
		exit(1);
	}
}

Tree::Tree(istream& is)	{

	root = 0;
	taxset = 0;
	ReadFromStream(is);

	if (! CheckRootDegree())	{
		cerr << "error: root should be of degree tree\n";
		cerr << "i.e. tree should have following format: (A,B,C) and not (A,(B,C));\n";
		exit(1);
	}

	if (! RecursiveCheckDegree(GetRoot()))	{
		cerr << "error: input tree is not bifurcating\n";
		exit(1);
	}
}

void Tree::ReadFromStream(istream& is)	{

	RecursiveDelete(GetRoot());
	string expr = "";
	int cont = 1;
	while (cont)	{
		string s;
		is >> s;
		if (s.length() == 0)	{
			cerr << "in tree: null string\n";
			exit(1);
		}
		unsigned int k = 0;
		while ((k<s.length()) && (s[k] != ';'))	k++;
		expr += s.substr(0,k);
		cont = (! is.eof()) && (k==s.length());
	}
	SetRoot(ParseGroup(expr,0));
	if (taxset)	{
		RegisterWith(taxset);
	}
}

Link* Tree::ParseList(string input, Node* node)	{
	
	try	{

		// parse input as a list of strings separated by ','
		list<string> lst;
		int n = input.size();
		int k = 0;
		int brack = 0;
		int b = 0;
		while (k<n)	{
			char c = input[k];
			if (c == '(') brack++;
			if (c == ')') brack--;
			if ((!brack) && (c == ','))	{
				lst.push_back((string) (input.substr(b,k-b)));
				b = k+1;
			}

			if (brack<0)	{
				cout << "in parse list : too many )\n";
				cout << input.substr(0,k) << '\n';
				cout << input << '\n';
				throw;
			}
			k++;
		}
		if (brack)	{
			cout << "in parse list : too many (\n";
			cout << input << '\n';
			throw;
		}
		lst.push_back(input.substr(b,k-b));

		// make a circular single link chain around the node
		// with one link for each term of the list
		// and call parse group on each term
		Link* firstlink = new Link;
		Link* prevlink = firstlink;
		firstlink->SetNode(node);
		for (csit i=lst.begin(); i!=lst.end(); i++)	{ 
			Link* link = new Link;
			link->SetNode(node);
			link->AppendTo(prevlink);
			ParseGroup(*i,link);
			prevlink = link;
		}
		firstlink->AppendTo(prevlink);
		return firstlink;
	}
	catch(...)	{
		cout << "exit in parse list\n";
		exit(1);
	}
}

Link* Tree::ParseGroup(string input, Link* from)	{
	
	try	{

		// parse input as (body)nodeval:branchval
		
		string body = "";
		unsigned int k = 0;
		if (input[0] == '(')	{
			int brack = 1;
			k = 1;
			while ((k<input.length()) && brack)	{
				char c = input[k];
				if (c == '(') brack++;
				if (c == ')') brack--;
				k++;
			}
			if (brack)	{
				cout << "in parse group: too many (\n";
				cout << input << '\n';
				throw;
			}
			body = input.substr(1,k-2);
		}
		
		int b = k;
		while ((k<input.length()) && (input[k]!=':')) k++;
		string nodeval = input.substr(b,k-b);

		string branchval = "";
		if (k<input.length())	{
			branchval = input.substr(k+1, input.length() - k);
		}

		// make a new node and a new branch
		Node* node = new Node(nodeval);

		// call parse body
		Link* link = 0;
		if (body != "")	{
			link = ParseList(body,node);
		}
		else	{
			link = new Link;
			link->SetNode(node);
		}
		if (from)	{
			Branch* branch = new Branch(branchval);
			link->SetBranch(branch);
			from->SetBranch(branch);
			link->InsertOut(from);
		}
		return link;
	}
	catch(...)	{
		cout << "exit in parse group\n";
		exit(1);
	}
}

void Tree::Subdivide(Link* from, int Ninterpol)	{

	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		Subdivide(link->Out(),Ninterpol);
	}
	// if ((! from->isLeaf()) && (! from->isRoot()))	{
	if (! from->isRoot())	{
		double l = atof(from->GetBranch()->GetName().c_str());
		if (l <= 0)	{
			cerr << "warning : non strictly positive branch length : " << l << '\n';
			cerr << "correcting and setting to 0.001\n";
			l = 0.001;
		}

		ostringstream s;
		s << l / Ninterpol;
		
		delete from->GetBranch();
		
		Link* current = from;
		Link* final = from->Out();
		int i = 0;
		while (i < Ninterpol-1)	{
			Link* link1 = new Link;
			Link* link2 = new Link;
			Branch* newbranch = new Branch(s.str());
			Node* newnode = new Node();
			current->SetBranch(newbranch);
			link1->SetNext(link2);
			link2->SetNext(link1);
			link1->SetBranch(newbranch);
			link1->SetNode(newnode);
			link2->SetNode(newnode);
			current->SetOut(link1);
			link1->SetOut(current);
			current = link2;
			i++;
		}
		current->SetOut(final);
		final->SetOut(current);
		Branch* newbranch = new Branch(s.str());
		final->SetBranch(newbranch);
		current->SetBranch(newbranch);
	}
}

/*
void Tree::Print(ostream& os, const Link* from)	const {

	if (!from->isLeaf())	{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			Print(os,link->Out());
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	os << GetNodeName(from);
	if (!from->isRoot())	{
		string brval = GetBranchName(from);
		if (brval != "")	{
			os << ':' << brval;
		}
	}
}
		
void Tree::Print(ostream& os)	const {
	Print(os,GetRoot());
	os << ";\n";
}
*/

bool Tree::CheckRootDegree(int testdegree)	{
	bool ret = true;
	int degree = 0;
	const Link* from = GetRoot();
	for (const Link* link = from->Next(); link!= from; link=link->Next())	{
		degree++;
	}
	if (degree != testdegree)	{
		ret = false;
	}
	return ret;
}

bool Tree::RecursiveCheckDegree(const Link* from, int testdegree)	{

	bool ret = true;
	int degree = 0;
	if (! from->isRoot())	{
		degree++;
	}
	for (const Link* link = from->Next(); link!= from; link=link->Next())	{
		degree++;
	}
	if (degree != testdegree)	{
		ret = false;
	}
	for (const Link* link = from->Next(); link!= from; link=link->Next())	{
		if (! link->Out()->isLeaf())	{
			ret &= RecursiveCheckDegree(link->Out(), testdegree);
		}
	}
	return ret;
}

Link* Tree::Detach(Link* down, Link* up)	{

	bool foundup = false;
	Link* fromdown = 0;
	Link* downout = down->Out();
	int degree = 0;
	for (Link* link=downout->Next(); link!=downout; link=link->Next())	{
		degree++;
		if (link == up)	{
			foundup = true;
		}
		else if (! link->isRoot())	{
			fromdown = link->Out();
		}
		if (link->isRoot())	{
			cerr << "link is root!\n";
		}
	}
	if (degree != 2)	{
		cerr << "error in detach: node not of degree 2\n";
		cerr << degree << '\n';
		// cerr << down->GetIndex() << '\t' << up->GetIndex() << '\n';
		exit(1);
	}
	if (! fromdown)	{
		cerr << "error in Detach: fromdown not found\n";
		exit(1);
	}
	if (!foundup)	{
		cerr << "error in Detach: dit not find up\n";
		exit(1);
	}

	Link* fromout = fromdown->Out();
	Link* upout = up->Out();
	Link* linkprev = 0;
	for (Link* link=upout->Next(); link!= upout; link=link->Next())	{
		if (link->Next() == upout)	{
			linkprev = link;
		}
	}
	linkprev->SetNext(fromout);
	fromout->SetNext(upout->Next());
	fromout->SetNode(upout->Next()->GetNode());
	upout->SetNext(upout);
	upout->SetNode(0);
	up->SetNext(downout);
	downout->SetNext(up);
	return fromdown;
}

void Tree::Attach(Link* down, Link* up, Link* todown, Link* toup)	{

	bool found = false;
	Link* linkprev = toup;
	Link* upout = up->Out();
	Link* downout = down->Out();
	for (Link* link=toup->Next(); ((!found) && (link!=toup)); link=link->Next())	{
		if (link == todown->Out())	{
			linkprev->SetNext(upout);
			upout->SetNext(link->Next());
			found = true;
		}
		linkprev = link;
	}
	if (! found)	{
		cerr << "error in attach: not found\n";
		exit(1);
	}
	if (downout->Next() != up)	{
		cerr << "error in attach: downoutnext != up\n";
		exit(1);
	}
	up->SetNext(todown->Out());
	todown->Out()->SetNext(downout);
	todown->Out()->SetNode(up->GetNode());
	up->Out()->SetNode(toup->GetNode());
}


void Tree::NNIturn(Link* from){

	Link* up = from->Next()->Out();
	Link* down = up->Next()->Out();
	Detach(down, up);
	Link* todown = from->Next()->Next()->Out();
	from->Knit();
	Attach(down, up, todown, from);
}

int Tree::CountInternalNodes(const Link* from)	{
	int total = 0;
	if (! from->isLeaf())	{
	// if ((! from->isLeaf()) && (! from->isRoot()))	{
		total = 1;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += CountInternalNodes(link->Out());
		}
	}
	return total;
}

Link* Tree::ChooseInternalNode(Link* from, Link*& fromup, int& n)	{
	if (from->isLeaf())	{
		return 0;
	}
	Link* ret = 0;
	if (!n)	{
		ret = from;
	}
	else	{
		n--;
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			if (!ret)	{
				Link* tmp = ChooseInternalNode(link->Out(), fromup, n);
				if (tmp)	{
					ret = tmp;
				}
			}
		}
		if (ret && (! fromup))	{
			fromup = from;
		}
	}
	return ret;
}

int Tree::CountNodes(const Link* from)	{
	int total = 1;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += CountNodes(link->Out());
	}
	return total;
}

Link* Tree::ChooseNode(Link* from, Link*& fromup, int& n)	{
	Link* ret = 0;
	if (!n)	{
		ret = from;
	}
	else	{
		n--;
		for (Link* link=from->Next(); (!ret && (link!=from)); link=link->Next())	{
			Link* tmp = ChooseNode(link->Out(), fromup, n);
			if (tmp)	{
				ret = tmp;
			}
		}
		if (ret && (! fromup))	{
			fromup = from;
		}
	}
	return ret;
}

Link* Tree::ChooseLinkAtRandom()	{
	int n = CountInternalNodes(GetRoot());
	int choose = (int) (n * rnd::GetRandom().Uniform());
	Link* tmp = 0;
	Link* newrootnext = ChooseInternalNode(GetRoot(),tmp,choose);
	return newrootnext;
}

void Tree::RootAt(Link* newrootnext)	{
	if (newrootnext->GetNode() != GetRoot()->GetNode())	{
		Link* prev = 0;
		for (Link* link=root->Next(); link!= root; link=link->Next())	{
			if (link->Next() == root)	{
				prev = link;
			}
		}
		prev->SetNext(root->Next());
		Link* newprev = 0;
		for (Link* link=newrootnext->Next(); link!= newrootnext; link=link->Next())	{
			if (link->Next() == newrootnext)	{
				newprev = link;
			}
		}
		newprev->SetNext(root);
		root->SetNext(newrootnext);
		root->SetNode(newrootnext->GetNode());
	}
}

void Tree::RootAtRandom()	{
	Link* newrootnext = ChooseLinkAtRandom();
	RootAt(newrootnext);
}

int Tree::DrawSubTree(Link*& down, Link*& up)	{

	int nodestatus[Nnode];
	for (int i=0; i<Nnode; i++)	{
		nodestatus[i] = 1;
	}
	int m = Nnode;
	nodestatus[root->GetNode()->GetIndex()] = 0;
	m--;
	for (const Link* link= root->Next(); link!=root; link=link->Next())	{
		nodestatus[link->Out()->GetNode()->GetIndex()] = 0;
		m--;
	}
	int choose = (int) (m * rnd::GetRandom().Uniform());
	int k = 0;
	while ((choose < Nnode) && (k < choose))	{
		if (!nodestatus[k])	{
			choose ++;
		}
		k++;
	}
	while ((choose < Nnode) && (! nodestatus[choose]))	{
		choose++;
	} 
	if (choose == Nnode)	{
		cerr << "error in draw sub tree: overflow\n";
		exit(1);
	}

	down = 0;
	up = 0;
	for (Link* link=root->Next(); link!=root; link=link->Next())	{
		if (! up)	{
			GrepNode(link,down,up,choose);
		}
	}

	if ((! down) || (down->isRoot()))	{
		cerr << "error in tree::drawsubtree: down\n";
		cerr << down << '\n';
		exit(1);
	}
	if ((! up) || (up->isLeaf()))	{
		cerr << "error in tree::drawsubtree: up\n";
		cerr << up << '\n';
		exit(1);
	}
	return 0;
}

void Tree::GrepNode(Link* from, Link*& down, Link*& up, int choose)	{
	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! up)	{
			if (link->Out()->GetNode()->GetIndex() == choose)	{
				up = from;
				down = link->Out();
			}
			GrepNode(link->Out(),down,up,choose);
		}
	}
}

/*
int Tree::DrawSubTree(Link*& down, Link*& up)	{

	int n = 0;
	for (const Link* link= root->Next(); link!=root; link=link->Next())	{
		n += CountNodes(link->Out());
		n-=4;
	}
	int choose = (int) (n * rnd::GetRandom().Uniform());
	down = 0;
	up = 0;
	for (const Link* link= root->Next(); link!=root; link=link->Next())	{
		if (! down)	{
			choose++;
			down = ChooseNode(link->Out(),up,choose);
		}
	}
	if ((! down) || (down->isRoot()))	{
		cerr << "error in tree::drawsubtree: down\n";
		cerr << down << '\n';
		exit(1);
	}
	if ((! up) || (up->isLeaf()))	{
		cerr << "error in tree::drawsubtree: up\n";
		cerr << up << '\n';
		 exit(1);
	}
	return 0;
}
*/
