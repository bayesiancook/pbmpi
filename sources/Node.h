
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


class Node	{
	public :

	Node* up;
	Node* left;
	Node* right;
	int label;
	double branchLength;


	void _init()	{
		up = 0;
		left = 0;
		right = 0;
		label = 0;
		branchLength = 2;
	}

	Node()	{_init();}

	virtual ~Node (){};

	int isRoot() const	{ return (! up);	}

	int isLeaf() const	{ return (! left); 	}

	void swap()	{

		if (! this)	{
			cerr << " Node::swap()	: called on null Node\n";
		}
		Node* temp = left;
		left = right;
		right = temp;
	}

	void upToRight()	{
		if (! this) cerr << " Node::upToRight()	: called on null Node\n";
		if (isRoot()) cerr << " Node::upToRight()	: cannot do that on root !\n";
		if (up->isRoot())	{
			if (up->right == this)	right = up->left;
			else					right = up->right;
			right->up = this;
			up = 0;
			right->branchLength += branchLength;
			branchLength = 0;
			return;
		}
		if ((up->left) == this)		up->upToLeft();
		else 						up->upToRight();

		right = up;
		right->up = this;
		right->branchLength += branchLength;
		branchLength =0;

	}

	void upToLeft()	{
		if (! this) cerr << " Node::upToLeft()	: called on null Node\n";
		if (isRoot()) cerr <<  " Node::upToLeft()	: cannot do that on root !\n";
		if (up->isRoot())	{
			if (up->left == this)	left = up->right;
			else					left = up->left;
			left->up = this;
			up = 0;
			left->branchLength += branchLength;
			branchLength =0;
			return;
		}
		if ((up->right) == this)	up->upToRight();
		else 						up->upToLeft();
		left = up;
		left->up = this;
		left->branchLength += branchLength;
		branchLength =0;

	}


	void rootAt(Node* node)	{		// new root is just upstream node
		if (! node) cerr << " Node::rootAt()	: called on null Node\n";
		if (node->isRoot()) cerr << " Node::rootAt()	: cannot not reroot upstream root!\n";
		else 	{
			Node* Up = node->up;
			if (!Up->isRoot())	{
				if ( (Up->right) == node)		{
					Up->upToRight();
				}
				else	{
					Up->upToLeft();
				}
				this->left = Up;
				Up->up = this;
				this->right = node;
				node->up = this;
			}
		}

		// place the root at random on the corresponding internal edge
		double totalLength = left->branchLength + right->branchLength;
		double x = 0.5;
		left->branchLength = x;
		right->branchLength = totalLength - x;
	}

	void rootAt2(Node* node)	{		// new root is just upstream node
		if (! node) cerr << " Node::rootAt()	: called on null Node\n";
		if (node->isRoot()) cerr << " Node::rootAt()	: cannot not reroot upstream root!\n";
		else 	{
			Node* Up = node->up;
			if (!Up->isRoot())	{
				if ( (Up->right) == node)		{
					Up->upToRight();
				}
				else	{
					Up->upToLeft();
				}
				this->left = Up;
				Up->up = this;
				this->right = node;
				node->up = this;
			}
		}

		// place the root next to left
		double totalLength = left->branchLength + right->branchLength;
		left->branchLength = 0;
		right->branchLength = totalLength;
	}

	int getSize()	{			// number of leaves
		if (! this) cerr << " Node::getSize()	: called on null Node\n";
		if (isLeaf() ) 	{return 1;}
		else	{return (left->getSize() + right->getSize());}
	}

	int getFullSize()	{		// total number of nodes, root included
		if (! this) cerr << " Node::getFullSize()	: called on null Node\n";
		if (isLeaf() ) 	{return 1;}
		else	{return (1+ left->getFullSize() + right->getFullSize());}
	}

	double getLength()	{		// basal branch length  included
		if (! this) cerr << " Node::getLength()	: called on null Node\n";
		if (isLeaf())	{return branchLength;}
		else	{return (branchLength + left->getLength() + right->getLength());}
	}


	void Phylip(ostream& os)	{

		if (isLeaf())	{
			os << label ;
		}
		else	{
			os  << '(';
			left->Phylip(os);
			os << ',';
			right->Phylip(os);
			os <<  ')';
		}
	}

	void PhylipWithBranchLengths(ostream& os, int Nsite)	{

		if (isLeaf())	{
			os << label ;
		}
		else	{
			os  << '(';
			left->PhylipWithBranchLengths(os, Nsite);
			os << ',';
			right->PhylipWithBranchLengths(os, Nsite);
			os <<  ')';
		}
		if (! isRoot()) {
			os << ':' << branchLength ;
		}
	}

	void PhylipWithBranchLengthsAndSpeciesNames(ostream& os, string* names, int Nsite)	{

		if (isLeaf())	{
			os << names[label] ;
		}
		else	{
			os  << '(';
			left->PhylipWithBranchLengthsAndSpeciesNames(os, names, Nsite);
			os << ',';
			right->PhylipWithBranchLengthsAndSpeciesNames(os, names, Nsite);
			os <<  ')';
		}
		if (! isRoot()) {
			os << ':' << branchLength ;
		}
	}


	void PhylipWithSpeciesNames(ostream& os, string* names, int nameLength = -1)	{

		if (isLeaf())	{
			if (nameLength == -1)	{
				os <<  names[label];
			}
			else	{
				os << setw(nameLength) << names[label];
			}
		}
		else	{
			os  << '(';
			left->PhylipWithSpeciesNames(os,names);
			os << ',';
			right->PhylipWithSpeciesNames(os,names);
			os <<  ')';
		}
	}

};
