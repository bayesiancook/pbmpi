
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "phylo.h"


inline double approx(double f)	{
	return ((double) ((int) (100 * f))) / 100;
}
		
// ---------------------------------------------------------------------------------
//		FinishCreation()
// ---------------------------------------------------------------------------------

void PBTree::FinishCreation()	{

	mRoot->SetSuper(this);
	SetLabelOffset(0);
	if (mParam)	{
		SetNames();
	}
	// mRoot->SortLeavesAlphabetical();
}

void PBTree::SetNames()	{
	if (! mParam)	{
		cerr << "error in PBTree:SetNames\n";
		exit(1);
	}
	mRoot->SetNames();
}

	
// ---------------------------------------------------------------------------------
//		PBTree(TaxaParameters* inParam)
// ---------------------------------------------------------------------------------

PBTree::PBTree(TaxaParameters* inParam)	{

	Name = "T";
	mParam = inParam;
	mNodeList = 0;
	fromPB = 0;
	mRoot = 0;
	mLabelOffset = 0;
}


// ---------------------------------------------------------------------------------
//		PBTree(string filename)
// ---------------------------------------------------------------------------------

PBTree::PBTree(string filename)	{

	mParam = 0;
	Name = "T";
	ifstream is(filename.c_str());
	if (! is)	{
		cerr << "error : non existing file : " << filename << '\n';
		exit(1);
	}
	ReadFromStream(is,1);
	FinishCreation();

}


// ---------------------------------------------------------------------------------
//		PBTree(PhyloBayes*)
// ---------------------------------------------------------------------------------

PBTree::PBTree(PhyloBayes* pb, TaxaParameters* inparam)	{

	/*
	Name = "T";
	fromPB = 1;
	if (inparam)	{
		mParam = inparam;
	}
	else	{
		mParam = new TaxaParameters(pb->mParam);
	}

	mNodeList = new PolyNode[pb->mParam->Nnode];
	PolyNode* p = mNodeList;
	AmphiNode* fromp = pb->tree;

	for (int i=0; i< pb->mParam->Nnode; i++)	{

		p->mPBTree = this;
		p->label = i;
		p->branchLength = pb->BL[fromp->label];

		p->up = (fromp->up) ? &mNodeList[fromp->up->label] : 0;
		p->down = (fromp->left) ? &mNodeList[fromp->left->label] : 0;

		if (p->down)	{
			p->down->next = (fromp->right) ? &mNodeList[fromp->right->label] : 0;
			if (! p->down)	{
				cerr << "error in PBTree:PBTree(PhyloBayes*) : right is null, but not left\n";
				exit(1);
			}
			p->down->prev = p->down->next;
			p->down->next->prev = p->down;
			p->down->next->next = p->down;
		}
		else	{
			p->name = pb->mParam->SpeciesNames[fromp->label];
		}
		p++;
		fromp++;
	}
	mRoot = &mNodeList[pb->root->label];
	FinishCreation();
	*/
}


// ---------------------------------------------------------------------------------
//		PBTree(PBTree* from)
// ---------------------------------------------------------------------------------

PBTree::PBTree(PBTree* from)	{

	mParam = from->mParam;
	Name = from->Name;
	mLabelOffset = from->mLabelOffset;
	fromPB = 0;
	
	mRoot = new PolyNode();
	Clone(mRoot, from->mRoot);
	mRoot->SetSuper(this);
}
	

// ---------------------------------------------------------------------------------
//		Clone()
// ---------------------------------------------------------------------------------

void PBTree::Clone(PolyNode* node, PolyNode* fromNode)	{

	node->label = fromNode->label;
	node->branchLength = fromNode->branchLength;
	node->name = fromNode->name;

	if (fromNode->IsLeaf())	{
		node->down = 0;
	}
	else	{
		PolyNode* fromfirst = fromNode->down;
		PolyNode* first = new PolyNode();
		first->up = node;
		node->down = first;
		Clone(first, fromfirst);
		
		PolyNode* fromnext = fromfirst->next;
		PolyNode* thiscurrent = first;

		while (fromnext  != fromfirst)	{
			PolyNode* thisnext = new PolyNode();
			thisnext->up = node;
			thisnext->prev = thiscurrent;
			thiscurrent->next = thisnext;
			Clone(thisnext, fromnext);
			thiscurrent = thisnext;
			fromnext = fromnext->next;
		}
		thiscurrent->next = first;
		first->prev = thiscurrent;
	}
}
	

// ---------------------------------------------------------------------------------
//		~PBTree()
// ---------------------------------------------------------------------------------

PBTree::~PBTree()	{
	
	if (fromPB)	{
		delete[] mNodeList;
	}
	else	{
		mRoot->deleteRecursive();
	}
}


// ---------------------------------------------------------------------------------
//		GetLength()
// ---------------------------------------------------------------------------------

double PBTree::GetLength()	{
	return mRoot->GetLength();
}

// ---------------------------------------------------------------------------------
//		GetSize()
// ---------------------------------------------------------------------------------

int PBTree::GetSize()	{
	return mRoot->GetSize();
}


// ---------------------------------------------------------------------------------
//		GetSpeciesNames
// ---------------------------------------------------------------------------------

void PBTree::GetSpeciesNames(string* name)	{
	int index = 0;
	mRoot->GetSpeciesNames(name, index);
}


// ---------------------------------------------------------------------------------
//		Dichotomise
// ---------------------------------------------------------------------------------

int PBTree::Dichotomise(int arrangelabels)	{

	int returnValue = 0;
	if (! IsDichotomous())	{
		returnValue = mRoot->Dichotomise();
	}
	if (arrangelabels)	{ // this is for normal approximation
		int found = 0;
		PolyNode* node = mRoot->down;
		do	{
			if (node->label == -1)	{
				if (found)	{
					cerr << "error in PolyNode::Dichotomise: found 2 nodes with label -1\n";
				// 	exit(1);
				}
				found = 1;
				node->label = mRoot->label;
				mRoot->label++;
			}
			node = node->next;
		} while (node != mRoot->down);
	}
	return returnValue;
}

// ---------------------------------------------------------------------------------
//		IsDichotomous
// ---------------------------------------------------------------------------------

int PBTree::IsDichotomous()	{

	return mRoot->IsDichotomous();

}

// ---------------------------------------------------------------------------------
//		Trichotomise
// ---------------------------------------------------------------------------------

void PBTree::Trichotomise()	{

	double lengthbefore = GetLength();
	PolyNode* node = mRoot->down;
	if (node->next->next == node)	{
		if (!node->IsLeaf())	{
			node = node->next;
		}
		else	{
			if (node->next->IsLeaf())	{
				cerr << "error : cannot trichotomise a two species tree\n";
				exit(1);
			}
		}

		PolyNode* temp = node->next;
		node->Detach();
		// PolyNode* temp = mRoot->down;
		mRoot->down = 0;
		// delete mRoot;
		mRoot = temp;
		mRoot->up = 0;
		node->branchLength += mRoot->branchLength;
		mRoot->branchLength = 0;
		node->AttachTo(mRoot);
	}
	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in trichotomise: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		exit(1);
	}
}

// ---------------------------------------------------------------------------------
//		Trichotomise2
// ---------------------------------------------------------------------------------

/*
void PBTree::Trichotomise()	{

	double lengthbefore = GetLength();
	PolyNode* node = mRoot->down;
	if (node->next->next == node)	{
		if (!node->IsLeaf())	{
			node = node->next;
		}
		else	{
			if (node->next->IsLeaf())	{
				cerr << "error : cannot trichotomise a two species tree\n";
				exit(1);
			}
		}

		PolyNode* temp = node->next;
		node->Detach();
		// PolyNode* temp = mRoot->down;
		mRoot->down = 0;
		// delete mRoot;
		mRoot = temp;
		mRoot->up = 0;
		node->branchLength += mRoot->branchLength;
		mRoot->branchLength = 0;
		node->AttachTo(mRoot);
	}
	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in trichotomise: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		exit(1);
	}
}
*/

// ---------------------------------------------------------------------------------
//		SetFormalBranchLengths
// ---------------------------------------------------------------------------------

void	PBTree::SetFormalBranchLengths(double scale)	{

	int order = mRoot->GetIntDepth();
	mRoot->SetFormalBranchLengths(scale, order+1);
}


// ---------------------------------------------------------------------------------
//		operator!=
// ---------------------------------------------------------------------------------

Boolean	PBTree::operator!=(PBTree& from)	{

	return ! operator==(from);
}

// ---------------------------------------------------------------------------------
//		operator==
// ---------------------------------------------------------------------------------

Boolean	PBTree::operator==(PBTree& from)	{

	if (! mParam)	{
		cerr << "error : operator == called on a tree without taxa param\n";
		exit(1);
	}
	
	BipartitionList bblist(mParam);
	BipartitionPruning(&bblist);
	BipartitionList frombblist(mParam);
	from.BipartitionPruning(&frombblist);

	int n = bblist.GetSize();
	int i= 0;
	int cont = 1;
	while ((i<n) && cont)	{
		Bipartition& bp = bblist[i];
		if (bp.IsInformative())	{
			int j = 0;
			while ( (j < n) && (bp != frombblist[j]) )	{
				j++;
			}
			if (j == n)	{
				cont = 0;
			}
		}
		i++;
	}
	return cont;
}


// ---------------------------------------------------------------------------------
//		 FindNode()
// ---------------------------------------------------------------------------------

PolyNode* PBTree::FindNode(Bipartition& inPartition, Boolean& Orientation)	{

	if (! mParam)	{
		cerr << "error in PBTree::FindNode: taxa param is nil\n";
		exit(1);
	}
	return mRoot->FindNode(inPartition , Orientation);
}


// ---------------------------------------------------------------------------------
//		 FindNode()
// ---------------------------------------------------------------------------------

PolyNode* PBTree::FindNode(Bipartition& inPartition)	{

	if (! mParam)	{
		cerr << "error in PBTree::FindNode: taxa param is nil\n";
		exit(1);
	}
	Boolean temp;
	return FindNode(inPartition, temp);
}

// ---------------------------------------------------------------------------------
//		 BipartitionPruning
// ---------------------------------------------------------------------------------

void PBTree::BipartitionPruning(BipartitionList* inList)	{

	Trichotomise();
	mRoot->BipartitionPruning(inList);
	inList->Modulo();

}

// ---------------------------------------------------------------------------------
//		 GetRootBipartition
// ---------------------------------------------------------------------------------

Bipartition PBTree::GetRootBipartition()	{

	return mRoot->down->GetBipartition();
}


// ---------------------------------------------------------------------------------
//		 GetLeaf()
// ---------------------------------------------------------------------------------

PolyNode* PBTree::GetLeaf(int inLabel)	{

	return &mNodeList[inLabel];
}

// ---------------------------------------------------------------------------------
//		RegisterWithData
// ---------------------------------------------------------------------------------

int PBTree::RegisterWithParam(TaxaParameters* taxaparam)	{

	return RegisterWithData(taxaparam->SpeciesNames, taxaparam->Ntaxa);
}

int PBTree::RegisterWithData(string* SpeciesNames, int Ntaxa)	{


	int size = GetSize();
	string* name = new string[size];
	GetSpeciesNames(name);
	
	for (int i=0; i<size; i++)	{
		int k = 0;
		while ((k<Ntaxa) && (name[i] != SpeciesNames[k])) k++;
		if (k == Ntaxa)	{
			// cerr << "eliminating " << name[i] << '\n';
			Eliminate(name[i]);
		}
	}

	int runningLabel = Ntaxa;
	int* found = new int[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		found[i] = 0;
	}
	int returnValue = mRoot->RegisterWithData(runningLabel, found, SpeciesNames, Ntaxa);
	for (int i=0; i<Ntaxa; i++)	{
		if (! found[i])	{
			cerr << "error when registering tree with data: did not find taxon " << SpeciesNames[i] << " in tree\n";
			exit(1);
		}
	}
	return returnValue;
}


// ---------------------------------------------------------------------------------
//		SetLabelOffset
// ---------------------------------------------------------------------------------

void	PBTree::SetLabelOffset(int offset)	{

	int currentoffset = mRoot->ComputeMinLeafLabel();
	int shift = offset - currentoffset;
	if (shift)	{
		TranslateLabels(shift);
	}
	mLabelOffset = offset;
}

// ---------------------------------------------------------------------------------
//		 TranslateLabels(int offset)
// ---------------------------------------------------------------------------------

void	PBTree::TranslateLabels(int offset)	{

	mRoot->TranslateLabels(offset);
}


// ---------------------------------------------------------------------------------
//		ReadFromStream(istream&, int header)
// ---------------------------------------------------------------------------------

int PBTree::ReadFromStream(istream& is, int header)	{

	mRoot = new PolyNode();
	PolyNode* currentNode = mRoot;
	/*
	char c = is.peek();
	while (c != '(')	{
		c = is.get();
		if (c == '[')	{
			do {
				c = is.get();
			}
			while (c != ']');
		}
		c = is.peek();
	}
	*/
	int ok = ReadPhylip(is,currentNode);
	SetLabelOffset(0);
	mRoot->SetSuper(this);
	return (ok && (GetSize() > 1));
}


// ---------------------------------------------------------------------------------
//		WriteToStream(ostream&, int header)
// ---------------------------------------------------------------------------------

void PBTree::Phylip(ostream& os, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{
	WriteToStream(os, 0, withLengths, withProbs, withSpeciesNames, withInternalLabels);
}


void PBTree::WriteToStream(ostream& os, int header, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{
	mRoot->Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);
	if (withInternalLabels && IsDichotomous())	{
		os << mRoot->label;
	}
	os << ";\n";
}


// ---------------------------------------------------------------------------------
//		ReadPhylip
// ---------------------------------------------------------------------------------

	
int PBTree::ReadPhylip(istream& is, PolyNode* currentNode)	{

// 	int verbose = 1;
	PolyNode* node = 0;
	double d;

	int END = 0;

	char c = is.peek();
	while ((c == ' ') || (c == '\n') || (c == '\t'))	{
		is.get();
		c = is.peek();
	}


	while( (! is.eof()) && (! END) )	{

		string s;
		char c = is.peek();

		// cerr << c ;

		switch (c)	{

			case '\n'  :
			case ' '   :
			case '\t'  :

				is >> c;

				if (verbose)	{
					cerr << c;
				}

				break;


			case '('   :

				is >> c;
				if (verbose)	{
					cerr << c;
				}

				// make a new polynode
				node = new PolyNode();
				node->AttachTo(currentNode);
				currentNode = node;

				break;


			case ','   :

				is >> c;
				if (verbose)	{
					cerr << c;
				}
				node = new PolyNode();
				if (currentNode->IsRoot())	{
					cerr << "error in reading tree file : found a coma not within brackets\n";
					exit(1);
				}
				node->AttachTo(currentNode->Up());
				currentNode = node;

				break;


			case ')' :

				is >> c;
				if (verbose)	{
					cerr << c;
				}
				if (currentNode->IsRoot())	{
					cerr << "error in reading tree file: found more closing than opening brackets\n";
					exit(1);
				}

				currentNode = currentNode->Up();

				// check whether is followed by a digit, a ':', a ')' or  a ',' (anything else will throw exception)
				c = is.peek();

				switch(c)	{

					case '0':
					case '1':
					case '2':
					case '3':
					case '4':
					case '5':
					case '6':
					case '7':
					case '8':
					case '9':

						// probability
						is >> d;
						if (verbose)	{
							cerr << d;
						}

						currentNode->SetProb(d);

						// then check whether is further followed by ':' or ',' or ')' (anything else exits)

						c = is.peek();

						switch(c)	{

							case ':' 	:

								is >> c;
								if (verbose)	{
									cerr << c;
								}

								// read branchLength
								// and check that it is followed either by ',' or by ')'

								is >> d;
								if (verbose)	{
									cerr << "BL" << d;
								}
								currentNode->SetBranchLength(d);

								c = is.peek();

								if ((c != ',') && (c != ')') && (c != ';'))	{
									cerr << "error\n";
									cerr << c << '\n';
									exit(1);
								}

							break;

							case ',' :
							case ')' :

							break;

							default	:

								cerr << "error in reading tree file : after reading a branchlength of an internal node\n";
								exit(1);

							break;
						}

					break;


					case ':' 	:

						is >> c;
						if (verbose)	{
							cerr << c;
						}

						// read branchLength
						// and check that it is followed either by ',' or by ')'
						is >> d;
						if (verbose)	{
							cerr << "BL" << d;
						}

						currentNode->SetBranchLength(d);

						c = is.peek();

						if ((c != ',') && (c != ')') && (c != ';'))	{
							cerr << "error\n";
							cerr << c << '\n';
							exit(1);
						}

					break;

					case ';' :
					case ',' :
					case ')' :
					break;

					default	:

						cerr << "error in reading tree file\n";
						cerr << "character : " << c << '\n';
						exit(1);

					break;

				}

			break;


			case ';'	:

				is >> c;
				if (verbose)	{
					cerr << "END" << c;
				}

				// close the tree
				// current node should be root
				if (! currentNode->IsRoot())	{
					cerr << "error in reading tree file : lacking closing brackets overall\n";
					exit(1);
				}

				END = 1;

			break;


			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':

				// species label

				s == "";
				do	{
					s += c;
					is >> c;
					if (is.eof())	{
						cerr << "error in reading treefile : unexpected end of file in string\n";
						exit(1);
					}
					c = is.peek();
				}	while ((c != ',') && (c != ':') && (c != ')'));
				if (IsInt(s))	{
					int i = atoi(s.c_str());
					// is >> i;
					if (verbose)	{
						cerr << i;
					}
					node->SetLabel(i);
				}
				else	{
					node->SetName(s);
				}
				
				// check whether is followed by : or , or ) (anything else will throw exception)
				c = is.peek();

				switch(c)	{

					case ':' :

						is >> c;
						if (verbose)	{
							cerr << c;
						}

						// branchlength
						double d;
						is >> d;
						if (verbose)	{
							cerr << "BL" << d;
						}

						node->SetBranchLength(d);

						c = is.peek();

						if ((c != ',') && (c != ')'))	{
							cerr << "error\n";
							cerr << c << '\n';
							exit(1);
						}

					break;

					case ',' :
					case ')' :

					break;

					default	:

						cerr << "error in reading tree file\n";
						cerr << "character : " << c << '\n';
						exit(1);

					break;

				}

			break;

			default :

				// assume this is a species name

				s == "";
				do	{
					s += c;
					is >> c;
					if (is.eof())	{
						return 0;
					}
					c = is.peek();
				}	while ((c != ',') && (c != ':') && (c != ')'));

				currentNode->SetName(s);
				if (verbose)	{
					cerr << '\"' << s << '\"';
				}

				// check whether is followed by : or , (anything else will throw exception)
				c = is.peek();

				switch(c)	{

					case ':' :

						is >> c;
						if (verbose)	{
							cerr << c;
						}

						// branchlength
						double d;
						is >> d;
						if (verbose)	{
							cerr << "BL" << d;
						}

						node->SetBranchLength(d);
					break;



					case ',' :
					case ')':

					break;

					default	:

						cerr << "error in reading tree file : after reading species name\n";
						exit(1);

					break;

				}


			break;

		}
	}
	if (verbose)	{
		cerr << "tree read : ok " << '\n' << '\n';
		cerr.flush();
	}
	return 1;
}


// ---------------------------------------------------------------------------------
//		Eliminate(string species)
// ---------------------------------------------------------------------------------

void PBTree::Eliminate(string species)	{

	Eliminate(mRoot, species);	
	Simplify();
}


// ---------------------------------------------------------------------------------
//		Eliminate(int label)
// ---------------------------------------------------------------------------------

void PBTree::Eliminate(int label)	{

	Eliminate(mRoot, label);
	Simplify();
}


// ---------------------------------------------------------------------------------
//		Eliminate(PolyNode* node, string species)
// ---------------------------------------------------------------------------------

int  PBTree::Eliminate(PolyNode* node, string species)	{

	int returnValue = 1;
	if (node->IsLeaf())	{
		if (node->name == species)	{
			node->DetachRecursive();
			returnValue = 0;
		}
	}
	else	{
		PolyNode* dnode = node->down;
		do 	{
			returnValue &= Eliminate(dnode, species);
			dnode = dnode->next;
		}	while (returnValue && (dnode != node->down));
	}
	return returnValue;
}
	

void PBTree::Simplify()	{

	PolyNode* temp = mRoot->Simplify();
	if (temp)	{
		mRoot = temp;
	}
}


// ---------------------------------------------------------------------------------
//		Eliminate(PolyNode* node, int label)
// ---------------------------------------------------------------------------------

void PBTree::Eliminate(PolyNode* node, int label)	{

	if (node->label == label)	{
		if (! node->IsLeaf())	{
			cerr << "error in Eliminate : cannot eliminate a non leaf node\n";
			exit(1);
		}
		node->Detach();
	}
	else	{
		if (! node->IsLeaf())	{
			PolyNode* dnode = node->down;
			do 	{
				Eliminate(dnode, label);
				dnode = dnode->next;
			}	while (dnode != node->down);
		}
	}
}
	


// ---------------------------------------------------------------------------------
//		 RootAt(Bipartition)
// ---------------------------------------------------------------------------------

void PBTree::RootAt(Bipartition inPartition)	{

	PolyNode* node = FindNode(inPartition);
	if (node)	{
		RootAt(node);	
	}
	else	{
		cerr << "error in PBTree::RootAt(Bipartition) : could not find this outgroup\n";
	}
}

// ---------------------------------------------------------------------------------
//		 RootAt
// ---------------------------------------------------------------------------------

void PBTree::RootAt0Alphabetical()	{

	if (!mParam)	{
		cerr << "error: parameter is not defined\n";
		exit(1);
	}
	Bipartition bp(mParam);
	int min = 0;
	string name = mParam->SpeciesNames[min];
	for (int i=1; i<mParam->Ntaxa; i++)	{
		if (name > mParam->SpeciesNames[i])	{
			name = mParam->SpeciesNames[i];
			min = i;
		}
	}
	bp.mArray[min] = 1;
	bp.Modulo();
	RootAt(bp);
}


// ---------------------------------------------------------------------------------
//		 RootAt
// ---------------------------------------------------------------------------------

void
PBTree::RootAt(PolyNode* inNode)	{

	double lengthbefore = GetLength();
	// inNode should not be root
	if (inNode->IsRoot())	{
		cerr << "error in PolyNode::RootAt : inNode is root\n";
		exit(1);
	}

	if (inNode)	{

		PolyNode* theUp = inNode->Up();
		PolyNode* newRoot = 0;


		if (theUp->IsRoot())	{
			inNode->Detach();
			if (theUp->down->next != theUp->down)	{
				newRoot = new PolyNode();
				newRoot->SetTree(this);
				theUp->AttachTo(newRoot);
				theUp->mProb = inNode->mProb;
			}
			else	{
				newRoot = theUp;
			}
		}

		else	{
			newRoot = theUp->FlipFlop();
			if ( ! newRoot)	{
				newRoot = new PolyNode();
				newRoot->SetTree(this);
			}
			inNode->Detach();
			theUp->AttachTo(newRoot);
			theUp->mProb = inNode->mProb;
		}
		inNode->AttachTo(newRoot);
		mRoot = newRoot;
	}

	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in rerooting: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		// exit(1);
	}

	double totallength = mRoot->down->branchLength + mRoot->down->next->branchLength;
	mRoot->down->branchLength  = totallength / 2;
	mRoot->down->next->branchLength  = totallength / 2;
	
}




