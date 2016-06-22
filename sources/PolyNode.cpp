
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
//		 PolyNode()
// ---------------------------------------------------------------------------------

PolyNode::PolyNode()	{

	up = 0;
	down = 0;
	next = this;
	prev = this;

	label = -1;
	degree = 0;

	size = 0;
	intdepth = 0;
	depth = 0;

	name = "";
	branchLength = 0;
	mProb = 1;
	mTree = 0;

}


// ---------------------------------------------------------------------------------
//		 ~PolyNode()
// ---------------------------------------------------------------------------------

PolyNode::~PolyNode()	{
}


// ---------------------------------------------------------------------------------
//		 deleteRecursive()
// ---------------------------------------------------------------------------------

void PolyNode::deleteRecursive()	{

	if (! IsLeaf())	{

		PolyNode* node = down;
		do	{
			PolyNode* theNext = node->next;
			node->deleteRecursive();
			node = theNext;
		} while (node != down);

	}
	delete this;
}


// ---------------------------------------------------------------------------------
//		 SetSuper(Tree* inTree)
// ---------------------------------------------------------------------------------

void 	PolyNode::SetSuper(PBTree* inTree)	{

	mTree = inTree;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			node->SetSuper(inTree);
			node = node->next;
		} while (node != down);
	}
}
	

// ---------------------------------------------------------------------------------
//		 GetMaxDegree()
// ---------------------------------------------------------------------------------

int PolyNode::GetMaxDegree()	{

	int degree = 0;
	if ((! IsLeaf()) && (! IsRoot()))	{
		PolyNode* node = down;
		do	{
			degree++;
			node = node->next;
		}	while (node != down);
	}
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			int tmp = node->GetMaxDegree();
			if (tmp > degree)	{
				degree = tmp;
			}
			node = node->next;
		}	while (node != down);
	}
	return degree;
}

// ---------------------------------------------------------------------------------
//		 ComputeDegree()
// ---------------------------------------------------------------------------------

void PolyNode::ComputeDegree()	{

	degree = 0;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			degree++;
			node->ComputeDegree();
			node = node->next;
		}	while (node != down);
	}
}

// ---------------------------------------------------------------------------------
//		 CheckDegree()
// ---------------------------------------------------------------------------------

int PolyNode::CheckDegree()	{

	PolyNode* node = down;
	int temp = 0;
	do	{
		temp++;
		node = node->next;
	}	while (node != down);

	return (degree == temp);
}


// Attach and Detach() do not involve news nor deletes

// ---------------------------------------------------------------------------------
//		 AttachTo(PolyNode* inNode)
// ---------------------------------------------------------------------------------


PolyNode* PolyNode::AttachTo(PolyNode* inNode)	{

	// check that this is a root node
	if (! IsRoot())	{
		cerr << "error : trying to attach non root node\n";
		cerr.flush();
	}

	if (inNode->IsLeaf())	{

		inNode->down = this;
		up = inNode;
		next = this;
		prev = this;
	}

	else	{
		PolyNode* temp = inNode->down->prev;
		inNode->down->prev = this;
		prev = temp;
		temp->next = this;
		next = inNode->down;
		up = inNode;
	}

	inNode->degree++;

	return inNode;
}

// ---------------------------------------------------------------------------------
//		 Detach()
// ---------------------------------------------------------------------------------

// at the end : up = 0, and prev = next = this
//

void PolyNode::Detach()	{

	if (up != 0)	{

		up->degree --;

		if (prev == this)	{
			up->down = 0;
			up = 0;
		}

		else	{
			if (up->down == this)	{
				up->down = prev;
			}
			prev->next = next;
			next->prev = prev;
			next = this;
			prev = this;
			up = 0;
		}
	}
	else	{
		cerr << "error : detaching a root node \n";
	}

}



// ---------------------------------------------------------------------------------
//		 DetachRecursive()
// ---------------------------------------------------------------------------------

// at the end : up = 0, and prev = next = this
//

void PolyNode::DetachRecursive()	{

	if (up != 0)	{

		up->degree --;

		if (prev == this)	{
			up->DetachRecursive();
			up->down = 0;
			up = 0;
		}

		else	{
			if (up->down == this)	{
				up->down = prev;
			}
			prev->next = next;
			next->prev = prev;
			next = this;
			prev = this;
			up = 0;
		}
	}
	else	{
		cerr << "error : detaching a root node \n";
	}

}


PolyNode* PolyNode::Simplify()	{

	PolyNode* returnNode = 0;
	if (! IsLeaf())	{
		PolyNode* node = down;

		int temp = 0;
		do	{
			temp ++;
			node = node->next;
		}	while (node != down);

		for (int i=0; i<temp; i++)	{
			PolyNode* nextnode = node->next;
			node->Simplify();
			node = nextnode;
		}
	
		if (down->next == down)	{
			if (IsRoot())	{
				down->up = 0;
				returnNode = down;
			}
			else	{
				PolyNode* mydown = down;
				PolyNode* myup = up;
				Detach();
				mydown->Detach();
				mydown->AttachTo(myup);	
				mydown->branchLength += branchLength;
			}
		}		
	}
	return returnNode;
}

// ---------------------------------------------------------------------------------
//		 RegisterWithData(int& currentLabel, string* SpeciesNames)
// ---------------------------------------------------------------------------------

// first dichotomise the tree
// then recursive traversal of the tree
// 		check that each non leaf node has no name
// 		and give it a label between Ntaxa and 2*Ntaxa-1
//
//		likewise, check that each leaf node has the name of a species present in datafile
//		and them give it the label corresponding to this species' rank order in data file


int PolyNode::RegisterWithData(int& currentLabel, int* found, string* SpeciesNames, int Ntaxa)	{

	if (! IsLeaf())	{

		label = currentLabel;
		currentLabel ++;
		int cont = 1;
		PolyNode* node = down;

		do	{
			cont = node->RegisterWithData(currentLabel, found, SpeciesNames, Ntaxa);
			node = node->next;
		}	while( cont && (node!=down) );
		return cont;
	}
	else	{

		if (name == "")	{
			if ((0>label) || (label>=Ntaxa))	{
				cerr << "assuming labels in tree file, but leaf labels out of range\n";
				cerr << "label equals " << label << " (should be between 0 and " << Ntaxa<< ")\n";

				return 0;
			}
			else	{
				return 1;
			}
		}

		int i=0;
		while ( (i<Ntaxa) && (name != SpeciesNames[i]) )	{
			i++;
		}
		if (i == Ntaxa)	{
			cerr << "found a leaf node with name " << name << " but no corresponding species name in data file \n";
			return 0;
		}
		label = i;
		found[i] = 1;
		return 1;
	}
}

		
// ---------------------------------------------------------------------------------
//		 SortLeaves()
// ---------------------------------------------------------------------------------

int PolyNode::SortLeaves()	{

	if (IsLeaf())	{
		
		if ((!mTree->mParam) ||(name == ""))	{
			return label;
		}
		int k = 0;
		while ((k<mTree->mParam->Ntaxa) && (name != mTree->mParam->SpeciesNames[k])) k++;
		if (k == mTree->mParam->Ntaxa)	{
			cerr << "error in PolyNode::SortLeaves: did not find taxon : " << name << '\n';
			exit(1);
		}
		return k;
	}
	else	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			node = node->next;
			degree++;
		}	while(node!=down);
		int retval[degree];
		PolyNode* nodelist[degree];
		int k = 0;
		do	{
			nodelist[k] = node;
			retval[k] = node->SortLeaves();
			k++;
			node = node->next;
		}	while(node!=down);
		for (int i=0; i<degree; i++)	{
			for (int j=degree-1; j>i; j--)	{
				if (retval[i] > retval[j])	{
					int tmp = retval[i];
					retval[i] = retval[j];
					retval[j] = tmp;
					PolyNode* tmpnode = nodelist[i];
					nodelist[i] = nodelist[j];
					nodelist[j] = tmpnode;
				}
			}
		}
		for (int i=0; i<degree-1; i++)	{
			nodelist[i]->next = nodelist[i+1];
		}
		nodelist[degree-1]->next = nodelist[0];
		for (int i=1; i<degree; i++)	{
			nodelist[i]->prev = nodelist[i-1];
		}
		nodelist[0]->prev = nodelist[degree-1];
		down = nodelist[0];
		return retval[0];
	}	
	return 0;
}

string PolyNode::SortLeavesAlphabetical()	{

	if (IsLeaf())	{
		return name;
	}
	else	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			node = node->next;
			degree++;
		}	while(node!=down);
		string* retval = new string[degree];
		PolyNode* nodelist[degree];
		int k = 0;
		do	{
			nodelist[k] = node;
			retval[k] = node->SortLeavesAlphabetical();
			k++;
			node = node->next;
		}	while(node!=down);
		for (int i=0; i<degree; i++)	{
			for (int j=degree-1; j>i; j--)	{
				if (retval[i] > retval[j])	{
					string tmp = retval[i];
					retval[i] = retval[j];
					retval[j] = tmp;
					PolyNode* tmpnode = nodelist[i];
					nodelist[i] = nodelist[j];
					nodelist[j] = tmpnode;
				}
			}
		}
		for (int i=0; i<degree-1; i++)	{
			nodelist[i]->next = nodelist[i+1];
		}
		nodelist[degree-1]->next = nodelist[0];
		for (int i=1; i<degree; i++)	{
			nodelist[i]->prev = nodelist[i-1];
		}
		nodelist[0]->prev = nodelist[degree-1];
		down = nodelist[0];
		string ret = retval[0];
		delete[] retval;
		return ret;
	}	
	return "";
}



// ---------------------------------------------------------------------------------
//		 SetNames()
// ---------------------------------------------------------------------------------

void PolyNode::SetNames()	{

	if (IsLeaf())	{
		name = mTree->mParam->SpeciesNames[label];
	}
	else	{
		PolyNode* node = down;
		do	{
			node->SetNames();
			node = node->next;
		}	while(node!=down);
	}		
}

// ---------------------------------------------------------------------------------
//		 IsDichotomous()
// ---------------------------------------------------------------------------------

int PolyNode::IsDichotomous()	{

	int answer = 1;
	if (! IsLeaf())	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			answer &= node->IsDichotomous();
			degree ++;
			node = node->next;
		}	while (node != down);
		answer &= (degree == 2);
	}
	return answer;
}


// ---------------------------------------------------------------------------------
//		 Dichotomise()
// ---------------------------------------------------------------------------------

int PolyNode::Dichotomise()	{

	if (! IsLeaf())	{

		degree = 0;
		PolyNode* node = down;
		int cont = 1;
		do	{
			cont = node->Dichotomise();
			node = node->next;
			degree ++;
		}	while ( cont && (node != down) );

		if ( (! cont) || (degree == 1))	{
			return 0;
		}
		else	{

			down = down->prev;
			while (down->next != down->prev)	{

				PolyNode* newNode = new PolyNode();
				newNode->SetTree(mTree);
				PolyNode* node1 = down->prev;
				PolyNode* node2 = down->prev->prev;
				node1->Detach();
				node2->Detach();
				node1->AttachTo(newNode);
				node2->AttachTo(newNode);
				newNode->AttachTo(this);

			}
		}

		return 1;
	}
	return 1;
}


// ---------------------------------------------------------------------------------
//		ComputeMinLeafLabel
// ---------------------------------------------------------------------------------

int PolyNode::ComputeMinLeafLabel()	{

	int min = -1;
	if (IsLeaf())	{
		min = label;
	}
	else	{
		PolyNode* node = down;
		do	{
			int a = node->ComputeMinLeafLabel();
			if ((min == -1) || (min > a))	{
				min = a;
			}
			node = node->next;
		}	while (node != down);
	}
	return min;
}

					

// ---------------------------------------------------------------------------------
//		GetSize
// ---------------------------------------------------------------------------------

int PolyNode::GetSize()	{

	if (IsLeaf())	{
		size = 1;
	}
	else	{
		PolyNode* node = down;
		size = 0;
		do	{
			size += node->GetSize();
			node = node->next;
		}
		while (node != down);
	}
	return size;
}


// ---------------------------------------------------------------------------------
//		GetSpeciesNames(string* name, int index)
// ---------------------------------------------------------------------------------

void PolyNode::GetSpeciesNames(string* names, int& index)	{

	if (IsLeaf())	{
		names[index++] = name; 
	}
	else	{
		PolyNode* node = down;
		size = 0;
		do	{
			node->GetSpeciesNames(names, index);
			node = node->next;
		}
		while (node != down);
	}
}


// ---------------------------------------------------------------------------------
//		GetSpeciesNames(string* name)
// ---------------------------------------------------------------------------------

void PolyNode::GetSpeciesNames(string* names)	{

	if (IsLeaf())	{
		names[label] = name; 
	}
	else	{
		PolyNode* node = down;
		size = 0;
		do	{
			node->GetSpeciesNames(names);
			node = node->next;
		}
		while (node != down);
	}
}


// ---------------------------------------------------------------------------------
//		GetIntDepth
// ---------------------------------------------------------------------------------

int PolyNode::GetIntDepth()	{

	intdepth = 0;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			int tmp = node->GetIntDepth();
			if (intdepth < tmp)	{
				intdepth = tmp;
			}
			node = node->next;
		}
		while (node != down);
		intdepth ++;
	}
	return intdepth;
}

// ---------------------------------------------------------------------------------
//		GetLength
// ---------------------------------------------------------------------------------

double PolyNode::GetLength()	{

	double total = branchLength;

	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			total +=node->GetLength();
			node = node->next;
		}
		while (node != down);
	}
	return total;
}



// ---------------------------------------------------------------------------------
//		GetDepth
// ---------------------------------------------------------------------------------

double PolyNode::GetDepth()	{

	depth = branchLength;

	if (! IsLeaf())	{
		PolyNode* node = down;
		double temp = 0;
		do	{
			if (temp < node->GetDepth())	{
				temp = node->depth;
			}
			node = node->next;
		}
		while (node != down);
		depth += temp;
	}
	return depth;
}


// ---------------------------------------------------------------------------------
//		ComputeSizeAndDepth
// ---------------------------------------------------------------------------------

void PolyNode::ComputeSizeAndDepth()	{

	depth = branchLength;

	if (IsLeaf())	{
		size = 1;
	}
	else	{
		size = 0;
		PolyNode* node = down;
		double temp = 0;
		do	{
			node->ComputeSizeAndDepth();
			if (temp < node->depth)	{
				temp = node->depth;
			}
			size += node->size;
			node = node->next;
		}
		while (node != down);
		depth += temp;
	}
}

// ---------------------------------------------------------------------------------
//		SetFormalBranchLengths
// ---------------------------------------------------------------------------------

void PolyNode::SetFormalBranchLengths(double scale, int order)	{

	// assumes that intdepth are updated

	if (IsLeaf())	{
		branchLength = order * scale;
	}
	else	{

		if (order == 1)	{
			cerr << "error in PolyNode::SetFormalBranchLength : order = 0 \n";
			exit(1);
		}

		branchLength = scale;
		PolyNode* node = down;
		do	{
			node->SetFormalBranchLengths(scale, order-1);
			node = node->next;
		}
		while (node != down);
	}
}




// ---------------------------------------------------------------------------------
//		 GetLeafSet
// ---------------------------------------------------------------------------------


void PolyNode::GetLeafSet(Bipartition& bp)	{

	if (IsLeaf())	{
		bp.mArray[label] = 0;
	}
	else	{
		PolyNode* node = down;
		do	{
			node->GetLeafSet(bp);
			node = node->next;
		}
		while (node != down);
	}
}

// ---------------------------------------------------------------------------------
//		 BipartitionPruning
// ---------------------------------------------------------------------------------

Bipartition PolyNode::BipartitionPruning(BipartitionList* inBList)	{

	if (! mTree)	{
		cerr << "error : BipartitionPruning called on polynode not connected to a consensus\n";
		exit(1);
	}

	Bipartition temp(mTree->GetParameters());
	if (IsLeaf())	{
		temp.SetTaxon(label);
		if (inBList)	{
			inBList->FastInsert(temp, 1.0, branchLength);
		}
	}

	else	{
		PolyNode* node = down;
		do	{
			temp |= node->BipartitionPruning(inBList);
			node = node->next;
		} while (node != down);
		if ((! IsRoot()) && inBList)	{
			inBList->FastInsert(temp, 1.0, branchLength);
		}
	}
	return temp;
}


// ---------------------------------------------------------------------------------
//		 BipartitionPruningWithSupports
// ---------------------------------------------------------------------------------

Bipartition PolyNode::BipartitionPruningWithSupports(BipartitionList* inBList)	{

	if (! mTree)	{
		cerr << "error : BipartitionPruning called on polynode not connected to a consensus\n";
		exit(1);
	}

	Bipartition temp(mTree->GetParameters());
	if (IsLeaf())	{
		temp.SetTaxon(label);
		inBList->FastInsert(temp, mProb, mProb * branchLength);
	}

	else	{
		PolyNode* node = down;
		do	{
			temp |= node->BipartitionPruningWithSupports(inBList);
			node = node->next;
		} while (node != down);
		if ((! IsRoot()) && inBList)	{
			inBList->FastInsert(temp, mProb, mProb * branchLength);
		}
	}
	return temp;
}

// ---------------------------------------------------------------------------------
//		 GetBipartition
// ---------------------------------------------------------------------------------

Bipartition PolyNode::GetBipartition()	{
	return BipartitionPruning(0);
}


// ---------------------------------------------------------------------------------
//		 FindNode
// ---------------------------------------------------------------------------------

PolyNode* PolyNode::FindNode(Bipartition& inPartition, Boolean& Orientation)	{

	Bipartition temp(mTree->GetParameters());
	return FindNodePruning(inPartition, Orientation, temp);
}


// ---------------------------------------------------------------------------------
//		 FindNodePruning
// ---------------------------------------------------------------------------------

PolyNode* PolyNode::FindNodePruning(Bipartition& inPartition, Boolean& Orientation, Bipartition& outPartition)	{

	PolyNode* returnNode = 0;
	if (! mTree)	{
		cerr << "Error : FindNodePruning called on node with no mTree\n";
		exit(1);
	}
	else	{
		if (IsLeaf())	{
			outPartition.SetTaxon(label);
		}
		else	{
			PolyNode* node = down;

			do	{
				Bipartition temp(mTree->GetParameters());
				returnNode = node->FindNodePruning(inPartition, Orientation, temp);
				outPartition |= temp;
				node = node->next;
			}	while ((node != down) && ! returnNode);
		}

		if (! returnNode)	{
			if (outPartition == inPartition)	{
				Orientation = true;
				returnNode = this;
			}
			else	{
				if (outPartition == ! inPartition)	{
					Orientation = false;
					returnNode = this;
				}
			}
		}
	}
	return returnNode;
}

// ---------------------------------------------------------------------------------
//		 Analyse
// ---------------------------------------------------------------------------------

int PolyNode::Analyse(Bipartition& inPartition)	{

	value = 2;
	if (IsLeaf())	{
		if (inPartition.GetTaxonStatus(label) == -1)	{
			cerr << "error in PolyNode::Analyse\n";
			exit(1);
		}
		if (inPartition.GetTaxonStatus(label) == -3)	{
			cerr << "error in PolyNode::Analyse\n";
			exit(1);
		}
		if (inPartition.GetTaxonStatus(label) == -10)	{
			cerr << "error in PolyNode::Analyse\n";
			exit(1);
		}
		value = inPartition.GetTaxonStatus(label);
	}
	else	{

		PolyNode* node = down;
		do	{
			node->Analyse(inPartition);

			if (node->value <= -1)	{
				value = -2;
			}
			if (value != -2)	{
				if (value ==  2)	{
					value = node->value;
				}
				else	{
					if (value != node->value)	{
						value = -1;
					}
				}
			}
			 node = node->next;
		} while (node != down);
	}

	return value;
}

// ---------------------------------------------------------------------------------
//		 Insert(Bipartition&, double prob)
// ---------------------------------------------------------------------------------

void PolyNode::Insert(Bipartition& inPartition, double prob, double length)	{

	if (IsLeaf())	{
		cerr << "error in PolyNode::insert : can not do that on leaf ! \n";
		exit(1);
	}

	else	{

		if (value == -2)	{
			PolyNode* node = down;
			while (node->value != -2 && node->value != -1)	{
				node = node->next;

				if (node == down)	{
					cerr << "error in PolyNode::Insert : turning around ! \n";
					exit(1);
				}
			}

			node->Insert(inPartition, prob, length);
		}
		else	{
			// value should be -1
			if (value != -1)	{
				cerr  << "error in PolyNode::Insert : value should be -1\n";
				inPartition.WriteToStream(cerr,1);
				cerr << '\n';
				exit(1);
			}

			int neighbor = 0;

			if (! IsRoot())	{
				neighbor = next->value;
			}


			// try to gather all the nodes in down whose value is 1-neighbor

			PolyNode* insert = new PolyNode();
			insert->SetTree(GetTree());

			PolyNode* node = down;

			while(node->value != 3)		{
				if (node->value == 1 - neighbor)	{

					node->value = 3;
					PolyNode* nextNode = node->next;

					node->Detach();
					node->AttachTo(insert);

					node = nextNode;

				}
				else	{
					node->value = 3;
					node = node->next;
				}
			}

			if (insert->down->next == insert->down)	{
				PolyNode* temp = insert->down;
				temp->Detach();
				temp->AttachTo(this);
				delete insert;
				if (prob != 1)	{
					cerr << "error : probability is not 1\n";
					cerr << prob << '\n';
					exit(1);
				}
				temp->branchLength = length;
			}
			else	{
				insert->AttachTo(this);
				insert->mProb = prob;
				insert->branchLength = length;
			}
		}
	}
}



// ---------------------------------------------------------------------------------
//		 FlipFlop
// ---------------------------------------------------------------------------------

PolyNode* PolyNode::FlipFlop()	{

	PolyNode* returnNode = 0;
	if (IsRoot())	{
		cerr << "error in polynode flipflop : called on root\n";
	}
	else	{

		PolyNode* theUp = up;
		if (! theUp->IsRoot())	{
			returnNode = theUp->FlipFlop();
		}

		Detach();
		if (theUp->IsRoot())	{
			theUp->mProb = mProb;
			// theUp->branchLength = branchLength;
		}
		PolyNode* insert = 0;

		if (theUp->down->next == theUp->down)	{
			insert = theUp->down;
			theUp->down = 0;
			theUp->up = 0;
			theUp->mProb = 0;
			theUp->branchLength = 0;
			returnNode = theUp;
		}
		else	{
			insert = theUp;
		}
		if (IsLeaf())	{
			cerr << "error in polynode flip flop : called on leaf\n";
		}
		else	{
			insert->up = 0;
			insert->AttachTo(this);
			up = 0;
		}
		insert->branchLength += branchLength;
		branchLength = 0;
	}
	return returnNode;
}


// ---------------------------------------------------------------------------------
//		 Phylip(ostream& os)
// ---------------------------------------------------------------------------------

void PolyNode::Phylip(ostream& os, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{

	if (IsLeaf())	{
		if (withSpeciesNames)	{
			os << name;
		}
		else	{
			os << label;
		}
		if (withLengths)	{
			os << ':' << Decimal(branchLength,6);
			// os  << ':' << ((double) ((int) (precision * branchLength))) / precision ;
		}
	}
	else	{
		os  << '(';

		PolyNode* node = down->prev;
		do	{
			node->Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);
			if (node->prev != down->prev)	{
				os << ',';
			}
			node = node->prev;
		}	while (node != down->prev);

		/*
		PolyNode* node = down;
		do	{
			node->Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);
			if (node->next!= down)	{
				os << ',';
			}
			node = node->next;
		}	while (node != down);
		*/
		os <<  ')';
		if (! IsRoot())	{
			if (withInternalLabels)	{
				os << label;
			}
			if (withProbs)	{
				os  << Decimal(mProb, 2);
				// os  << ((double) ((int) (precision * mProb))) / precision;
			}
			if (withLengths)	{
				os << ':' << Decimal(branchLength,6);
				// os << ':' << ((double) ((int) (precision * branchLength))) / precision ;
			}
		}
	}
}


void PolyNode::ResetLabels()	{

	label = -1;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			node->ResetLabels();
			node = node->next;
		} while (node!=down);
	}
}


// ---------------------------------------------------------------------------------
//		 TranslateLabels()
// ---------------------------------------------------------------------------------

void 	PolyNode::TranslateLabels(int offset)	{

	if (IsLeaf())	{
		label += offset;
	}
	else	{
		PolyNode* node = down;
		do	{
			node->TranslateLabels(offset);
			node= node->next;
		}	while (node != down);
	}
}
