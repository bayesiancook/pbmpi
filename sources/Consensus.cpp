
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


// ---------------------------------------------------------------------------------
//		 Consensus
// ---------------------------------------------------------------------------------

Consensus::Consensus(TreeList* inTreeList, double* probarray, double cutoff)	{

	mTreeList = inTreeList;
	mParam = mTreeList->GetParameters();
	mTreeList = inTreeList;
	mBList = new BipartitionList(mTreeList,probarray);
	BipartitionList* reduced = MakeConsensus(cutoff);
	SetNames();
	RootAt0Alphabetical();
	mRoot->SortLeavesAlphabetical();
	Trichotomise();

	BipartitionList* test = new BipartitionList(this);
	if (!test->IsEqualTo(reduced))	{
		cerr << "error in consensus making: non matching bp lists\n";
		exit(1);
	}
}


// ---------------------------------------------------------------------------------
//		 Consensus
// ---------------------------------------------------------------------------------

Consensus::Consensus(BipartitionList* inBList, double cutoff)	{

	mBList = inBList;
	mParam = mBList->GetParameters();
	mTreeList = 0;
	BipartitionList* reduced = MakeConsensus(cutoff);
	SetNames();
	RootAt0Alphabetical();
	mRoot->SortLeavesAlphabetical();
	Trichotomise();

	BipartitionList* test = new BipartitionList(mParam);
	test->PruneWithSupports(this);
	if (!test->IsEqualTo(reduced))	{
		cerr << "error in consensus making: non matching bp lists\n";
		exit(1);
	}
}

// ---------------------------------------------------------------------------------
//		 ~Consensus
// ---------------------------------------------------------------------------------

Consensus::~Consensus()	{
	
}

// ---------------------------------------------------------------------------------
//		 MakeConsensus
// ---------------------------------------------------------------------------------


BipartitionList*
Consensus::MakeConsensus(double cutoff)	{

	// sort the list, and apply threshold
	mBList->Sort();
	mBList->Truncate(cutoff);

	// make the reduced list
	// contain only the major mutually compatible bipartitions
	mReducedList = new BipartitionList(mParam);
	mReducedList->mWeight = mBList->mWeight;

	for (int i=0; i<mBList->GetSize(); i++)	{
		Bipartition& temp = (*mBList)[i];
		if (mReducedList->IsCompatibleWith(temp))	{
			mReducedList->Insert(temp,mBList->mWeightArray[i], mBList->mLengthArray[i]);
		}
	}
	if (mReducedList->GetSize() > (2*mParam->Ntaxa - 3))	{
		cerr << "error : too many partitions to build a consensus : " << mReducedList->GetSize() << '\n';
		exit(1);
	}
	if (mReducedList->GetSize() == 0)	{
		cerr << "error : cannot build a consensus with an empty list\n";
		exit(1);
	}

	// make the star tree
	mRoot = new PolyNode();
	mRoot->SetTree(this);
	
	PolyNode** LeafTable = new PolyNode*[mParam->Ntaxa];

	for (int i=0; i<mParam->Ntaxa; i++)	{
		PolyNode* node = new PolyNode();
		node->SetTree(this);
		node->SetLabel(i);
		node->AttachTo(mRoot);
		LeafTable[i] = node;
	}


	// make the consensus out of bipartition list
	for (int i=0; i<mReducedList->GetSize(); i++)	{
		Insert((*mReducedList)[i],mReducedList->GetProb(i), mReducedList->GetLength(i));
	}
	return mReducedList;
}

// ---------------------------------------------------------------------------------
//		 Insert
// ---------------------------------------------------------------------------------

void
Consensus::Insert(Bipartition& inPartition, double prob, double length)	{
	mRoot->Analyse(inPartition);
	mRoot->Insert(inPartition, prob, length);
}
