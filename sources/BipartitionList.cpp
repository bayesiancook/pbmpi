
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
//		 BasicAllocation()
// ---------------------------------------------------------------------------------


/*
double BPCompProgressive(string name1, string name2, int burnin, int every, int until)	{

	ifstream is1(name1.c_str());
	ifstream is2(name1.c_str());

	string temp;
	is1 >> temp;
	is2 >> temp;
	// make tree and grep param

	// make 2 lists
	BipartitionList* list1 = new BipartitionList(param);
	BipartitionList* list2 = new BipartitionList(param);

	int i = 1;

}
*/

double BPCompare(string* ChainName, int P, string reftreename, int burnin, int every, int until, int ps, int verbose, int mergeallbp, string OutFile, double cutoff, double conscutoff, bool rootonly, bool bench)	{


	if (!P)	{
		cerr << "error : should specify chains to read\n";
		exit(1);
	}

	if (OutFile == "")	{
		if (P == 1)	{
			OutFile = ChainName[0];
		}
		else	{
			OutFile = "bpcomp";
		}
	}
	BipartitionList** bplist = new BipartitionList*[P];
	TaxaParameters* taxaparam = 0; 
	BipartitionList* refbplist = 0;

	if (verbose)	{
		cout << '\n';
	}
	int Q = 0;
	for (int p=0; p<P; p++)	{
		string name = ChainName[p];
		int ok = 1;
		if (! ifstream((ChainName[p]).c_str()))	{
			name = ChainName[p] + ".treelist";
			if (! ifstream(name.c_str()))	{
				cerr << "Warning : cannot find " << ChainName[p] << " nor " << name << "\n";
				ok = 0;
			}
		}
		if (reftreename != "")	{
			if (! ifstream(reftreename.c_str()))	{
				cerr << "did not find ref tree\n";
				ok = 0;
			}
		}
		if (ok)	{
			if (verbose)	{
				cout << name << " ";
			}
			bplist[Q] = new BipartitionList(name,burnin,every,until, 0, rootonly);
			if (verbose)	{
				cout << ": " << bplist[Q]->Ntree << " trees\n";
			}
			if (bplist[Q]->Ntree)	{
				if (Q)	{
					if (!bplist[Q]->RegisterWith(taxaparam))	{
						cerr << "Warning: taxon set does not match with that of the first tree list. skipped\n";
						delete bplist[Q];
					}
					else	{
						Q++;
					}
				}
				else	{
					taxaparam = bplist[0]->mParam;
					Q++;
				}
			}
			else	{
				delete bplist[Q];
			}

			if (reftreename != "")	{
				refbplist = new BipartitionList(reftreename,0,1,1,0,rootonly);
				if (refbplist->Ntree != 1)	{
					cerr << "error: ref tree file should contain exactly one tree\n";
					exit(1);
				}
				refbplist->RegisterWith(taxaparam);
			}
		}
	}

	if (!Q)	{
		cerr << "\n";
		cerr << "empty tree collection\n";
		cerr << "sorry\n";
		cerr << "\n";
		exit(1);
	}

	P = Q;

	ofstream osp((OutFile + ".bplist").c_str());

	// merge the bipartition lists into mergedbplist
	BipartitionList* mergedbplist = new BipartitionList(bplist[0]->mParam);
	for (int p=0; p<Q; p++)	{
		mergedbplist->Append(bplist[p]);
	}
	if (refbplist)	{
		mergedbplist->Append(refbplist);
	}

	int bpsize = mergedbplist->GetSize();

	double diff[bpsize];
	double refdiff[bpsize];
	double prob[Q][bpsize];
	double length[Q][bpsize];
	double meanprob[bpsize];
	double refprob[bpsize];
	double reflength[bpsize];
	for (int k=0; k<bpsize; k++)	{
		meanprob[k] = 0;
		refprob[k] = 0;
	}

	for (int k=0; k<bpsize; k++)	{
		Bipartition bp = mergedbplist->GetBipartition(k);
		for (int p=0; p<Q; p++)	{
			int i = 0;
			while ( (i < bplist[p]->GetSize()) && (bp != bplist[p]->GetBipartition(i)) )	{
				i++;
			}
			if (i == bplist[p]->GetSize())	{
				prob[p][k] = 0;
				length[p][k] = 0;
			}
			else	{
				prob[p][k] = bplist[p]->GetProb(i);
				length[p][k] = bplist[p]->GetLength(i);
			}
		}
		if (refbplist)	{
			int i = 0;
			while ( (i < refbplist->GetSize()) && (bp != refbplist->GetBipartition(i)) )	{
				i++;
			}
			if (i == refbplist->GetSize())	{
				refprob[k] = 0;
				reflength[k] = 0;
			}
			else	{
				refprob[k] = refbplist->GetProb(i);
				reflength[k] = refbplist->GetLength(i);
			}
		}
	}

	double maxdiff = 0;
	if (refbplist)	{
		// compute diffs
		double meandiff = 0;
		double weight = 0;
		for (int k=0; k<bpsize; k++)	{
			double max = 0;
			double var = 0;
			double mean = 0;
			double meand = 0;
			int count = 0;
			for (int p=0; p<P; p++)	{
				for (int q=p+1; q<P; q++)	{
					double tmp = fabs(prob[p][k] - prob[q][k]);
					if (max < tmp)	{
						max = tmp;
					}
					meand += tmp;
					count++;
				}
				mean += prob[p][k];
				var += prob[p][k] * prob[p][k];
			}
			mean /= P;
			var /= P;
			var -= mean * mean;
			meand /= count;
			meandiff += mean * meand;
			weight += mean;
			diff[k] = max;
			refdiff[k] = fabs(refprob[k] - mean);
			if (maxdiff < max)	{
				maxdiff = max;
			}
			meanprob[k] = mean;
		}
		meandiff /= weight;
		
		if (! bench)	{
			osp << "total number of bp : " << bpsize << '\n';
			osp << "cutoff : " << cutoff << '\n';
			taxaparam->WriteToStream(osp);
			osp << '\n';
		}

		int permut[bpsize];
		for (int k=0; k<bpsize; k++)	{
			permut[k] = k;
		}
		for (int i=0; i<bpsize; i++)	{
			for (int j = bpsize-1; j>i; j--)	{
				if (refdiff[permut[j]] > refdiff[permut[j-1]])	{
					int tmp = permut[j];
					permut[j] = permut[j-1];
					permut[j-1] = tmp;
				}
			}
		}
		if (! bench)	{
			osp << "bipartition\trefdiff\tprobs\n";
			osp << '\n';
		}
		for (int i=0; i<bpsize; i++)	{
			int k = permut[i];
			if (refdiff[k] > 0.01)   {
				mergedbplist->GetBipartition(k).WriteToStream(osp,verbose);
				// mergedbplist->GetBipartition(k).WriteToStream(osp);
				osp << '\t' << (int) (100 * refdiff[k]);
				osp << '\t' << (int) (100 * diff[k]);
				osp << '\t';
				osp << '\t' << (int) (refprob[k]);
				for (int p=0; p<P; p++)	{
					osp << '\t' << (int) (100 * prob[p][k]);
				}
				osp << '\t';
				osp << '\t' << reflength[k];
				for (int p=0; p<P; p++)	{
					osp << '\t' << length[p][k];
				}
				osp << '\n';
			}
		}
		if (! bench)	{
			osp << '\n' << '\n';
		}
	}
	else if (P > 1)	{
		// compute diffs
		double meandiff = 0;
		double weight = 0;
		for (int k=0; k<bpsize; k++)	{
			double max = 0;
			double var = 0;
			double mean = 0;
			double meand = 0;
			int count = 0;
			for (int p=0; p<P; p++)	{
				for (int q=p+1; q<P; q++)	{
					double tmp = fabs(prob[p][k] - prob[q][k]);
					if (max < tmp)	{
						max = tmp;
					}
					meand += tmp;
					count++;
				}
				mean += prob[p][k];
				var += prob[p][k] * prob[p][k];
				/*
				mean += prob[p][k];
				var += prob[p][k] * prob[p][k];
				*/
			}
			mean /= P;
			var /= P;
			var -= mean * mean;
			meand /= count;
			meandiff += mean * meand;
			// meandiff += mean * sqrt(var);
			weight += mean;
			diff[k] = max;
			if (maxdiff < max)	{
				maxdiff = max;
			}
			meanprob[k] = mean;
		}
		meandiff /= weight;
		// meandiff = sqrt(meandiff);
		
		if (! bench)	{
			osp << "total number of bp : " << bpsize << '\n';
			osp << "cutoff : " << cutoff << '\n';
			taxaparam->WriteToStream(osp);
			osp << '\n';
		}

		if (rootonly)	{
			// sort bipartitions by decreasing diff
			int permut[bpsize];
			for (int k=0; k<bpsize; k++)	{
				permut[k] = k;
			}
			for (int i=0; i<bpsize; i++)	{
				for (int j = bpsize-1; j>i; j--)	{
					if (meanprob[permut[j]] > meanprob[permut[j-1]])	{
						int tmp = permut[j];
						permut[j] = permut[j-1];
						permut[j-1] = tmp;
					}
				}
			}
			if (! bench)	{
				osp << "bipartition\tmaxdiff\tmeanprob\tprobs\tlengths\n";
				osp << '\n';
			}
			for (int i=0; i<bpsize; i++)	{
				int k = permut[i];
				if ((! bench) || (meanprob[k] > conscutoff))	{
					mergedbplist->GetBipartition(k).WriteToStream(osp,verbose);
					// mergedbplist->GetBipartition(k).WriteToStream(osp);
					osp << '\t' << (int) (100 * diff[k]);
					osp << '\t' << (int) (100 * meanprob[k]);
					osp << '\t';
					for (int p=0; p<P; p++)	{
						osp << '\t' << (int) (100 * prob[p][k]);
					}
					osp << '\t';
					for (int p=0; p<P; p++)	{
						osp << '\t' << length[p][k];
					}
					osp << '\n';
				}
			}
			if (! bench)	{
				osp << '\n' << '\n';
			}
		}
		else	{
			int permut[bpsize];
			for (int k=0; k<bpsize; k++)	{
				permut[k] = k;
			}
			for (int i=0; i<bpsize; i++)	{
				for (int j = bpsize-1; j>i; j--)	{
					if (diff[permut[j]] > diff[permut[j-1]])	{
						int tmp = permut[j];
						permut[j] = permut[j-1];
						permut[j-1] = tmp;
					}
				}
			}
			if (! bench)	{
				osp << "bipartition\tmaxdiff\tprobs\n";
				osp << '\n';
			}
			for (int i=0; i<bpsize; i++)	{
				int k = permut[i];
				if ((! bench) || (meanprob[k] > conscutoff))	{
					mergedbplist->GetBipartition(k).WriteToStream(osp,verbose);
					// mergedbplist->GetBipartition(k).WriteToStream(osp);
					osp << '\t' << (int) (100 * diff[k]);
					osp << '\t';
					for (int p=0; p<P; p++)	{
						osp << '\t' << (int) (100 * prob[p][k]);
					}
					osp << '\t';
					for (int p=0; p<P; p++)	{
						osp << '\t' << length[p][k];
					}
					osp << '\n';
				}
			}
			if (! bench)	{
				osp << '\n' << '\n';
			}
		}


		if (verbose)	{
			cout << '\n';
			cout << "maxdiff     : " << maxdiff << '\n';
			cout << "meandiff    : " << meandiff << '\n';
			cout << '\n';
		}
	
		ofstream sos((OutFile + ".bpdiff").c_str());
		sos << '\n';
		sos << "maxdiff     : " << maxdiff << '\n';
		sos << "meandiff    : " << meandiff << '\n';
		sos << '\n';
	}
	else	{
		bplist[0]->WriteToStream(osp);
		if (verbose)	{
			cout << '\n';
		}
	}

	if (verbose)	{
		cout << "bipartition list in : " << OutFile <<  ".bplist\n";
	}

	if (!rootonly)	{
		Consensus* cons = new Consensus(mergedbplist,conscutoff);
		ofstream cons_os((OutFile + ".con.tre").c_str());
		cons->Phylip(cons_os, 1, 1, 1, 0);
		if (verbose)	{
			cout << "consensus in        : " << OutFile << ".con.tre\n";
		}
		if (ps)	{
			cerr << "in bipartition list: ps deprecated\n";
			exit(1);
		}
		delete cons;
	}
	if (verbose)	{
		cout << '\n';
	}

	delete mergedbplist;
	for (int p=0; p<P; p++)	{
		delete bplist[p]->mParam;
		delete bplist[p];
	}
	delete[] bplist;
	return maxdiff;
}

void BipartitionList::CompareWith(BipartitionList* with)	{

	mCompatibleArray = new int[mSize];
	for (int i=0; i<mSize; i++)	{
		mCompatibleArray[i] = mBipartitionArray[i]->CompareWith(with);
	}
}

void BipartitionList::SupportCheck(BipartitionList* with)	{

	if (! mCompatibleArray)	{
		mCompatibleArray = new int[mSize];
	}
	for (int i=0; i<mSize; i++)	{
		mCompatibleArray[i] = mBipartitionArray[i]->SupportCheck(with);
	}
}

void BipartitionList::BasicAllocation()	{

	mAllocatedSize = basicsize;
	mSize = 0;
	mWeight = 0;
	
	mCompatibleArray = 0;

	mWeightArray = new double[mAllocatedSize];
	for (int i=0; i<mAllocatedSize; i++)	{
		mWeightArray[i] = 0;
	}

	mLengthArray = new double[mAllocatedSize];
	for (int i=0; i<mAllocatedSize; i++)	{
		mLengthArray[i] = 0;
	}

	mBipartitionArray = new Bipartition*[mAllocatedSize];
	for (int i=0; i<mAllocatedSize; i++)	{
		mBipartitionArray[i] = 0;
	}
}

// ---------------------------------------------------------------------------------
//		 BipartitionList(TaxaParameters* inParam)
// ---------------------------------------------------------------------------------

BipartitionList::BipartitionList(TaxaParameters* inParam)	{

	mParam = inParam;
	BasicAllocation();
}

// ---------------------------------------------------------------------------------
//		 BipartitionList(string filename)
// ---------------------------------------------------------------------------------

BipartitionList::BipartitionList(string filename, int burnin, int every, int until, double cutoff, bool rootonly)	{

	BasicAllocation();
	
	if (burnin <0)	{
		if (burnin == -1)	{
			burnin = -5;
		}
		ifstream is(filename.c_str());
		int n = 0;
		string tmp;
		while (!is.eof())	{
			is >> tmp;
			n++;
		}
		burnin = - n / burnin;
	}
			
	ifstream is(filename.c_str());
	PBTree tree;	
	int size = 0;
	Ntree = 0;
	if (! tree.ReadFromStream(is))	{
		return;
	}
	size = 1;
	mParam = new TaxaParameters(&tree);	
	int cycle = 0;
	if (size > burnin)	{
		if (rootonly)	{
			Bipartition tempbp = tree.GetRootBipartition();
			Append2(tempbp,1,1);
		}
		else	{
			tree.Trichotomise();
			Prune(&tree);
		}
		Ntree++;
		cycle ++;
		if (cycle == every)	{
			cycle = 0;
		}
	}

	while (! is.eof())	{
		PBTree tree(mParam);
		int noteof = tree.ReadFromStream(is);
		if (noteof)	{
			size++;
			if (size > burnin)	{
				cycle ++;
				if (cycle == every)	{
					cycle = 0;
				}
				if (((until == -1) || (size<=until)) && (! cycle) && (tree.GetSize() == mParam->Ntaxa))	{
					tree.RegisterWithParam(mParam);
					if (rootonly)	{
						Bipartition tempbp = tree.GetRootBipartition();
						Append2(tempbp,1,1);
					}
					else	{
						BipartitionList tempList(mParam);
						tree.Trichotomise();
						tempList.Prune(&tree);
						Append(&tempList);
					}
					Ntree++;
				}
			}
		}
	}
	Sort();
	Truncate(cutoff);
	CheckForDuplicates();
}

// ---------------------------------------------------------------------------------
//		 BipartitionList( Tree* inTree)
// ---------------------------------------------------------------------------------


BipartitionList::BipartitionList(PBTree* inTree, double weight)	{

	mParam = inTree->GetParameters();
	BasicAllocation();
	Prune(inTree);
	// PruneWithSupports(inTree);
}


// ---------------------------------------------------------------------------------
//		 IsEqualTo
// ---------------------------------------------------------------------------------

int BipartitionList::IsEqualTo(BipartitionList* comp)	{
	int equal = 1;
	if (GetSize() != comp->GetSize())	{
		equal = 0;
		cerr << "comparing bp lists : non matching sizes: " << GetSize() << '\t' << comp->GetSize() << "\n";
	}
	for (int i=0; i<comp->GetSize(); i++)	{
		Bipartition bp = GetBipartition(i);
		int k = 0;
		while ((k<GetSize()) && (bp != comp->GetBipartition(k))) k++;
		if (k == GetSize())	{
			cerr << "comparing bplists: non matching bipartition\n";
			equal = 0;
		}
		else	{
			if (fabs(GetProb(i) - comp->GetProb(k)) > 1e-6)	{
				cerr << "comparing bp lists: non matching probs: " << GetProb(i) << '\t' << comp->GetProb(k) << '\n';
				equal = 0;
			}
			if (fabs(GetLength(i) -comp->GetLength(k)) > 1e-6)	{
				cerr << "comparing bp lists: non matching lengths: " << GetLength(i) << '\t' << comp->GetLength(k) << '\n';
				equal = 0;
			}
		}
	}
	return equal;
}
// ---------------------------------------------------------------------------------
//		 Prune(Tree* inTree)
// ---------------------------------------------------------------------------------


void BipartitionList::Prune(PBTree* inTree)	{
	
	if (inTree->IsDichotomous())	{
		cerr << "in BipartitionList::Prune : warning, tree is dichotomous\n";
	}	
	Flush();
	inTree->GetRoot()->BipartitionPruning(this);
	Modulo();
	mWeight = 1;
}


void BipartitionList::PruneWithSupports(PBTree* inTree)	{
	
	Flush();
	inTree->GetRoot()->BipartitionPruningWithSupports(this);
	Modulo();
	mWeight = 1;
}


// ---------------------------------------------------------------------------------
//		 BipartitionList(TreeList*)
// ---------------------------------------------------------------------------------


BipartitionList::BipartitionList( TreeList* inTreeList, double* weightarray, double cutoff)	{

	mParam = inTreeList->GetParameters();
	BasicAllocation();
	
	for (int i=0; i<inTreeList->GetSize(); i++)	{
		BipartitionList tempList(mParam);
		inTreeList->GetTree(i)->Trichotomise();
		tempList.Prune(inTreeList->GetTree(i));
		if (weightarray)	{
			tempList.Reweight(weightarray[i]);
		}
		Append(&tempList);
	}
	Sort();
	Truncate(cutoff);
	CheckForDuplicates();
}

// ---------------------------------------------------------------------------------
//		 ~BipartitionList()
// ---------------------------------------------------------------------------------

BipartitionList::~BipartitionList()	{

	for (int i=0; i<mAllocatedSize; i++)	{
		delete mBipartitionArray[i];
	}
	delete[] mBipartitionArray;
	delete[] mWeightArray;
	delete[] mLengthArray;
	delete[] mCompatibleArray;
}

		
// ---------------------------------------------------------------------------------
//		 Sort();		// reverse order !
// ---------------------------------------------------------------------------------

void BipartitionList::Sort()	{

	for (int i=0; i< mSize-1; i++)	{
		for (int j= mSize-1; i<j ; j--)		{
			if (mWeightArray[j] > mWeightArray[j-1])	{

				// swap the probs
				double tmp = mWeightArray[j];
				mWeightArray[j] = mWeightArray[j-1];
				mWeightArray[j-1] = tmp;

				//swap the mean lengths
				tmp = mLengthArray[j];
				mLengthArray[j] = mLengthArray[j-1];
				mLengthArray[j-1] = tmp;

				// swap the bipartitions
				Bipartition* Btemp = mBipartitionArray[j];
				mBipartitionArray[j] = mBipartitionArray[j-1];
				mBipartitionArray[j-1] = Btemp;
			}
		}
	}
}
	
// ---------------------------------------------------------------------------------
//		 GetHisto();
// ---------------------------------------------------------------------------------

	
void BipartitionList::GetHisto(int* histo, int ncat)	{
	for (int i=0; i<ncat; i++)	{
		histo[i] = 0;
	}
	for (int i=0; i<GetSize(); i++)	{
		int k = (int) (GetProb(i)*ncat);
		if (k == ncat) k--;
		histo[k] ++;
	}
}
			
// ---------------------------------------------------------------------------------
//		 RegisterWith();
// ---------------------------------------------------------------------------------

		
int BipartitionList::RegisterWith(TaxaParameters* inParam)	{

	int ok = 1;
	if (! inParam)	{
		cerr << "error : non existing taxon set\n";
		exit(1);
	}
	if (mParam->Ntaxa != inParam->Ntaxa)	{
		ok = 0;
	}
	if (ok)	{
	int* permut = new int[mParam->Ntaxa];
	int k = 0;
	while (k<mParam->Ntaxa)	{
		int j = 0;
		while ((j<mParam->Ntaxa) && (mParam->SpeciesNames[k] != inParam->SpeciesNames[j])) j++;
		if (j == mParam->Ntaxa)	{
			ok = 0;
		}
		permut[k] = j;
		k++;
	}

	for (int i=0; i<mSize; i++)	{
		mBipartitionArray[i]->PermutTaxa(permut);
	}	
	delete[] permut;
	}
	return ok;
}

// ---------------------------------------------------------------------------------
//		 Truncate(double cutoff);
// ---------------------------------------------------------------------------------

void BipartitionList::Truncate(double cutoff)	{
	
	int i=0;
	while ( (i<mSize) && (GetProb(i) >= cutoff))	{
		i++;
	}
	for (int j=i; j<mSize; j++)	{
		delete mBipartitionArray[j];
		mBipartitionArray[j] = 0;
	}
	mSize = i;
}


// ---------------------------------------------------------------------------------
//		 Flush()
// ---------------------------------------------------------------------------------

void BipartitionList::Flush(){
	mSize = 0;
	mWeight = 0;
}
	

// ---------------------------------------------------------------------------------
//		 Pop()
// ---------------------------------------------------------------------------------

void BipartitionList::Pop(){
	
	if (! mSize)	{
		cerr << "error in BipartitionList : pop called on list of size 0\n";
		exit(1);
	}
	mSize--;
}
	

// ---------------------------------------------------------------------------------
//		 Modulo
// ---------------------------------------------------------------------------------

void BipartitionList::Modulo()	{

	for (int i=0; i<mSize; i++)	{
		mBipartitionArray[i]->Modulo();
	}
}


// ---------------------------------------------------------------------------------
//		 CheckForDuplicates()
// ---------------------------------------------------------------------------------

int BipartitionList::CheckForDuplicates()	{
	int returnvalue = 0;
	for (int i=0; i<mSize; i++)	{
		for (int j=i+1; j<mSize; j++)	{
			if ((*mBipartitionArray[i]) == (*mBipartitionArray[j]))	{
				cerr << "error: found duplicates\n";
				mBipartitionArray[i]->WriteToStream(cerr);
				cerr << '\n';
				mBipartitionArray[j]->WriteToStream(cerr);
				cerr << '\n';
				returnvalue = 1;
			}
			if ((*mBipartitionArray[i]) == !(*mBipartitionArray[j]))	{
				cerr << "error: found symmetrical duplicates\n";
				mBipartitionArray[i]->WriteToStream(cerr);
				cerr << '\n';
				mBipartitionArray[j]->WriteToStream(cerr);
				cerr << '\n';
				returnvalue = 1;
			}
		}
	}
	return returnvalue;
}

// ---------------------------------------------------------------------------------
//		 IsCompatibleWith
// ---------------------------------------------------------------------------------

Boolean BipartitionList::IsCompatibleWith(const Bipartition& inPartition)	{

	Boolean temp = true;
	int i=0;
	while ((i<mSize) && temp)	{
		temp &= mBipartitionArray[i]->IsCompatibleWith(inPartition);
		i++;
	}
	return temp;
}

Boolean BipartitionList::IsCompatibleWith(BipartitionList* bplist){

	Boolean temp = true;
	for (int i=0; i<bplist->GetSize(); i++)	{
		temp &= IsCompatibleWith(*(bplist->mBipartitionArray[i]));
	}
	return temp;
}


void BipartitionList::Suppress(const Bipartition& bp)	{

	for (int i=0; i<mSize; i++)	{
		mBipartitionArray[i]->Suppress(bp);
	}
}
		

// ---------------------------------------------------------------------------------
//		 GetIndex()
// ---------------------------------------------------------------------------------

int BipartitionList::GetIndex(const Bipartition& inPartition){
	
	int i=0;
	while ( (i<mSize) && (*mBipartitionArray[i] != inPartition))	{
		i++;
	}
	if (i == mSize)	{
		i = -1;	
	}
	return i;
}


// ---------------------------------------------------------------------------------
//		 Reweight()
// ---------------------------------------------------------------------------------

void BipartitionList::Reweight(double weight)	{

	for (int i=0; i<mSize; i++)	{
		mWeightArray[i] *= weight / mWeight;
	}
	mWeight = weight;
}

// ---------------------------------------------------------------------------------
//		 Append()
// ---------------------------------------------------------------------------------

void BipartitionList::Append(BipartitionList* inList)	{

	for (int j=0; j<inList->GetSize(); j++)	{

		Bipartition& inPartition = (*inList)[j];
		Append(inPartition, inList->mWeightArray[j], inList->mLengthArray[j]);
	}
	mWeight += inList->GetWeight();
}

void BipartitionList::Append(Bipartition inPartition, double inWeight, double inLength)	{

	int i=0;
	while ( (i<mSize) && (*mBipartitionArray[i] != inPartition))	{
		i++;
	}
	if (i == mSize)	{
		mSize++;
		if (mSize > mAllocatedSize)	{
			Reallocate();
		}
		mBipartitionArray[i] = new Bipartition(inPartition);
		mBipartitionArray[i]->mParam = mParam;
		mWeightArray[i] = 0;
		mLengthArray[i] = 0;
	}
	mWeightArray[i] += inWeight;
	mLengthArray[i] += inLength;
}

void BipartitionList::Append2(Bipartition inPartition, double inWeight, double inLength)	{

	int i=0;
	while ( (i<mSize) && (*mBipartitionArray[i] != inPartition))	{
		i++;
	}
	if (i == mSize)	{
		mSize++;
		if (mSize > mAllocatedSize)	{
			Reallocate();
		}
		mBipartitionArray[i] = new Bipartition(inPartition);
		mBipartitionArray[i]->mParam = mParam;
		mWeightArray[i] = 0;
		mLengthArray[i] = 0;
	}
	mWeightArray[i] += inWeight;
	mLengthArray[i] += inLength;
	mWeight += inWeight;
}



// ---------------------------------------------------------------------------------
//		 FastInsert()
// ---------------------------------------------------------------------------------
// do not check if inbipartition already present in the list
// (assumes it is absent)

void BipartitionList::FastInsert(Bipartition inPartition, double weight, double length)	{

	mSize++;
	if (mSize > mAllocatedSize)	{
		Reallocate();
	}
	mBipartitionArray[mSize-1] = new Bipartition(inPartition);
	mWeightArray[mSize-1] = weight;
	mLengthArray[mSize-1] = length;

}


// ---------------------------------------------------------------------------------
//		 Insert(Bipartition inPartition, double weight, double length)
// ---------------------------------------------------------------------------------
// check whether the bipartition is in the list
// if yes, send an error message

void BipartitionList::Insert(Bipartition inPartition, double weight, double length)	{

	if (CheckLevel)	{
		int i=0;
		while ( (i<mSize) && (*mBipartitionArray[i] != inPartition))	{
			i++;
		}
		if (i != mSize)	{
			cerr << "error : inserting an already existing bipartition\n";
			exit(1);
		}
	}
	
	mSize++;
	if (mSize > mAllocatedSize)	{
		Reallocate();
	}
	mBipartitionArray[mSize-1] = new Bipartition(inPartition);
	mWeightArray[mSize-1] = weight;
	mLengthArray[mSize-1] = length;

}


// ---------------------------------------------------------------------------------
//		 Reallocate()
// ---------------------------------------------------------------------------------

void BipartitionList::Reallocate()	{


	double* mWeightArray2 = new double[2*mAllocatedSize];
	for (int i=0; i<mAllocatedSize; i++)	{
		mWeightArray2[i] = mWeightArray[i];
	}
	for (int i=mAllocatedSize; i<2*mAllocatedSize; i++)	{
		mWeightArray2[i] = 0;
	}

	double* mLengthArray2 = new double[2*mAllocatedSize];
	for (int i=0; i<mAllocatedSize; i++)	{
		mLengthArray2[i] = mLengthArray[i];
	}
	for (int i=mAllocatedSize; i<2*mAllocatedSize; i++)	{
		mLengthArray2[i] = 0;
	}

	Bipartition** mBipartitionArray2 = new Bipartition*[2*mAllocatedSize];
	for (int i=0; i<mAllocatedSize; i++)	{
		mBipartitionArray2[i] = mBipartitionArray[i];
	}
	for (int i=mAllocatedSize; i<2*mAllocatedSize; i++)	{
		mBipartitionArray2[i] = 0;
	}

	delete[] mWeightArray;
	mWeightArray = mWeightArray2;
	delete[] mLengthArray;
	mLengthArray = mLengthArray2;
	delete[] mBipartitionArray;
	mBipartitionArray = mBipartitionArray2;
	mAllocatedSize *= 2;
}


// ---------------------------------------------------------------------------------
//		 WriteToStream(ostream& os, int header)
// ---------------------------------------------------------------------------------

void BipartitionList::WriteToStream(ostream& os, int header, int verbose)	{

	if (header)	{
		mParam->WriteToStream(os);
	}
	os << "BipartitionList\n";	
	os << "Size " << mSize << '\n';
	os << "Weight " << mWeight << '\n';
	os << '\n';
 	for (int i=0; i<mSize; i++)	{
		mBipartitionArray[i]->WriteToStream(os,verbose);
		os << '\t';
	       	if (GetProb(i) == 2)	{
			cerr << "error in bplist: crisse d'ostie de programme!\n";
			exit(1);
			os << 1;
		}
		else	{
			os << GetProb(i);
		}
		os << '\t';
		os << mBipartitionArray[i]->GetPriorProb() << '\t';
		os << GetLength(i) << '\n';
	}
}


// ---------------------------------------------------------------------------------
//		 ReadFromStream()
// ---------------------------------------------------------------------------------

void BipartitionList::ReadFromStream(istream& is)	{

	string temp;
	is >> temp;
	while (temp != "End")	{
		if (temp == "TaxaList")	{
			mParam = new TaxaParameters();
			mParam->ReadFromStream(is);
		}
		else if (temp == "BipartitionList")	{
		
			int size;
			is >> temp >> size >> temp >> mWeight;
			cerr << size << '\n';
			mSize = 0;
				
			for (int i=0; i<size; i++)	{
				double prob;
				double length;
				Bipartition bp(mParam);
				bp.ReadFromStream(is);
				is >> prob >> length ;
				Insert(bp,mWeight*prob,mWeight*length);
				// FastInsert(bp,mWeight*prob,mWeight*length);
			}
		}
		else	{
			cerr << "error in BipartitionList::ReadFromStream\n";
			exit(1);
		}
		is >> temp;
	}
}

