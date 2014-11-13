
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "BP2Stat.h"

// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------
// -------------------- BipartitionList -----------------------
// ------------------------------------------------------------
// ------------------------------------------------------------



BipartitionList::BipartitionList(string infilename, TaxaParametersBis* intaxaparameters, int burnin, int every, int until, int inverbose, double  inmixstatdown, double inmixstatup, double inCutOff) 
				: BipartitionContainer(intaxaparameters){
	
	filename = infilename;
	mixstatdown = inmixstatdown;
	mixstatup = inmixstatup;
	CutOff = inCutOff;

	// Opening the file
	if (! ifstream(filename.c_str())){
		string name = filename + ".treelist";
		if (! ifstream(name.c_str())){
			cerr << "Warning : cannot find " << filename << " nor " << name << "\n";
			exit(1);
		}
		filename = name;
	}
	ifstream is(filename.c_str());


	// Reading the forest
	string temp;
	int maxsize = (until-burnin)/every;
	for(int i = 0; i<burnin;i++){
		is >> temp;
	}
	int size = 0;
	while (is >> temp and ( (until == -1) or (size < maxsize) ) ){
		for(int i = 1; i<every;i++){
			is >> temp;
		}
		istringstream iss(temp);
		Tree* tree = new Tree(iss); 
		Trichotomise(tree->GetRoot());
		// Reading a tree
		RecursiveCreateList(tree->GetRoot(), size);
		size++;
		delete tree;
	}
	cerr << filename << " : " << size << " trees\n";

	// Normalization
	LocalIt local_it = map.begin();
	while(local_it != map.end()){
		local_it->second->Normalize(size);
		if(local_it->second->GetBootstrap() < CutOff){
			local_it = map.erase(local_it);
		}
		else{
			++local_it;
		}
	}

	// Create an array of bipartition sorted according to their bootstrap
	sortedarray = new ShellSort<Statistics>(map, &Statistics::GetBootstrap);

}


void BipartitionList::Write(BipartitionHashTable<CrossStatistics*>* otherMap){
	// Output
	ofstream os((filename + ".bplist2").c_str());
	os.precision(6);
	int totalswitch=0;
	double sumMixStat = 0;
	double minMixStat = 10;
	int NbMixStat = 0;
	for(int index=0; index<sortedarray->GetSize(); index++){
		string bp=sortedarray->Get(index);
		if( CountTaxaMin(bp) > 1){
			double bootstrap = map[bp]->GetBootstrap();
			os << fixed << BipartitionToString(bp) << '\t';
			os << bootstrap << '\t';
			if	( (otherMap != 0) and (otherMap->find(bp) != otherMap->end() ))	{
				os << otherMap[0][bp]->GetBootstrap()<< '\t';
				if((bootstrap < mixstatup) and (bootstrap > mixstatdown)){
					double ms = map[bp]->GetMixStat( otherMap[0][bp]->GetBootstrap() );
					sumMixStat += ms;
					if(ms<minMixStat){minMixStat=ms;}
					os << ms << '\t';
					NbMixStat++;
				}
				else{
					os << "    -     \t";
				}
			}
			else{
				os << "< cut off\t    -     \t";
			}
			totalswitch += map[bp]->GetNbSwitch();
			os << map[bp]->GetNbSwitch() << '\t';
			os << map[bp]->GetHistory() << '\n';
		}
	}
	os << "Bipartitions change :\t" << totalswitch << '\n';
	os << "Bipartitions number :\t" << map.size() << '\n';
	os << "Mean of MixStatistic :\t" << sumMixStat/NbMixStat << '\n';
	os << "Mini of MixStatistic :\t" << minMixStat << '\n';
	cout << "\nmin of MixStatistic :\t" << minMixStat << '\n';
	cout <<'\t'<< filename <<".bplist2 written\n";
	os.close();

	if(verbose){
		ofstream os((filename + ".details").c_str());
		for(int index=0; index<sortedarray->GetSize(); index++){
			string bp=sortedarray->Get(index);
			double bootstrap = map[bp]->GetBootstrap();
			if((bootstrap < 0.8)and (bootstrap > 0.2)){
				os << fixed << BipartitionToString(bp) << '\t';
				os << bootstrap << '\t';
				os << BipartitionToTaxaList(bp) << '\n';
			}
		}
		cout <<'\t'<< filename <<".details written\n";
		os.close();
	}

}

Statistics* BipartitionList::GetStatistics(string bp){
	if(map.find(bp) == map.end()){
		return 0;
	}
	else{
		return map[bp];
	}
}
	
void BipartitionList::NewBipartition(string bp, double length, int size){
	if(map.find(bp) == map.end()){
		map[bp] = new Statistics(size, length);
	}
	else{
		map[bp]->Add(size, length);
	}
}
	

string BipartitionList::RecursiveCreateList(Link* from, int size){
	if (from->isLeaf()){
		string bp = CreateOneTaxaBipartition(from->GetNode()->GetName());
		NewBipartition(bp, atof(from->GetBranch()->GetName().c_str()), size);
		return bp;
	}
	else{
		string bp = CreateEmptyBipartition();
		for (const Link* link=from->Next(); link!=from; link=link->Next()){
			string bp2 = RecursiveCreateList(link->Out(),size);
			bp = MergeBipartition(bp, bp2);
		}
		if(! from->isRoot()){
			NewBipartition(bp, atof(from->GetBranch()->GetName().c_str()), size);
		}
		return bp;
	}
}




// Unroot the tree when it is rooted
void BipartitionList::Trichotomise(Link* root){
	if(root->Next()->Next()->Next() == root){
		cerr << " BipartitionList::Trichotomise was not tested ! \n";
		if (root->Next()->Out()->isLeaf()){
			root->Knit();
		}
		if (root->Next()->Out()->isLeaf()){
			cerr << "error : cannot trichotomise a two Taxa tree\n";
			exit(1);
		}
		Link* next = root->Next();
		next->Insert(root->Next()->Out()->Next()->Next());
		next->Insert(root->Next()->Out()->Next());
		root->SetNext(root->Next()->Next());
		delete next->Out()->GetNode();
		delete next->Out();
		delete next->GetBranch();
		delete next;
	}
}


// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------
// -------------------- ConsensusTree -------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------


ConsensusTree::ConsensusTree(TaxaParametersBis* intaxaparameters,
				BipartitionHashTable<CrossStatistics*> inmap)
					: BipartitionContainer(intaxaparameters){

	map = inmap;

	bipartitionlist = new ShellSort<CrossStatistics>(map, &CrossStatistics::GetBootstrap);



	//MultifurcateTree
	root = new Link();
	root->SetNode(new Node());
	for(int i = 0; i < GetNtaxa(); i++){
		string taxaname = taxaparameters->GetIndexToTaxa(i);
		Node* node = new Node(taxaname);
		Link* leaf = new Link;
		Link* leafout = new Link;
		Branch* branch = new Branch();

		branch->SetName( MakeBranchName(1,map[CreateOneTaxaBipartition(taxaname)]->GetLength()) );

		leaf->SetNode(node);
		leaf->SetBranch(branch);
		leaf->InsertOut(leafout);
		leafout->SetBranch(branch);
		root->Insert(leafout);
	}



	//Compatible N-3 partitions
	string* compatibleBipartitionList = new string[GetNtaxa() - 3];
	int it_from=0;
	for(int it_to=0; it_to < GetNtaxa()-3;it_from++){
		if(it_from == bipartitionlist->GetSize()) {
			cerr << "bipartitionlist is not big enough to construct ConsensusTree\n";
			exit(1);
		}
		bool ok = true;
		if( CountTaxaMin(bipartitionlist->Get(it_from)) > 1){ // If there is more than one taxa
			for(int i = 0; i < it_to ; i++){
				ok = ok and (isCompatible(compatibleBipartitionList[i], bipartitionlist->Get(it_from)));
			}
		}
		else{ok=false;}

		if(ok){compatibleBipartitionList[it_to]=bipartitionlist->Get(it_from);it_to++;}
	}

	//Adding the branches
	for(int j = 0; j < GetNtaxa()-3 ; j++){
		Link* link = FindLinkToBifurcate(GetRoot(), compatibleBipartitionList[j]);
		if(!link){cerr << "problem constructing ConsensusTree\n";exit(1);}
		CreateBifurcation(link, compatibleBipartitionList[j]);
	}
		
}


string ConsensusTree::MakeBranchName(double bootstrap, double length){
	ostringstream s;
	int b = bootstrap*100;
	if(b!=100){s.precision(0);s << "0." << b;}
	s.precision(6);
	s << ':' << length;
	return s.str();
}

void ConsensusTree::CreateBifurcation(Link* from, string bp){
	Link* bip = new Link();
	bip->SetNode(new Node());
	Link* bipout = new Link;
	bip->InsertOut(bipout);
	Branch* bipbranch = new Branch();
	bipbranch->SetName( MakeBranchName(map[bp]->GetBootstrap(), map[bp]->GetLength()) );
	bipout->SetBranch(bipbranch);
	bip->SetBranch(bipbranch);
	bipout->SetNode(from->GetNode());

	Link* previous=from;
	for (Link* link=from->Next(); link!=from; link=previous->Next()){
		int type = TaxaType(link->Out(), bp);
		if(type == 2){
			cerr << "CreateBranch should not be call here\n";
			exit(1);
		}
		if(type){
			previous->SetNext(link->Next());
			bip->Insert(link);
		}
		else{
			previous=link;
		}
	}

	if(bip->Next()->Next() == bip){
		cerr <<"bip unaire\n";
		from->Insert(bip->Next());
		bip->SetNext(bip);
		for (Link* link=from->Next()->Next(); link!=from; link=from->Next()){
			from->Next()->SetNext(link->Next());
			bip->Insert(link);
		}
	}

	from->Insert(bipout);
}


// From a given link, explore all the subtree.
// Return 1 if all the taxa are present in the bipartition, 0 if not one is present
// Return 2 if at least one taxa is present and one is not
int ConsensusTree::TaxaType(Link* from, string bp){
	if(from->isLeaf()){
		return TaxaIsHere(bp, from->GetNode()->GetName());
	}
	else{
		int ret = TaxaType(from->Next()->Out(), bp);
		for (const Link* link=from->Next()->Next(); link!=from; link=link->Next()){
			if(ret != 2){
				int t = TaxaType(link->Out(), bp);
				if(t != ret){
					ret = 2;
				}
			}
		}
		return ret;
	}
}

Link* ConsensusTree::FindLinkToBifurcate(Link* from, string bp){
	Link* ret = 0;
	const Link* link=from->Next();
	while(!ret and link!=from){
		ret = FindLinkToBifurcate(link->Out(), bp);
		link=link->Next();
	}
	if(ret){return ret;}
	else{
		if( TaxaType(from, bp) == 2 ){
			return from;
		}
	}
	return 0;
}
