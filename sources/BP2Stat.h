
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Tree.h"
#include "BP2util.h"
#include <tr1/unordered_map>
#include <fstream>
#include <sstream>
#include <cmath>


// In order to read a file were a lot of trees with the same taxa set
class BipartitionList: public BipartitionContainer {

	int verbose;
	double CutOff;	
	double mixstatdown;
	double mixstatup;

	// hashtable from bipartition to statistic
	BipartitionHashTable<Statistics*> map;
	typedef BipartitionHashTable<Statistics*>::iterator LocalIt;

	// Array of sorted bipartitions
	StringArray* sortedarray;

	string filename;
	

	public:

	// Constructor take the name of a file having a set of trees containing the same set of taxa
	BipartitionList(string infilename, TaxaParametersBis* intaxaparameters, int burnin, int every, int until, int verbose, double  inmixstatdown, double inmixstatup, double inCutOff =0);

	void Write(BipartitionHashTable<CrossStatistics*>* map = 0);
	~BipartitionList(){
		delete[] sortedarray;
	}


	double GetBootstrap(string bp){
		Statistics* s = GetStatistics(bp);
		if(s==0){
			return 0;
		}
		else{
			return s->GetBootstrap();
		}
	}

	double GetLength(string bp){
		Statistics* s = GetStatistics(bp);
		if(s==0){
			return 0;
		}
		else{
			return s->GetLength();
		}
	}

	StringArray* GetSortedArray(){
		return sortedarray;
	}

	// Return a pointer to an instance of class Statistics
	Statistics* GetStatistics(string bp);

	private :


	void NewBipartition(string bp, double length, int size);

	string RecursiveCreateList(Link* from, int size);


	// Unroot the tree when it is rooted
	// Not used 
	void Trichotomise(Link* root);

};

class ConsensusTree : public BipartitionContainer, public Tree{


	BipartitionHashTable<CrossStatistics*> map;
	StringArray* bipartitionlist;

	public : 

	ConsensusTree(TaxaParametersBis* intaxaparameters, BipartitionHashTable<CrossStatistics*> inmap);

	private :
	

	string MakeBranchName(double bootstrap, double length);

	void CreateBifurcation(Link* from, string bp);


	// From a given link, explore all the subtree.
	// Return 1 if all the taxa are present in the bipartition, 0 if not one is present
	// Return 2 if at least one taxa is present and one is not
	int TaxaType(Link* from, string bp);

	Link* FindLinkToBifurcate(Link* from, string bp);

};



class BipartitionCompare : public BipartitionContainer{
	
	string outFile;


	double CutOff;
	int NbBipartitionList;
	BipartitionList** list;
	StringArray* sortedarray;

	BipartitionHashTable<CrossStatistics*> map;
	typedef BipartitionHashTable<CrossStatistics*>::iterator LocalIt;

	public :

	// Constructor take the name of a file having a set of trees containing the same set of taxa
	BipartitionCompare(BipartitionList** inlist, int inNbBipartitionList, double inCutOff, string inoutFile) : BipartitionContainer(inlist[0]->GetTaxaParameters()){
		list=inlist;
		NbBipartitionList = inNbBipartitionList;
		CutOff=inCutOff;
		outFile=inoutFile;

		// Run over the list
		for(int i=0; i<NbBipartitionList; i++){
			for( int j=0; j < list[i]->GetSortedArray()->GetSize() ; j++ ){
				string bp = list[i]->GetSortedArray()->Get(j);
				if(map.find(bp) == map.end()){
					map[bp] = new CrossStatistics();
				}
				map[bp]->Add(list[i]->GetStatistics(bp));
			}
		}

		// Normalize and remove when bootstrap is under the CutOff
		LocalIt local_it = map.begin();
		while(local_it != map.end()){
			local_it->second->Normalize(NbBipartitionList);
			cerr << local_it->second->GetBootstrap() << ' ';
			if(local_it->second->GetBootstrap() < CutOff){
				local_it = map.erase(local_it);
			}
			else{
				++local_it;
			}
		}


		for(int i=0; i<NbBipartitionList; i++){
			list[i]->Write(&map);
		}

		// Create an array of bipartition sorted according to their bootstrap
		sortedarray = new ShellSort<CrossStatistics>(map, &CrossStatistics::GetMaxDiff);
		if(NbBipartitionList>1){
			cout <<"\n\n\tmaxdiff\t:"<< map[sortedarray->Get(0)]->GetMaxDiff() <<'\n';
		}
		Write();
	}

	void Write(){
		// OutputStringArray* bipartitionlist
		if(outFile==""){outFile="bpcomp2";}
		ofstream osbplist((outFile+".bplist").c_str());
		osbplist.precision(6);
		osbplist << "total number of bp : " << sortedarray->GetSize() << '\n';
		osbplist << "cutoff : " << CutOff << '\n';
		osbplist << "Ntaxa "<< GetNtaxa() << "\nTaxaList :\n\n";
		for(int i = 0; i<GetNtaxa();i++){
			osbplist << GetTaxaParameters()->GetIndexToTaxa(i) << '\n';
		}
		osbplist <<'\n';
		osbplist << "bipartition" << string(GetNtaxa()-11, ' ');
		osbplist << "\tMaxDiff\t\tBootstrapVar\tBootstrapMean\t";
		for(int i=0; i<NbBipartitionList; i++){
			osbplist << "    file " << i << '\t';
		}
		osbplist <<'\n';


		for(int index=0; index < sortedarray->GetSize(); index++){
			string bp = sortedarray->Get(index);
			if( CountTaxaMin(bp) > 1){
				osbplist << BipartitionToString(bp);
				osbplist << fixed << '\t' << map[bp]->GetMaxDiff();
				osbplist << '\t' << map[bp]->GetVarBootstrap();
				osbplist << '\t' << map[bp]->GetBootstrap();
				for(int i=0; i<NbBipartitionList; i++){
					osbplist << '\t' << list[i]->GetBootstrap(bp);
				}
				osbplist << '\n';
			}
		}
		cout << '\t'<<outFile<<".bplist written\n";
		osbplist.close();

		Tree* t = new ConsensusTree(taxaparameters,map);
		ofstream osct((outFile+".con.tre").c_str());
		t->Print(osct);
		osct.close();
		cout << '\t'<<outFile<<".con.tre written\n";
	}

};

