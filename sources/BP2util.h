
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
#include <fstream>
#include <sstream>
#include <cmath>


#include <tr1/unordered_map>
template<class T>
class BipartitionHashTable : public tr1::unordered_map<string,T> {};



class  StringArray {

	int size;

	protected :

	string* stringarray;

	public :

	StringArray(int insize){
		size=insize;
		stringarray = new string[GetSize()];
	}

	~StringArray(){
		delete[] stringarray;
	}

	string Get(int i){
		if(i>GetSize()){
			cerr << "ShellSort::GetElement called with " << i << " when there is only " << GetSize() <<" elements\n";
			exit(1);
		}
		return stringarray[i];
	}


	int GetSize(){
		return size;
	}

};



template<class T>
class ShellSort : public StringArray {

	public :

	ShellSort( BipartitionHashTable<T*> mymap, double (T::*ptr)(), bool decreasing=true) : StringArray(mymap.size()){

		// Copy of all the keys in the array
		int i = 0;
		typename BipartitionHashTable<T*>::iterator local_it = mymap.begin();
		while(local_it != mymap.end()){
			stringarray[i]=local_it->first;
			i++;
			++local_it;
		}

		// Sorting
		for(int gap = GetSize()/2; gap!=0; gap/=2){
			for(int index=gap; index<GetSize(); index++){
				for(int index2 = index; index2 >= gap ;index2-=gap){
					if( decreasing ^ (mymap[stringarray[index2]]->*ptr)() < (mymap[stringarray[index2-gap]]->*ptr)()){
						string swap = stringarray[index2];
						stringarray[index2] = stringarray[index2-gap];
						stringarray[index2-gap] = swap;
					}
					else{break;}
				}
			}
		}
	}

};


// A TaxaParametersBis instance contain a map from TaxaName to an index
// Index are continuous between 0 and GetNtaxa()
class TaxaParametersBis {

	int Ntaxa;

	BipartitionHashTable<int> map;

	// Recursive function used by constructor
	void RecursiveAddTaxaName(Link* from, int& i);

	public:

	// Constructor, take the filename of a newick tree.
	// Tree should have uniq taxa names at each leaf.
	TaxaParametersBis(string filename);

	// Use this method to have the number of taxa
	int GetNtaxa(){
		return Ntaxa;
	}

	// Use this method to have the index corresponding to the taxa name
	// Mostly for tests
	int GetTaxaIndex(string taxa){
		if(map.find(taxa) == map.end()){
			cerr << taxa << " is not present in TaxaParametersBis\n";
			exit(1);
		}
		else{
			return map[taxa];
		}
	}

	string GetIndexToTaxa(int sp){
		typedef BipartitionHashTable<int>::iterator LocalIt;
		for (LocalIt local_it = map.begin(); local_it!= map.end(); ++local_it ){
			if(local_it->second == sp){
				return local_it->first;
			}
		}
	}

};



// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------
// ---- Those function are specific to string/bipartition -----
// ------------------------------------------------------------
// ------------------------------------------------------------

class BipartitionContainer {

	protected :

	int Nchar;
	TaxaParametersBis* taxaparameters;

	public :

	BipartitionContainer(TaxaParametersBis* intaxaparameters){
		taxaparameters = intaxaparameters;
		if(!taxaparameters){
			cerr << "TaxaParametersBis not define\n"; exit(1);
		}
		Nchar=GetNtaxa()/8+1;
	}

	TaxaParametersBis* GetTaxaParameters(){
		return taxaparameters;
	}

	int GetNtaxa(){
		return taxaparameters->GetNtaxa();
	}

	protected:

	void CheckBipartition(string bp);

	// Create a bipartition containing/excluding everything
	string CreateEmptyBipartition();

	// Change one taxa of bipartition
	string SwitchTaxa(string bp, int index);


	// Change one taxa of bipartition
	bool isCompatible(string a, string b);

	// Create a bipartition containing/excluding everything except one taxa
	string CreateOneTaxaBipartition(string Taxa);

	// Union of 2 bipartition were there is no overlap
	string MergeBipartition(string a, string b);

	string ReverseBipartition(string bp);

	bool TaxaIsHere(string bp, int index);
	bool TaxaIsHere(string bp, string taxa){
		return TaxaIsHere(bp, taxaparameters->GetTaxaIndex(taxa));
	}

	int CountTaxaMin(string bp);

	string BipartitionToString(string bp);
	string BipartitionToTaxaList(string bp);

};



// A TaxaParametersBis Statistic aggregate the presence of something (here bipartitions)
// Add(index) have to be call for each index were the bipartition is present.
// Normalize(lastindex) have to be call at the end
class Statistics {

	int lastPresence;
	double bootstrap;
	double NbSwitch;
	double mixstat;
	string history;

	double length;

	public :

	Statistics(int index, double inlength){
		NbSwitch = (index!=0);
		history = string(index,'.');
		history += 'X';
		bootstrap = 1;
		length = inlength;
		lastPresence = index;
	}


	void Add(int index, double inlength){
		if(index != lastPresence+1){
			NbSwitch+=2;
			history += string(index-lastPresence-1,'.');
		}
		history += 'X';
		bootstrap++;
		lastPresence = index;
		length +=inlength;
	}

	void Normalize(int index){
		if(lastPresence != index){
			NbSwitch++;
			history += string(index-lastPresence,'.');
		}
		length/=bootstrap;
		bootstrap/=index;
		lastPresence=index; //Okay now it s the total size ...

	}

	// Return the bootstrap when Normalize have been call
	double GetBootstrap(){
		return bootstrap;
	}

	string GetHistory(){
		return history;
	}


	double GetMixStat(double globalbootstrap){
		if(globalbootstrap==-1){
			return (NbSwitch)/(2 * globalbootstrap * (1-globalbootstrap) * lastPresence);
		}
		else{
			return (NbSwitch)/(2 * bootstrap * (1-bootstrap) * lastPresence);
		}
	}


	double GetLength(){
		return length;
	}

	int GetNbSwitch(){
		return NbSwitch;
	}

};


// A TaxaParametersBis Statistic to aggregate some Statistics*
class CrossStatistics {

	double meanBootstrap;
	double varBootstrap;
	double length;
	double min;
	double max;
	int nb;

	public :

	CrossStatistics(){
		meanBootstrap = 0;
		varBootstrap = 0;
		max = 0;
		min = 1;
		length = 0;
		nb = 0;
	}


	void Add(Statistics* s){
		length += s->GetLength();
		double b = s->GetBootstrap();
		meanBootstrap += b;
		varBootstrap += b*b;
		if(min>b){min=b;}
		if(max<b){max=b;}
		nb++;
	}

	double GetBootstrap(){
		return meanBootstrap;
	}

	double GetVarBootstrap(){
		return varBootstrap;
	}

	double GetMaxDiff(){
		return max-min;
	}

	double GetLength(){
		return length;
	}

	void Normalize(int size){
		if(size>nb){
			min=0;
		}
		length/=nb;
		nb=size;
		double sum = meanBootstrap;
		meanBootstrap/=nb;
		varBootstrap = sqrt(meanBootstrap*(meanBootstrap*nb-2*sum) + varBootstrap)/nb;
	}


};




