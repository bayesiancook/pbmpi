
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
#include <tr1/unordered_map>
#include <fstream>
#include <sstream>



// A TaxaParametersBis instance contain a map from TaxaName to an index
// Index are continuous between 0 and GetNtaxa()
class TaxaParametersBis {

	int Ntaxa;

	tr1::unordered_map<string,int> map;

	// Recursive function used by constructor
	void RecursiveAddSpeciesName(Link* from, int& i){
		if (from->isLeaf()){
			map[from->GetNode()->GetName()] = i;
			i++;
		}
		else{
			for (const Link* link=from->Next(); link!=from; link=link->Next()){
				RecursiveAddSpeciesName(link->Out(), i);
			}
		}
	}

	public:

	// Constructor, take the filename of a newick tree.
	// Tree should have uniq taxa names at each leaf.
	TaxaParametersBis(string filename){
		Tree* tree = new Tree(filename);
		Ntaxa = tree->GetSize();
		int i = 0;
		RecursiveAddSpeciesName(tree->GetRoot(), i);
		delete tree;
	}
	
	// Use this method to have the number of taxa
	int GetNtaxa(){
		return Ntaxa;
	}

	// Use this method to have the index corresponding to the taxa name
	// Mostly for tests
	int GetSpeciesIndex(string taxa){
		if(map.find(taxa) == map.end()){
			cerr << taxa << " is not present in TaxaParametersBis\n";
			exit(1);
		}
		else{
			return map[taxa];
		}
	}

	string GetIndexToSpecies(int sp){
		typedef tr1::unordered_map<string,int>::iterator LocalIt;
		for (LocalIt local_it = map.begin(); local_it!= map.end(); ++local_it ){
			if(local_it->second == sp){
				return local_it->first;
			}
		}
	}

};
