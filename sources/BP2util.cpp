

/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "BP2util.h"



// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------
// -------------------- TaxaParameeters ----------------------
// ------------------------------------------------------------
// ------------------------------------------------------------

TaxaParametersBis::TaxaParametersBis(string filename){
	Tree* tree = new Tree(filename);
	Ntaxa = tree->GetSize();
	int i = 0;
	RecursiveAddTaxaName(tree->GetRoot(), i);
	delete tree;
}


void TaxaParametersBis::RecursiveAddTaxaName(Link* from, int& i){
	if (from->isLeaf()){
		map[from->GetNode()->GetName()] = i;
		i++;
	}
	else{
		for (const Link* link=from->Next(); link!=from; link=link->Next()){
			RecursiveAddTaxaName(link->Out(), i);
		}
	}
}



// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------
// ---- Those function are specific to string/bipartition -----
// ------------------------------------------------------------
// ------------------------------------------------------------



void BipartitionContainer::CheckBipartition(string bp){
	for(int i = GetNtaxa(); i<Nchar*8; i++){
		if(TaxaIsHere(bp,i)){
			cerr << "bit " << i << " should be off\n";
			cerr << "with " << GetNtaxa() << " taxa\n";
			exit(1);
		}
	}
}



string BipartitionContainer::SwitchTaxa(string bp, int index){
	if(!index){
		bp = ReverseBipartition(bp);
	}
	bp[index/8] ^= 1 <<  (index%8);
	return bp;
}



string BipartitionContainer::CreateEmptyBipartition(){
	return string(Nchar, '\0');	
}



string BipartitionContainer::CreateOneTaxaBipartition(string Taxa){
	string bp = SwitchTaxa(CreateEmptyBipartition(), taxaparameters->GetTaxaIndex(Taxa));
	return bp;
}



string BipartitionContainer::MergeBipartition(string a, string b){
	CheckBipartition(a);
	CheckBipartition(b);
	for(int i =0;i< Nchar;i++){
		a[i] = a.c_str()[i] ^ b.c_str()[i];
	}
	return a;
}



bool BipartitionContainer::isCompatible(string a, string b){
	string m = MergeBipartition(a,b);
	string m2 = CreateEmptyBipartition();
	for(int i =0;i< Nchar;i++){
		m2[i] = m.c_str()[i] & a.c_str()[i];
	}
	bool ok = ( (m2 == a) or (m2 == m) or (m2 == CreateEmptyBipartition()));
	return ok ;		
}



string BipartitionContainer::ReverseBipartition(string bp){
	CheckBipartition(bp);
	for(int i =0;i< Nchar;i++){
		bp[i] = ~bp[i];
	}
	for(int i = GetNtaxa(); i<Nchar*8; i++){
		bp = SwitchTaxa(bp, i);
	}
	return bp;
}



bool BipartitionContainer::TaxaIsHere(string bp, int index){
	return (bp[index/8] >> (index%8)) & 1;
}



string BipartitionContainer::BipartitionToString(string bp){
	CheckBipartition(bp);
	string ret = "";
	for(int i =0;i< GetNtaxa();i++)
		ret += (TaxaIsHere(bp, i) ? '*' : '.');
	return ret;
}

int BipartitionContainer::CountTaxaMin(string bp){
	CheckBipartition(bp);
	//Count the number of taxa setted to one !
	int count=0;
	for (unsigned char mask = 0x80; mask != 0; mask >>= 1)
		for(int i =0;i< Nchar;i++)
  			if(bp[i] & mask) count++;

	// Decide to print the 1 or 0 taxas
	if(count>(GetNtaxa()/2)){return GetNtaxa() - count;}
	else{return count;}
}



string BipartitionContainer::BipartitionToTaxaList(string bp){
	CheckBipartition(bp);
	//Count the number of taxa setted to one !
	int count=0;
	for (unsigned char mask = 0x80; mask != 0; mask >>= 1)
		for(int i =0;i< Nchar;i++)
  			if(bp[i] & mask) count++;

	// Decide to print the 1 or 0 taxas
	int reverse=0;
	if(count>(GetNtaxa()/2)){reverse=1;}

	string ret = "";
	for(int i =0;i< GetNtaxa();i++){
		if( (TaxaIsHere(bp, i) xor reverse ) ){
			ret += taxaparameters->GetIndexToTaxa(i);
			ret += ' ';
		}
	}
	return ret;
}
