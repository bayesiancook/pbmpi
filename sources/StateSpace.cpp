
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <iostream>
#include <cstdlib>
using namespace std;

#include <sstream>
#include "BiologicalSequences.h"
#include "StateSpace.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 StateSpace
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


static inline int EquivalentStrings(string a, string b)	{

	if (a.length() != b.length())	{
		return 0;
	}
	unsigned int k = 0;
	int cont = 1;
	while ((k < a.length()) && (cont))	{
		char ca = a[k];
		char cb = b[k];
		if ((ca >=65) && (ca <= 90))	{
			ca += 32;
		}
		if ((cb >=65) && (cb <= 90))	{
			cb += 32;
		}
		if (ca != cb)	{
			cont = 0;
		}
		k++;
	}
	return cont;	
}

string SimpleStateSpace::GetState(int state)	{
	ostringstream s;
	if (state == unknown)	{
		s << "-";
	}
	else	{
		s << Alphabet[state];
	}
	return s.str();
}

int SimpleStateSpace::GetState(string from) {
	if (from.length() != 1)	{
		cerr << "error in SingleLetterStateSpace\n";
		exit(1);
	}
	char c = from[0];
	int p = 0;
	while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
	if (p == NAlphabetSet)	{
		cout << "error: does not recognise character " << c << '\n';
		exit(1);
	}
	if (p >= 2*Nstate)	{
		return unknown;
	}
	else	{
		int k = 0;
		for (int l=0; l<Nstate; l++)		{
			if ((c == Alphabet[l]) || (c == Alphabet[l]+32))	{
				k = l;
			}
		}
		return k;
	}
	return 0;
}

DNAStateSpace::DNAStateSpace()	{
	Nstate = 4;
	Alphabet = new char[Nstate];
	for (int i=0; i<Nstate; i++)	{
		Alphabet[i] = DNAletters[i];
	}
	NAlphabetSet = DNAN;
	AlphabetSet = new char[NAlphabetSet];
	for (int i=0; i<NAlphabetSet; i++)	{
		AlphabetSet[i] = DNAset[i];
	}
}

DNAStateSpace::~DNAStateSpace()	{
	delete[] Alphabet;
	delete[] AlphabetSet;
}


RNAStateSpace::RNAStateSpace()	{
	Nstate = 4;
	Alphabet = new char[Nstate];
	for (int i=0; i<Nstate; i++)	{
		Alphabet[i] = RNAletters[i];
	}
	NAlphabetSet = RNAN;
	AlphabetSet = new char[NAlphabetSet];
	for (int i=0; i<NAlphabetSet; i++)	{
		AlphabetSet[i] = RNAset[i];
	}
}

RNAStateSpace::~RNAStateSpace()	{
	delete[] Alphabet;
	delete[] AlphabetSet;
}

ProteinStateSpace::ProteinStateSpace()	{
	Nstate = 20;
	Alphabet = new char[Nstate];
	for (int i=0; i<Nstate; i++)	{
		Alphabet[i] = AminoAcids[i];
	}
	NAlphabetSet = AAN;
	AlphabetSet = new char[NAlphabetSet];
	for (int i=0; i<NAlphabetSet; i++)	{
		AlphabetSet[i] = AAset[i];
	}
}

ProteinStateSpace::~ProteinStateSpace()	{
	delete[] Alphabet;
	delete[] AlphabetSet;
}


