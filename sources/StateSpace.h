
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef STATESPACE_H
#define STATESPACE_H

#include "BiologicalSequences.h"

// pure interface
//
class StateSpace	{

	public:

	virtual ~StateSpace() {}

	virtual string GetState(int state) = 0;
	virtual int GetNstate() = 0;

	virtual int GetState(string from) = 0;

};

// simple state space: assumes that states are referred to using a one-letter code
//
class SimpleStateSpace : public StateSpace	{

	public:

	SimpleStateSpace() {}

	SimpleStateSpace(int inNstate, int inNAlphabetSet, char* inAlphabet, char* inAlphabetSet)	{

		Nstate = inNstate;
		NAlphabetSet = inNAlphabetSet;
		Alphabet = new char[Nstate];
		AlphabetSet = new char[NAlphabetSet];
		for (int k=0; k<Nstate; k++)	{
			Alphabet[k] = inAlphabet[k];
		}
		for (int k=0; k<NAlphabetSet; k++)	{
			AlphabetSet[k] = inAlphabetSet[k];
		}
	}

	int GetState(string from);

	int GetNstate() {
		return Nstate;
	}
	
	string GetState(int state);

	char GetCharState(int state) {return AlphabetSet[state];}

	protected:
	int Nstate;
	char* Alphabet;
	int NAlphabetSet;
	char* AlphabetSet;
};

class DNAStateSpace : public SimpleStateSpace	{

	public:

	DNAStateSpace();
	~DNAStateSpace();
};

class RNAStateSpace : public SimpleStateSpace	{

	public:

	RNAStateSpace();
	~RNAStateSpace();
};

class ProteinStateSpace : public SimpleStateSpace	{

	public:

	ProteinStateSpace();
	~ProteinStateSpace();
};

#endif // STATESPACE_H

