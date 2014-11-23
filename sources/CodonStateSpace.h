#ifndef CODONSTATESPACE_H
#define CODONSTATESPACE_H

#include "StateSpace.h"
#include <map>

class CodonStateSpace : public StateSpace	{

	public:

	static const int Npos = 3;

	// by default, codons always exclude stops
	// if a method takes or returns a codon stops INCLUDED, then this is made explicit in the method's name

	CodonStateSpace(GeneticCodeType incode);
	~CodonStateSpace();

	// -----
	// generic methods
	// exist for any state space, and should have a consistent meaning throughout

	int GetNstate() {return Nstate;}

	// give a three letter code, returns codon (if stop exits with error message)
	int GetState(string state);

	// give a codon (stops excluded), returns a three letter code
	string GetState(int state);


	// -----
	// codon specific methods

	DNAStateSpace* GetDNAStateSpace() {return nucstatespace;}

	ProteinStateSpace* GetProteinStateSpace() {return protstatespace;}

	// returns a codon based on three letters
	// returns -1 (== unknown) if at least one of the positions is unknown
	// if stop exits with error message...
	int GetCodonFromDNA(int pos1, int pos2, int pos3);

	/*
	string TranslateDNASequenceWithStops(string s);
	string GetStateWithStops(int state);
	int isCodingSequence(string s);
	*/

	// 2 codons excluding stops are compared
	// method returns -1 if identical
	// returns 3 is codons differ at more than one position
	// otherwise, returns the position at which codons differ (i.e. returns 0,1 or 2 if the codons differ at position 1,2 or 3)
	int GetDifferingPosition(int codon1, int codon2);

	// return the integer encoding for the base at requested position
	// stops excluded
	int GetCodonPosition(int pos, int codon)	{
		if ((pos<0) || (pos>=Npos))	{
			cerr << "GetCodonPosition: pos out of bound\n";
			cerr << pos << '\n';
			exit(1);
		}
		if (codon == -1)	{
			return -1;
		}
		if ((codon<0) || (codon>=Nstate))	{
			cerr << "GetCodonPosition: codon out of bound\n";
			cerr << codon << '\n';
			exit(1);
		}
		return CodonPos[pos][codon];
	}

	int IsNonCTNearest(int aminoacid1, int aminoacid2);

	// translation stops excluded
	int Translation(int codon)	{
		return CodonCode[codon];
	}

	// stops excluded
	bool Synonymous(int codon1, int codon2)	{
		return (CodonCode[codon1] == CodonCode[codon2]);
	}

	// returns -1 if stop codon
	// otherwise returns integer in [0,19] standing for an amino-acid (one letter code, alphabetical order)
	int TranslationWithStops(int codon)	{
		return CodonCodeWithStops[codon];
	}

	bool CheckStop(int pos1, int pos2, int pos3);

	/*
	cannot exist: indexing system excludes stop codons anyway...
	bool isStop(int codon)	{
		int n = 0;
		while ((n < Nstop) && (codon != StopCodons[n]))	{
			n++;
		}
		return (n != Nstop);
	}
	*/

	int GetDegeneracy(int codon);

	int GetNstop()	{
		return Nstop;
	}

	const int* GetStopPos1()	{
		return StopPos1;
	}

	const int* GetStopPos2()	{
		return StopPos2;
	}

	const int* GetStopPos3()	{
		return StopPos3;
	}

	private:

	void MakeDegeneracyMap();

	GeneticCodeType code;
	DNAStateSpace* nucstatespace;
	ProteinStateSpace* protstatespace;
	// number of codons, not including stops (61 in general)
	int Nstate;

	// and array of size Ncodon = 64
	// whose entries are between -1 and 19
	// -1 : stop codon
	// 0..19 : amino acid encoded (1 letter code, alphabetical order)
	int* CodonCodeWithStops;
	int* CodonCode;
	int** CodonPos;
	int* StopCodons;
	int Nstop;
	int* StopPos1;
	int* StopPos2;
	int* StopPos3;

	map<int,int> degeneracy;

};

#endif

