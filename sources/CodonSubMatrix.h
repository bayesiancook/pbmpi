
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONSUBMATRIX_H
#define CODONSUBMATRIX_H
#include "CodonStateSpace.h"
#include "SubMatrix.h"

static const double TOOSMALL = 1e-30;
static const double TOOLARGE = 500;
static const double TOOLARGENEGATIVE = -50;

class CodonSubMatrix : public SubMatrix	{

	public:

	CodonSubMatrix(CodonStateSpace* instatespace, double* innucrr, double* innucstat, bool innormalise) :
		SubMatrix(instatespace->GetNstate(),innormalise), 
		statespace(instatespace),
		nucrr(innucrr),
		nucstat(innucstat),
		normalise(innormalise)	{}

	CodonStateSpace* GetCodonStateSpace() {return statespace;}
	// this is just to avoid repeating "GetCodonStateSpace()->" all the time...
	int GetNstate() {return statespace->GetNstate();}
	bool Synonymous(int codon1, int codon2) {return statespace->Synonymous(codon1,codon2);}
	int GetCodonPosition(int pos, int codon) {return statespace->GetCodonPosition(pos,codon);}
	int GetDifferingPosition(int codon1, int codon2) {return statespace->GetDifferingPosition(codon1,codon2);}
	int GetNucRRIndex(int i, int j)	{return (i<j) ? (2 * Nnuc - i - 1) * i / 2 + j - i - 1 : (2 * Nnuc - j - 1) * j / 2 + i - j - 1 ;}			

	protected:

	void ComputeArray(int state);
	void ComputeStationary();

	CodonStateSpace* statespace;
	double* nucrr;
	double* nucstat;
	bool normalise;		
};


class AAMutSelProfileSubMatrix : public CodonSubMatrix	{

	public:

	AAMutSelProfileSubMatrix(CodonStateSpace* instatespace, double* innucrr, double* innucstat, double* inaaprofile, bool innormalise) :
		CodonSubMatrix(instatespace,innucrr,innucstat,innormalise),
		aaprofile(inaaprofile) {}

	double* GetAAProfile() {return aaprofile;}

	protected:

	void ComputeArray(int state);
	void ComputeStationary();
	double GetRate();
	double* aaprofile;


};

class AACodonMutSelProfileSubMatrix : public CodonSubMatrix	{

	public:

	AACodonMutSelProfileSubMatrix(CodonStateSpace* instatespace, double* innucrr, double* innucstat, double* incodonprofile, double* inaaprofile, double* inomega, bool innormalise) :
		CodonSubMatrix(instatespace,innucrr,innucstat,innormalise),
		codonprofile(incodonprofile),
		aaprofile(inaaprofile),
		omega(inomega) {}

	double* GetAAProfile() {return aaprofile;}
	double* GetCodonProfile() {return codonprofile;}

	protected:

	void ComputeArray(int state);
	void ComputeStationary();
	double GetRate();
	double* codonprofile;
	double* aaprofile;
	double* omega;
};

class CodonMutSelProfileSubMatrix : public CodonSubMatrix	{

	public:

	CodonMutSelProfileSubMatrix(CodonStateSpace* instatespace, double* innucrr, double* innucstat, double* incodonprofile, bool innormalise) :
		CodonSubMatrix(instatespace,innucrr,innucstat,innormalise),
		codonprofile(incodonprofile) {}

	double* GetCodonProfile() {return codonprofile;}

	protected:

	void ComputeArray(int state);
	void ComputeStationary();
	double GetRate();
	double* codonprofile;

};

#endif
