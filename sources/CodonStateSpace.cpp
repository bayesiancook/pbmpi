
#include <iostream>
#include<sstream>
#include <cstdlib>
using namespace std;

#include "CodonStateSpace.h"

CodonStateSpace::CodonStateSpace(GeneticCodeType type)	{

	nucstatespace = new DNAStateSpace;
	protstatespace = new ProteinStateSpace;

	code = type;
	if (code == Universal)	{
		Nstate = Ncodon - UniNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = UniCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}

		Nstop = UniNStopCodons;
		StopCodons = new int[Nstop];
		StopPos1 = new int[Nstop];
		StopPos2 = new int[Nstop];
		StopPos3 = new int[Nstop];
		for (int i=0; i<Nstop; i++)	{
			StopCodons[i] = UniStopCodons[i];
			StopPos1[i] = UniStopPos1[i];
			StopPos2[i] = UniStopPos2[i];
			StopPos3[i] = UniStopPos3[i];
		}

	}
	else if (code == MtInv)	{
		Nstate = Ncodon - MtInvNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = MtInvCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}

		Nstop = MtInvNStopCodons;
		StopCodons = new int[Nstop];
		StopPos1 = new int[Nstop];
		StopPos2 = new int[Nstop];
		StopPos3 = new int[Nstop];
		for (int i=0; i<Nstop; i++)	{
			StopCodons[i] = MtInvStopCodons[i];
			StopPos1[i] = MtInvStopPos1[i];
			StopPos2[i] = MtInvStopPos2[i];
			StopPos3[i] = MtInvStopPos3[i];
		}
	}
	else if (code == MtMam)	{
		Nstate = Ncodon - MtMamNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = MtMamCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}

		Nstop = MtMamNStopCodons;
		StopCodons = new int[Nstop];
		StopPos1 = new int[Nstop];
		StopPos2 = new int[Nstop];
		StopPos3 = new int[Nstop];
		for (int i=0; i<Nstop; i++)	{
			StopCodons[i] = MtMamStopCodons[i];
			StopPos1[i] = MtMamStopPos1[i];
			StopPos2[i] = MtMamStopPos2[i];
			StopPos3[i] = MtMamStopPos3[i];
		}
	}
	else 	{
		cerr << "genetic code not recognised\n";
		cerr << type << '\n';
		exit(1);
	}

}

CodonStateSpace::~CodonStateSpace()	{

	delete[] CodonCode;
	delete[] CodonCodeWithStops;
	for (int pos=0; pos<Npos; pos++)	{
		delete[] CodonPos[pos];
	}
	delete[] CodonPos;

	delete nucstatespace;
	delete protstatespace;
}

string CodonStateSpace::GetState(int codon)	{
	ostringstream s;
	if (codon == -1)	{
		s << "---";
	}
	else	{
		s << DNAletters[GetCodonPosition(0,codon)] << DNAletters[GetCodonPosition(1,codon)] << DNAletters[GetCodonPosition(2,codon)];
	}
	if (s.str().length() != 3)	{
		cerr << "error in translation\n";
		exit(1);
	}
	return s.str();
}

int CodonStateSpace::GetState(string word)	{
	return GetCodonFromDNA(GetDNAStateSpace()->GetState(word.substr(0,1)),GetDNAStateSpace()->GetState(word.substr(1,1)),GetDNAStateSpace()->GetState(word.substr(2,1)));
}

bool CodonStateSpace::CheckStop(int pos1, int pos2, int pos3)	{
	if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown))	{
		return false;
	}
	int l = 0;
	while ((l < Nstop) && ((pos1 != StopPos1[l]) || (pos2 != StopPos2[l]) || (pos3 != StopPos3[l])))	{
		l++;
	}
	return (l < Nstop);
}

int CodonStateSpace::GetCodonFromDNA(int pos1, int pos2, int pos3)	{
	if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown))	{
		return unknown;
	}
	int l = 0;
	while ((l<GetNstate()) && ((pos1 != GetCodonPosition(0,l)) || (pos2 != GetCodonPosition(1,l)) || (pos3 != GetCodonPosition(2,l))))	{
		l++;
	}
	if (l == GetNstate())	{
		// cerr << "warning in CodonStateSpace::GetCodonFromDNA : out of bound : " << GetDNAStateSpace()->GetState(pos1) << GetDNAStateSpace()->GetState(pos2) << GetDNAStateSpace()->GetState(pos3) << '\n';
		return -1;
		cerr << "warning in CodonStateSpace::GetCodonFromDNA : out of bound : " << GetDNAStateSpace()->GetState(pos1) << GetDNAStateSpace()->GetState(pos2) << GetDNAStateSpace()->GetState(pos3) << '\n';
		if (code == Universal)	{
			cerr << "universal\n";
		}
		else if (code == MtMam)	{
			cerr << "mt mam\n";
		}
		else if (code == MtInv)	{
			cerr << "mt inv\n";
		}
		cerr << "code not recognized\n";
		exit(1);
	}
	return l;

}

int CodonStateSpace::GetDifferingPosition(int i, int j)	{

	// identical
	if ((GetCodonPosition(0,i) == GetCodonPosition(0,j)) && (GetCodonPosition(1,i) == GetCodonPosition(1,j)) && (GetCodonPosition(2,i) == GetCodonPosition(2,j)))	{
		return -1;
	}
	if (GetCodonPosition(0,i) != GetCodonPosition(0,j))	{
		if ((GetCodonPosition(1,i) == GetCodonPosition(1,j)) && (GetCodonPosition(2,i) == GetCodonPosition(2,j)))	{
			return 0;
		}
		else	{
			return 3;
		}
	}
	if (GetCodonPosition(1,i) != GetCodonPosition(1,j))	{
		if ((GetCodonPosition(0,i) == GetCodonPosition(0,j)) && (GetCodonPosition(2,i) == GetCodonPosition(2,j)))	{
			return 1;
		}
		else	{
			return 3;
		}
	}
	if (GetCodonPosition(2,i) != GetCodonPosition(2,j))	{
		if ((GetCodonPosition(1,i) == GetCodonPosition(1,j)) && (GetCodonPosition(0,i) == GetCodonPosition(0,j)))	{
			return 2;
		}
		else	{
			return 3;
		}
	}
	return 3;
}

int CodonStateSpace::GetDegeneracy(int codon)	{

	if (!degeneracy.size())	{
		MakeDegeneracyMap();
	}
	if (codon == -1)	{
		return -1;
	}
	return degeneracy[codon];
};

void CodonStateSpace::MakeDegeneracyMap()	{

	for (int codon =0; codon < Nstate; codon++)	{
		cerr << codon << '\t' << GetState(codon) << '\t';
		int pos1 = GetCodonPosition(0,codon);
		int pos2 = GetCodonPosition(1,codon);
		int aa = Translation(codon);
		int d = 0;
		for (int n=0; n<Nnuc; n++)	{
			int cod = GetCodonFromDNA(pos1,pos2,n);
			if ((cod != -1) && (Translation(cod) == aa))	{
				d++;
			}
			cerr << GetState(cod);
			if (cod == -1)	{
				cerr << " $ ";
			}
			else	{
				cerr << ' ' << GetProteinStateSpace()->GetState(Translation(cod)) << ' ';
			}
		}
		cerr << '\t';
		cerr << d << '\n';
		degeneracy[codon] = d;
	}
	cerr << "ok\n";
};

int CodonStateSpace::IsNonCTNearest(int a, int b)	{

	int noct = 1;
	int nn = 0;
	for (int c1 = 0; c1 <GetNstate(); c1++)	{
		if (Translation(c1) == a)	{
			for (int c2 = 0; c2 <GetNstate(); c2++)	{
				if (Translation(c2) == b)	{
					int pos = GetDifferingPosition(c1,c2);
					if (pos < 3)	{
						nn = true;
						int n1 = GetCodonPosition(pos,c1);
						int n2 = GetCodonPosition(pos,c2);
						if (((n1 == 1) && (n2 == 3)) || ((n1 == 3) && (n2 == 1)))	{
							noct = false;
						}
					}
				}
			}
		}
	}
	if (! nn)	{
		return -1;
	}
	return noct;
}

/*
int CodonStateSpace::GetStateWithStops(string word)	{
	return GetCodonFromDNAWithStops(GetDNAStateSpace()->GetState(word.substr(0,1)),GetDNAStateSpace()->GetState(word.substr(1,1)),GetDNAStateSpace()->GetState(word.substr(2,1)));
}

int CodonStateSpace::GetCodonFromDNAWithStops(int pos1, int pos2, int pos3)	{
	if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown))	{
		return unknown;
	}
	int l = 0;
	while ((l<Ncodon) && ((pos1 != codonpos[0][l]) || (pos2 != codonpos[1][l]) || (pos3 != codonpos[2][l])))	{
		l++;
	}
	if (l == Ncodon)	{
		cerr << "error in CodonStateSpace::GetCodonFromDNAWithStops : out of bound : " << pos1 << '\t' << pos2 << '\t' << pos3 << '\n';
		exit(1);
	}
	return l;
}
string GetStateWithStops(int codon)	{

	if (codon == unknown)	{
		return "?";
	}
	int aa = TranslationWithStops(codon);
	if (aa == -1)	{
		return "O";
	}
	return GetProteinStateSpace()->GetState(aa);

}

string CodonStateSpace:TranslateDNASequenceWithStops(string dnaseq)	{

	if (dnaseq.length() % 3)	{
		cerr << "error in CodonStateSpace::Translate: dna sequence not multiple of three\n";
		cerr << "length is : " << dnaseq.length() << '\n';
		exit(1);
	}
	int N = dnaseq.length() / 3;

	ostringstream s;
	for (int i=0; i<N; i++)	{
		s << GetStateWithStops(TranslationWithStops(GetStateWithStops(dnaseq.substr(3*i,3))));
	}
	return s.str();
}
*/



