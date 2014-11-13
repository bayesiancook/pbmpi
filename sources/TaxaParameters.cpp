
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "phylo.h"


// ---------------------------------------------------------------------------------
//		 TaxaParameters()
// ---------------------------------------------------------------------------------

TaxaParameters::TaxaParameters()	{

	Ntaxa = 0;
	SpeciesNames = 0;
	mParam = 0;
}


TaxaParameters::TaxaParameters(string filename)	{

	Ntaxa = 0;
	SpeciesNames = 0;
	mParam = 0;
	
	ifstream is(filename.c_str());
	string temp;
	is >> temp;
	if (temp == "TaxaList")	{
		ReadFromStream(is);
	}
	else	{
		cerr << "error in TaxaParameters::constructor from file\n";
		exit(1);
	}
}


// ---------------------------------------------------------------------------------
//		 TaxaParameters(int, string*)
// ---------------------------------------------------------------------------------

TaxaParameters::TaxaParameters(int N, string* inNames)	{

	Ntaxa = N;
	mParam = 0;
	SpeciesNames = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		if (inNames)	{
			SpeciesNames[i] = inNames[i];
		}
		else	{
			SpeciesNames[i] = "";
		}
	}
}


// ---------------------------------------------------------------------------------
//		 TaxaParameters(MCParameters*)
// ---------------------------------------------------------------------------------

TaxaParameters::TaxaParameters(MCParameters* inParam)	{

	/*
	mParam = inParam;	// no copy
	// Ntaxa = mParam->Ntaxa;

	SpeciesNames = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		SpeciesNames[i] = mParam->SpeciesNames[i];
	}
	*/
	
}

// ---------------------------------------------------------------------------------
//		 TaxaParameters(const TaxaParameters& from)
// ---------------------------------------------------------------------------------


TaxaParameters::TaxaParameters(const TaxaParameters& from)	{

	mParam = from.mParam;
	Ntaxa = from.Ntaxa;
	
	SpeciesNames = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		SpeciesNames[i] = from.SpeciesNames[i];
	}
	
}

int TaxaParameters::RegisterWith(TaxaParameters* inParam)	{
	int ok = 1;
	if (! inParam)	{
		cerr << "error : non existing taxon set\n";
		exit(1);
	}
	if (Ntaxa != inParam->Ntaxa)	{
		ok = 0;
	}
	int k = 0;
	while (ok && (k<Ntaxa))	{
		int j = 0;
		while ((j<Ntaxa) && (SpeciesNames[k] != inParam->SpeciesNames[j])) j++;
		if (j == Ntaxa)	{
			ok = 0;
		}
		k++;
	}
	int j = 0;
	while (ok && (j<Ntaxa))	{
		int k = 0;
		while ((k<Ntaxa) && (SpeciesNames[k] != inParam->SpeciesNames[j])) k++;
		if (k == Ntaxa)	{
			ok = 0;
		}
		j++;
	}
	return ok;
}
		
	

// ---------------------------------------------------------------------------------
//		 TaxaParameters(Tree* inTree)
// ---------------------------------------------------------------------------------

TaxaParameters::TaxaParameters(PBTree* inTree)	{

	Ntaxa = inTree->GetSize();
	mParam = 0;
	SpeciesNames = new string[Ntaxa];
	inTree->GetSpeciesNames(SpeciesNames);
	inTree->mParam = this;
	inTree->RegisterWithParam(this);
}


// ---------------------------------------------------------------------------------
//		 ~TaxaParameters()
// ---------------------------------------------------------------------------------

TaxaParameters::~TaxaParameters()	{
	
	delete[] SpeciesNames;
}


// ---------------------------------------------------------------------------------
//		 GetSpeciesName
// ---------------------------------------------------------------------------------

string	TaxaParameters::GetSpeciesName(int index)	{
	return SpeciesNames[index];
}

// ---------------------------------------------------------------------------------
//		 Streams
// ---------------------------------------------------------------------------------
void TaxaParameters::ReadNexus(istream& is)	{

	/*
	string tempname[NtaxaMax];
	string temp;
	int N = 0;
	while ((! is.eof()) && (temp != "Translate"))	{
		is >> temp;
	}
	if (is.eof())	{
		cerr << "error when reading tree list from nexus filz: did not find 'Translate' keyword\n";
		exit(1);
	}
	
	do {
		is >> temp;
		if (temp != ";")	{
			is >> temp;
			tempname[N++] = temp;
		}
	}
	while (temp != ";");
	Ntaxa = N;
	SpeciesNames = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		SpeciesNames[i] = Filter(tempname[i], ',');
	}
	*/
}	
		

void TaxaParameters::ReadFromStream(istream& is)	{

	string temp;
	is >> temp;
	if (temp != "Ntaxa")	{
		cerr << "error in TaxaParameters::ReadFromStream\n";
		exit(1);
	}
	is >> Ntaxa;
	delete[] SpeciesNames;
	
	SpeciesNames = new string[Ntaxa];
	int SpeciesDone = 0;

	is >> temp;
	while (temp != "End")	{
		if (temp == "Names")	{
			for (int i=0; i<Ntaxa; i++)	{
				is >> SpeciesNames[i];
			}
			SpeciesDone = 1;
			is >> temp;
		}
		else	{
			cerr << "error in TaxaParameters::ReadFromStream : does not recognise " << temp << '\n';
			exit(1);
		}
	}

	if (! SpeciesDone)	{
		for (int i=0; i<Ntaxa; i++)	{
			ostringstream s;
			s << i;
			SpeciesNames[i] = s.str();
		}
	}
}

void TaxaParameters::WriteToStream(ostream& os)	{

	os << "TaxaList\n\n";
	os << "Ntaxa " << Ntaxa << '\n';
	os << "Names\n";
	for (int i=0; i<Ntaxa; i++)	{
		os << SpeciesNames[i] << '\n';
	}
	os << "End\n\n";	
}

