
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
const int NTreeMax = 10000;

// ---------------------------------------------------------------------------------
//		 TreeList(string filename)
// ---------------------------------------------------------------------------------

TreeList::TreeList(string filename)	{

	mParam = 0;
	mSize = 0;
	mTreeArray = 0;
	ifstream is(filename.c_str());
	ReadFromStream(is);
	SetParameters();
}


// ---------------------------------------------------------------------------------
//		 TreeList(TaxaParameters*, int)
// ---------------------------------------------------------------------------------

TreeList::TreeList(TaxaParameters* inParam, int inSize)	{

	mParam = inParam;
	mSize = inSize;

	mTreeArray = new PBTree*[mSize];

	for (int i=0; i<mSize; i++)	{
		mTreeArray[i] = 0;
	}
}


// ---------------------------------------------------------------------------------
//		 ~TreeList()
// ---------------------------------------------------------------------------------

TreeList::~TreeList()	{

	if (mTreeArray)	{	
		for (int i=0; i<mSize; i++)	{
			delete mTreeArray[i];
		}
	}

	delete[] mTreeArray;
}


// ---------------------------------------------------------------------------------
//		 SetParameters()
// ---------------------------------------------------------------------------------

void	TreeList::SetParameters()	{

	if (!mParam)	{
		cerr << "error in TreeList::SetParameters : mParam == 0\n";
		exit(1);
	}
	for (int i=0; i<mSize; i++)	{
		if (! mTreeArray[i])	{
			cerr << "error in TreeList::SetParameters : found a null tree\n";
			exit(1);
		}
		mTreeArray[i]->SetParameters(mParam);
		if (mTreeArray[i]->Name == "T")	{
			ostringstream s;
			s << i;
			mTreeArray[i]->Name = s.str();
		}
	}
}	
	
// ---------------------------------------------------------------------------------
//		 GetTree(int index)
// ---------------------------------------------------------------------------------

PBTree*	TreeList::GetTree(int index)	{

	if ((index < 0) || (index > mSize))	{
		cerr << "error in tree array : try to access a tree out of range\n";
		exit(1);
	}

	return mTreeArray[index];
}

// ---------------------------------------------------------------------------------
//		 WriteToStream(ostream& os, int header)
// ---------------------------------------------------------------------------------

void TreeList::WriteToStream(ostream& os, int header, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{

	if (header)	{
		mParam->WriteToStream(os);
	}
	for (int i=0; i<mSize; i++)	{
		mTreeArray[i]->WriteToStream(os,header,withLengths,withProbs,withSpeciesNames,withInternalLabels);
	}
}

// ---------------------------------------------------------------------------------
//		 ReadFromStream(istream& is)
// ---------------------------------------------------------------------------------

void TreeList::ReadFromStream(istream& is)	{

	string temp;
	is >> temp;

	if (temp == "#NEXUS")	{
		mParam = new TaxaParameters();
		mParam->ReadNexus(is);

		mTreeArray = new PBTree*[NTreeMax];
		int N = 0;
		do	{

			while ((! is.eof()) && (temp != "tree"))	{
				is >> temp;
			}
			if (is.eof())	{
				cerr << "error when reading tree list from nexus\n";
				exit(1);
			}
		
			is >> temp >> temp >> temp >> temp >> temp;

			mTreeArray[N] = new PBTree(mParam);
			mTreeArray[N]->ReadFromStream(is,1);
			mTreeArray[N]->SetParameters(mParam);
			mTreeArray[N]->FinishCreation();
			N++;
			is >> temp;
		}
		while (temp != "end;");

		cerr << "found " << N << "trees\n";
		mSize = N;
		PBTree** array = new PBTree*[mSize];
		for (int i=0; i<mSize; i++)	{
			array[i] = mTreeArray[i];
		}
		delete[] mTreeArray;
		mTreeArray = array;
		cerr << "tree list : ok\n";
	}

	else	{
		while (temp != "End")	{
			if (temp == "TaxaList")	{
				mParam = new TaxaParameters();
				mParam->ReadFromStream(is);
			}
			else if (temp == "TreeList")	{
				is >> temp >> mSize;
				cerr << temp << '\t' << mSize << '\n';
				delete[] mTreeArray;
				mTreeArray = new PBTree*[mSize];
				for (int i=0; i<mSize; i++)	{
					mTreeArray[i] = new PBTree(mParam);
					mTreeArray[i]->ReadFromStream(is);
					mTreeArray[i]->SetParameters(mParam);
				}
			}
			else	{
				cerr << "error in TreeList::ReadFromStream\n";
				exit(1);
			}
			is >> temp;
		}
	}
}


// ---------------------------------------------------------------------------------
//		 RootAt()
// ---------------------------------------------------------------------------------

void TreeList::RootAt(Bipartition outgroup)	{
	for (int i=0; i<mSize; i++)	{
		mTreeArray[i]->RootAt(outgroup);
	}
}

// ---------------------------------------------------------------------------------
//		 ToPS(string filename)
// ---------------------------------------------------------------------------------

void	TreeList::ToPS(string target, int every, double sizeX, double sizeY, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{

	/*
	// tex output ?
	string texfile = target + ".tex";
	string appel = "cp " + header + " " + texfile;
	system(appel.c_str());
	ofstream Tex_os(texfile.c_str(), IOS_APPEND);

	Tex_os << "\\begin{document}\n";
	for (int i=0; i<mSize; i+=every)	{
		Tex_os << "\\noindent\n";
		Tex_os << "tree " << i << '\n';
		Tex_os << "\\\\\n";
		mTreeArray[i]->ToLatex(Tex_os, sizeX, sizeY, withLengths, withProbs, withSpeciesNames, withInternalLabels);
		Tex_os << "\\newpage\n\n";
	}
	Tex_os << "\\end{document}\n";
	Tex_os.close();

	string latex = "latex " + texfile + " > tmp";
	string dvips = "dvips -o " + target + ".ps " + target + " 2> tmp";
	string rm = "rm -f tmp";
	string rm2 = "rm -f " + target + ".aux";
	string rm3 = "rm -f " + target + ".log";
	string rm4 = "rm -f " + target + ".dvi";
	system(latex.c_str());
	system(dvips.c_str());
	system(rm.c_str());
	system(rm2.c_str());
	system(rm3.c_str());
	system(rm4.c_str());
	*/

}


/*

// ---------------------------------------------------------------------------------
//		 ToMrBayes(ostream& os, int every)
// ---------------------------------------------------------------------------------

void	TreeList::ToMrBayes(ostream& MB_os, int every)	{

	int Ntaxa = mParam->Nspecies;
	//preparing mr bayes output

	MB_os << "#NEXUS\n";
	MB_os << "[ from phylobayes ]\n";
	MB_os << "begin trees;\n";
	MB_os << "translate\n";
	for (int i=0; i<Ntaxa; i++)	{
		MB_os << i+1 << ' ' << mParam->GetSpeciesName(i);
		if (i == (Ntaxa-1))	{
			MB_os << ";\n";
		}
		else	{
			MB_os << ",\n";
		}
	}

	int j = 0;
	int cycle = 0;
	for (int i=0; i<GetSize(); i++)	{
		if (! cycle)	{
			MB_os << "\ttree rep." << j+1 << " = ";
			GetTree(i)->Trichotomise();
			GetTree(i)->ToMrBayes(MB_os);
			MB_os << ";\n";
			j++;
		}
		cycle ++;
		if (cycle == every)	{
			cycle = 0;
		}
	}

	MB_os << "end;\n" ;
}

// ---------------------------------------------------------------------------------
//		 ReadMrBayes(string filename)
// ---------------------------------------------------------------------------------

void	TreeList::ReadMrBayes(string filename)	{

	ifstream is(filename.c_str());
	// first , compute size
	is.ReadLine();
	is.ReadLine();
	is.ReadLine();
	is.ReadLine();
	int Ntaxa = 0;
	string temp;
	is >> temp;
	while (temp != "tree")	{
		Ntaxa ++;
		is.ReadLine();
		is >> temp;
	}
	int size = 0;
	while (! is.eof())	{
		is.ReadLine();
		size++;
	}
	is.close();

	size -= 2;

	ifstream is2(filename.c_str());
	ReadMrBayes(is2, Ntaxa, size);

}

// ---------------------------------------------------------------------------------
//		 ReadMrBayes(istream& is)
// ---------------------------------------------------------------------------------

void	TreeList::ReadMrBayes(istream& is, int Ntaxa, int size)	{

	if (mTreeArray)	{
		for (int i=0; i<mSize; i++)	{
			delete mTreeArray[i];
		}
	}
	delete[] mTreeArray;
	delete mParam;

	is.ReadLine();
	is.ReadLine();
	is.ReadLine();
	is.ReadLine();

	string temp;
	string* names = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		is >> temp >> temp;
		if (! temp.length())	{
			cerr << "error in read tree list  from  mrbayes\n";
			exit(1);
		}
		names[i] = string(temp,0,temp.length() -1);
	}
	mParam = new TaxaParameters(Ntaxa, names);
	mSize = size;
	mTreeArray = new Tree*[mSize];
	for (int i=0; i<size; i++)	{
		is >> temp >> temp >> temp;
		mTreeArray[i] = new Tree(is);
		mTreeArray[i]->SetLabelOffset(0);
		mTreeArray[i]->mParam = mParam;
		mTreeArray[i]->Dichotomise();
		mTreeArray[i]->mRoot->SetSuper(mTreeArray[i]);
	}

	delete[] names;
}

*/
 
