
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
//		 Bipartition(TaxaParameters*)
// ---------------------------------------------------------------------------------

Bipartition::Bipartition(TaxaParameters* inParam)	{

	mParam = inParam;
	Ntaxa = mParam->Ntaxa;
	
	mArray = new int[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		mArray[i] = 0;
	}

}

// ---------------------------------------------------------------------------------
//		 Bipartition( const Bipartition& )
// ---------------------------------------------------------------------------------

Bipartition::Bipartition(const Bipartition& from)	{

	mParam = from.mParam;
	Ntaxa = mParam->Ntaxa;
	mArray = new int[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		mArray[i] = from.mArray[i];
	}
}


// ---------------------------------------------------------------------------------
//		 ~Bipartition()
// ---------------------------------------------------------------------------------

Bipartition::~Bipartition()	{

	delete[] mArray;
}

// ---------------------------------------------------------------------------------
//		 operator=
// ---------------------------------------------------------------------------------

Bipartition&	Bipartition::operator=(const Bipartition& from)	{

	if (this != & from)	{
		// assume they have the same TaxaParameters
		mParam = from.mParam;
		for (int i=0; i<Ntaxa; i++)	{
			mArray[i] = from.mArray[i];
		}
	}
	return *this;
}

// ---------------------------------------------------------------------------------
//		 operator=
// ---------------------------------------------------------------------------------

Bipartition&	Bipartition::operator=(const string& from)	{

	for (int i=0; i<Ntaxa; i++)	{
		if (from[i] == '.')	{		// inside
			mArray[i] = 0;
		}
		else if (from[i] == '*')	{	// outside
			mArray[i] = 1;
		}
		else if (from[i] == ' ')	{	// absent
			mArray[i] = -1;
		}
		else	{
			cerr << "error in Bipartition::operator=(string)\n";
			exit(1);
		}
	}
	return *this;
}

// ---------------------------------------------------------------------------------
//		 Bipartition( const Bipartition& )
// ---------------------------------------------------------------------------------

void Bipartition::AllAbsent()	{
	for (int i=0; i<Ntaxa; i++)	{
		mArray[i] = -1;
	}
}

void Bipartition::Suppress(const Bipartition& leafset)	{
	for (int i=0; i<Ntaxa; i++)	{
		if (leafset.mArray[i] == -1)	{
			if (mArray[i] != 0)	{
				cerr << "?? in Bipartition::Suppress: suppressing a non zero taxon\n";
				exit(1);
			}
			mArray[i] = -1;
		}
	}
}

void Bipartition::PermutTaxa(int* permut)	{

	int bk[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		bk[i] = mArray[i];
	}
	for (int i=0; i<Ntaxa; i++)	{
		mArray[permut[i]] = bk[i];
	}	
	Modulo();
}

	

// ---------------------------------------------------------------------------------
//		 CompareWith()
// ---------------------------------------------------------------------------------

int Bipartition::CompareWith(const Bipartition& with) {

	Bipartition bp(with.mParam);
	for (int i=0; i<with.Ntaxa; i++)	{
		bp.mArray[i] = -1;
	}
	for (int i=0; i<Ntaxa; i++)	{
		string name = mParam->SpeciesNames[i];
		int k = 0;
		while ((k<with.Ntaxa) && (name != with.mParam->SpeciesNames[k]))	{
			k++;
		}
		if (k == with.Ntaxa)	{
			// cerr << "error in Bipartition::CompareWith: overflow\n";
			// exit(1);
		}
		else	{
			bp.mArray[k] = mArray[i];
		}
	}
	return bp.IsCompatibleWith(with);
}

int Bipartition::SupportCheck(const Bipartition& with)	{

	Bipartition bp(mParam);
	for (int i=0; i<Ntaxa; i++)	{
		string name = mParam->SpeciesNames[i];
		int k = 0;
		while ((k<with.Ntaxa) && (name != with.mParam->SpeciesNames[k]))	{
			k++;
		}
		if (k == with.Ntaxa)	{
			cerr << "error in Bipartition::SupportCheck: overflow\n";
			exit(1);
		}
		bp.mArray[i] = with.mArray[k];
	}
	return bp == *this;
}

int Bipartition::CompareWith(BipartitionList* bplist)	{

	int compatible = 1;
	int i = 0;
	while ((i<bplist->GetSize()) && (compatible))	{
		compatible &= CompareWith(bplist->GetBipartition(i));
		i++;
	}
	return compatible;
}

int Bipartition::SupportCheck(BipartitionList* bplist)	{

	int found = 0;
	int i = 0;
	while ((i<bplist->GetSize()) && (! found))	{
		if (bplist->mBipartitionArray[i]->IsInformative())	{
			if (SupportCheck(bplist->GetBipartition(i)))	{
				found = 1;
			}
		}
		i++;
	}
	return found;
}

// ---------------------------------------------------------------------------------
//		 GetTaxonStatus(int index)
// ---------------------------------------------------------------------------------

int	Bipartition::GetTaxonStatus(int index) {

	return (index == -1) ? -3 : mArray[index];
}

// ---------------------------------------------------------------------------------
//		 SetTaxon(int index)
// ---------------------------------------------------------------------------------

void	Bipartition::SetTaxon(int index)	{
	if (mArray[index] == -1)	{
		cerr << "error in Bipartition::SetTaxon\n";
		exit(1);
	}
	else	{
		mArray[index] = 1;
	}
}

// ---------------------------------------------------------------------------------
//		 operator==
// ---------------------------------------------------------------------------------

Boolean	Bipartition::operator==( const Bipartition& inPartition)	{

	// assume they have the same TaxaParameters
	// assume they are oriented the same way
	Bipartition temp = inPartition;
	Boolean test = true;
	int i=0;
	while (test && (i<Ntaxa))	{
		test &= ((mArray[i] == -1) || (inPartition.mArray[i] == -1) || (mArray[i] == inPartition.mArray[i]));
		i++;
	}
	return test;
}


// ---------------------------------------------------------------------------------
//		 operator!=
// ---------------------------------------------------------------------------------

Boolean
Bipartition::operator!=( const Bipartition& inPartition)	{

	return ! operator==(inPartition);

}


double Bipartition::GetPriorProb()	{

	int n1 = 0;
	for (int k=0; k<Ntaxa; k++)	{
		if (mArray[k])	n1++;
	}
	int n2 = Ntaxa - n1;
	
	double logtotal = 0;
	for (int i=1; i<=2*n1-3; i+=2)	{
		logtotal += log(i);
	}
	for (int i=1; i<=2*n2-3; i+=2)	{
		logtotal += log(i);
	}
	for (int i=1; i<=2*Ntaxa-5; i+=2)	{
		logtotal -= log(i);
	}
	return exp(logtotal);
}

// ---------------------------------------------------------------------------------
//		 IsCompatibleWith
// ---------------------------------------------------------------------------------

Boolean	Bipartition::IsCompatibleWith( const Bipartition& inPartition)	{
	
	// assume they have the same TaxaParameters
	
	Boolean test = false;
	
	Bipartition temp = *this;
	temp &= inPartition;
	test |= (temp == inPartition) || (temp == *this);
	
	if (! test)	{
	
		Bipartition NotThis = ! *this;
		Bipartition temp2 = NotThis;
		temp2 &= inPartition;
		test |= (temp2 == inPartition) || (temp2 == NotThis);
	}
	
	return test;
}			


// ---------------------------------------------------------------------------------
//		 IsCompatibleWith
// ---------------------------------------------------------------------------------

Boolean	Bipartition::IsCompatibleWith( const Bipartition& inPartition, Boolean& Orientation)	{

	// assume they have the same TaxaParameters

	Boolean test = false;
	Orientation = true;

	Bipartition temp = *this;
	temp &= inPartition;
	test |= (temp == inPartition) || (temp == *this);

	if (! test)	{
	
		Orientation = false;
		Bipartition NotThis = ! *this;
		Bipartition temp2 = NotThis;
		temp2 &= inPartition;
		test |= (temp2 == inPartition) || (temp2 == NotThis);
	}
	
	return test;
}			

// ---------------------------------------------------------------------------------
//		 operator|=
// ---------------------------------------------------------------------------------

Bipartition& Bipartition::operator|=( const Bipartition& inPartition)	{
	
	// assumes they are oriented likewise
	for (int i=0; i<Ntaxa; i++)	{
		if (mArray[i] != -1)	{
			mArray[i]|=inPartition.mArray[i];
		}
	}
	
	return *this;
}

// ---------------------------------------------------------------------------------
//		 operator&=
// ---------------------------------------------------------------------------------

Bipartition& Bipartition::operator&=( const Bipartition& inPartition)	{
	
	// assumes they are oriented likewise
	for (int i=0; i<Ntaxa; i++)	{
		if (mArray[i] != -1)	{
			mArray[i]&=inPartition.mArray[i];
		}
	}
	return *this;
}

// ---------------------------------------------------------------------------------
//		 operator!
// ---------------------------------------------------------------------------------

Bipartition Bipartition::operator!()	{

	Bipartition temp  = *this;
	for (int i=0; i<Ntaxa; i++)	{
		if (mArray[i] != -1)	{
			temp.mArray[i] = 1 - mArray[i];
		}
	}
	return temp;
}


// ---------------------------------------------------------------------------------
//		 IsInformative
// ---------------------------------------------------------------------------------

Boolean	Bipartition::IsInformative()	{

	Boolean temp0 = false;
	Boolean ret0 = false;
	Boolean temp1 = false;
	Boolean ret1 = false;
	int i=0;
	while (i<Ntaxa && ! (ret0 && ret1) )	{
		int ind = mArray[i];
		if (ind != -1)	{
			if (temp0)	{
				ret0 |= ind == 0;
			}
			temp0 |= ind == 0;

			if (temp1)	{
				ret1 |= ind == 1;
			}
			temp1 |= ind == 1;
		}
		i++;
	}
	return (ret0 && ret1);
}

// ---------------------------------------------------------------------------------
//		 Modulo
// ---------------------------------------------------------------------------------

void
Bipartition::Modulo()	{

	int i=0;
	while ((i<Ntaxa) && (mArray[i] == -1))	{
		i++;
	}
	if (mArray[i] == 1)	{
		while (i<Ntaxa)	{
			if (mArray[i] != -1)	{
				mArray[i] = 1 - mArray[i];
			}
			i++;
		}
	}
}


// ---------------------------------------------------------------------------------
//		 WriteToStream()
// ---------------------------------------------------------------------------------

void	Bipartition::WriteToStream(ostream& os, int verbose)	{

	if ((verbose == 2) && mParam)	{
		os << '\n';
		int on = 0;
		int off = 0;
		for (int i=0; i<Ntaxa; i++)	{
			if (mArray[i] == 1)	{
				on++;
			}
			else if (mArray[i] == 0)	{
				off++;
			}
		}
		for (int i=0; i<Ntaxa; i++)	{
			if (on > off)	{
				if (mArray[i] == 0)	{
					os << mParam->SpeciesNames[i] << '\n';
				}
			}
			else	{
				if (mArray[i] == 1)	{
					os << mParam->SpeciesNames[i] << '\n';
				}
			}
		}
	}
	else	{
		for (int i=0; i<Ntaxa; i++)	{

			if (mArray[i] ==-1)	{
				os << ' ';
			}
			else if (mArray[i] == 1)	{
				os << '*';
			}
			else if (mArray[i] == 0)	{
				os << '.';
			}
			else	{
				cerr << "error in Bipartition::WriteToStream(ostream& os)\n";
				cerr << "found : " << mArray[i] << '\n'; 
			exit(1);
			}
		}
	}
}



// ---------------------------------------------------------------------------------
//		 ReadFromStream()
// ---------------------------------------------------------------------------------

void	Bipartition::ReadFromStream(istream& is)	{
	string temp;
	is >> temp;
	*this = temp;
}

