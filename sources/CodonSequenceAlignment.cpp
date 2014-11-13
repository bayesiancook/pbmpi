
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

#include "CodonSequenceAlignment.h"
// #include "Exception.h"

CodonSequenceAlignment::CodonSequenceAlignment(SequenceAlignment* from, bool force_stops, GeneticCodeType type)	{

		try	{
			DNAsource = from;

			if (from->Nsite % 3 != 0)	{
				cerr << "not multiple of three\n";
				exit(1);
			}		 
			Nsite = from->Nsite/3;
			Ntaxa = from->Ntaxa;
			statespace = new CodonStateSpace(type);

			taxset = DNAsource->GetTaxonSet();

			// make my own arrays
			// make translation
			Data = new int*[Ntaxa];
			for (int i=0; i<Ntaxa; i++)	{
				Data[i] = new int[Nsite];
				for (int j=0; j<Nsite; j++)	{
					try {
						Data[i][j] = GetCodonStateSpace()->GetCodonFromDNA(DNAsource->GetState(i, 3*j), DNAsource->GetState(i, 3*j+1), DNAsource->GetState(i, 3*j+2));
					}
					catch(...)	{
					// catch(Exception e)	{
						cerr << "in CodonSequenceAlignment: taxon " << i << " and codon " << j << " (site " << 3*j << ")\n";
						cerr << "taxon : " << taxset->GetTaxon(i) << '\n'; 
						if (force_stops)	{
							Data[i][j] = -1;
						}
						else	{
							throw;
						}
					}
				}
			}

		}
		catch(...)	{
		// catch(Exception)	{
			cerr << "Codon Sequence Alignment: failed to read the datafile\n";
			exit(1);
		}
	}


void CodonSequenceAlignment::ToStream(ostream& os)	{

	os << Ntaxa << '\t' << 3 * Nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}
	
	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			os << statespace->GetState(GetState(i,j));
		}
		os << '\n';
	}
	os << '\n';
}

/*void CodonSequenceAlignment::ToStream(ostream& os)	{

	os << Ntaxa << '\t' << Nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}
	
	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			int tempCodonState = Data[i][j];
			int tempAAState = GetCodonStateSpace()->Translation(tempCodonState);
			if (tempCodonState == unknown)	{
				os << '-';
			}
			else {
				os << AminoAcids[tempAAState];
			}
		}
		os << '\n';
	}
	os << '\n';
}*/
