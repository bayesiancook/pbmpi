
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

class TaxaParameters;

class BipartitionList; 

class Bipartition	{

	public:

	Bipartition(TaxaParameters* inParam);
	// empty by default
	Bipartition(const Bipartition& from);
	~Bipartition();

	Bipartition& 			operator=(const Bipartition& from);
	Bipartition& 			operator=(const string& from);

	Boolean				operator==( const Bipartition& inPartition); // no modulo !
	Boolean				operator!=( const Bipartition& inPartition); // no modulo !

	Bipartition&			operator|=( const Bipartition& inPartition); // assumes they are oriented likewise
	Bipartition&			operator&=( const Bipartition& inPartition);
	Bipartition			operator!();
	
	Boolean 			IsCompatibleWith( const Bipartition& inPartition);
	Boolean 			IsCompatibleWith( const Bipartition& inPartition, Boolean& Orientation);
	Boolean				IsInformative();

	int				CompareWith(const Bipartition& with);
	int				CompareWith(BipartitionList* bplist);
	int				SupportCheck(const Bipartition& with);
	int				SupportCheck(BipartitionList* bplist);

	void				AllAbsent();
	void				Suppress(const Bipartition& leafset);

	void				PermutTaxa(int* permut);
	int				GetTaxonStatus(int index) ;
	void				SetTaxon(int index);

	TaxaParameters*			GetParameters()	const;
	
	void				WriteToStream(ostream& os, int verbose = 0);
	void				ReadFromStream(istream& is);
	double				GetPriorProb();
	
	void				Modulo();
	
	int*				mArray;
			// -1 : species eliminated
			// 0 : species upstream
			// 1	: species downstream

	TaxaParameters* 		mParam;
	int				Ntaxa;
	
}
;
