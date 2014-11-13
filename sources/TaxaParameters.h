
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


class MCParameters;
class Bipartition;
class PBTree;

class TaxaParameters	{

	public:

	TaxaParameters();
	TaxaParameters(string filename);
	TaxaParameters(int N, string* inNames = 0);
	TaxaParameters(MCParameters* inParam);
	TaxaParameters(PBTree* inTree);
	TaxaParameters(const TaxaParameters& from);
	~TaxaParameters();

	void				ReadFromStream(istream& is);
	void				WriteToStream(ostream& os);
	void				ReadNexus(istream& is);

	
	string	GetSpeciesName(int index);
	Boolean Member(int index);

	int		RegisterWith(TaxaParameters* inParam);
	MCParameters*	GetParameters()	{return mParam;}

	Bipartition			GetOutGroup();
	void				SetOutGroup(Bipartition& inOutGroup);

	MCParameters*			mParam;
	int				Ntaxa;				// internal number of taxa
	string*				SpeciesNames;
		
}
;
