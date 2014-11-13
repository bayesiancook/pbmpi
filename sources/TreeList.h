
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
class PBTree;
class PolyNode;
class Bipartition;
class PhyloBayes;

class TreeList	{


	public:

				TreeList(TaxaParameters* inParam, int inSize);
				TreeList();
				TreeList(string filename);
				~TreeList();

	PBTree*			GetTree(int index);

	int			GetSize()	{return mSize;}

	TaxaParameters*		GetParameters()	const {return mParam;}
	void			SetParameters(); // set parameters for each tree of the list 

	void			RootAt(Bipartition outgroup);

	void			ReadFromStream(istream& is);
	void			WriteToStream(ostream& os, int header = 0, int withLengths = 1, int withProbs = 0, int withSpeciesNames = 1,int withInternalLabels = 0);
	void			ToPS(string target, int every = 1, double sizeX=12, double sizeY=20, int withLengths=0, int withProbs = 0, int withSpeciesNames = 1, int withInternalLabels = 0);

	TaxaParameters*		mParam;
	int 			mSize;
	PBTree**		mTreeArray;

}
;


