
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

class PolyNode;
class BipartitionList;
class Bipartition;
class TreeList;

class Consensus		:   public PBTree	{

	public:

	Consensus(TreeList* inTreeList, double* probarray = 0, double cutoff = 0.5);
	Consensus(BipartitionList* inBList, double cutoff = 0.5);

	~Consensus();

	BipartitionList*				MakeConsensus(double cutoff);

	BipartitionList*				GetBList()		const   {return mBList;}
	TreeList*					GetTreeList()	const {return mTreeList;}
	
	void						Insert( Bipartition& inPartition, double prob, double length);
	TreeList*					mTreeList;
	BipartitionList*				mBList;
	BipartitionList*				mReducedList;

}
;
