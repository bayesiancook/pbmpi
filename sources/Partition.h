
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <string>

struct PartitionScheme {

public:

	PartitionScheme(int Nsite = 0) : Npart(0)
	{
		sitePart.resize(Nsite);
	};

	void update(){
		partSites.clear();
		partSites.resize(Npart);

		for(int i = 0; i < sitePart.size(); i++)
		{
			partSites[sitePart[i]].push_back(i);
		}
	}

	size_t GetNsite(){ return sitePart.size(); }

	int Npart;
	std::vector<int> sitePart;
	std::vector<std::vector<int> > partSites;
	std::vector<std::string> partType;
};

class PartitionProcess {

public:

	PartitionProcess() {};
	virtual ~PartitionProcess() {}

	int GetNpart() {return scheme.Npart;}

	int GetSitePart(int site) {return scheme.sitePart[site];}

	std::vector<int> GetPartSites(int part) {return scheme.partSites[part];}

	int GetPartNsite(int part) {return scheme.partSites[part].size();}

	std::string GetPartType(int part){ return scheme.partType[part]; }

	protected:

	void Create(PartitionScheme& inscheme)
	{
		scheme = inscheme;
	}
	void Delete()
	{

	}

	std::vector<PartitionScheme> ReadSchemes(std::string schemefile, int Nsite);

	PartitionScheme scheme;
};

#endif

