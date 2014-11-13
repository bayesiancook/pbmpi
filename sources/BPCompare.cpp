
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
#include "Parallel.h"
MPI_Datatype Propagate_arg;

const double default_cutoff = 0.05;
const int MaxChain = 10;

int main(int argc, char* argv[])	{

	int mergeallbp = 0;
	string* ChainName = new string[argc];
	for (int i=0; i<argc; i++)	{
		ChainName[i] = "";
	}

	string OutFile = "";
	double cutoff = 0.05;
	int burnin = -1;
	double conscutoff = 0.5;
	int verbose = 1;

	bool bench = false;

	if (argc == 1)	{
		cerr << "bpcomp [-cox] ChainName1 ChainName2 ... \n";
		cerr << "\t-c <cutoff> : only partitions with max prob >  cutoff. (default 0.5)\n";
		cerr << "\t-o <output> : detailed output into file\n"; 
		cerr << "\t-ps         : postscript output (requires LateX)\n";
		cerr << "\t-x <burnin> [<every> <until>]. default burnin = 1/10 of the chain\n";
		cerr << '\n';
		cerr << "\t compare bipartition frequencies between independent chains\n";
		cerr << "\t and build consensus based on merged lists of trees\n";
		cerr << '\n';
		exit(1);
	}


	int i = 1;
	int P = 0; // number of chains to be compared

	int every = 1;
	int until = -1;
	int ps = 0;

	bool rootonly = false;

	while (i < argc)	{
		string s = argv[i];
		if (s == "-m")	{
			mergeallbp = 1;
		}
		else if (s == "-c")	{
			i++;
			s = argv[i];
			conscutoff = atof(argv[i]);
		}
		else if (s == "-r")	{
			rootonly = true;
		}
		else if ( (s == "-x") || (s == "-extract") )	{
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (! IsInt(s))	{
				throw(0);
			}
			burnin = atoi(argv[i]);
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (IsInt(s))	{
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					until = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else	{
				i--;
			}
		}
		else if (s == "-v")	{
			verbose = 2;
		}
		else if (s == "-o")	{
			i++;
			OutFile = argv[i];
		}
		else	{
			ChainName[P] = argv[i];
			P++;
		}
		i++;
	}

	BPCompare(ChainName, P, burnin, every, until, ps, verbose, mergeallbp, OutFile, cutoff, conscutoff, rootonly, bench);

}

