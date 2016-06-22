

/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "Parallel.h"
#include "BP2Stat.h"
#include <vector>
#include "StringStreamUtils.h"
MPI_Datatype Propagate_arg;

const double default_cutoff = 0.05;
const int MaxChain = 10;

int main(int argc, char* argv[])	{

	string* filenames = new string[argc];
	for (int i=0; i<argc; i++)	{
		filenames[i] = "";
	}

	string OutFile = "";
	double cutoff = 0.05;
	int burnin = -1;
	double conscutoff = 0.5;
	int verbose = 0;


	double mixstatdown = 0.5;
	double mixstatup = 0.9;

	bool bench = false;

	if (argc == 1)	{
		cerr << "bpcomp [-cox] ChainName1 ChainName2 ... \n";
		cerr << "\t-c <cutoff> : only partitions with max prob >  cutoff. (default 0.05)\n";
		cerr << "\t-x <burnin> [<every> <until>]. default burnin = 0\n";
		cerr << "\t-i <down> <up>. default interval for mix statistic is 0.2 - 0.8\n";
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

	while (i < argc)	{
		string s = argv[i];
		if (s == "-i")	{
			i++;
			mixstatdown = atof(argv[i]);
			i++;
			mixstatup = atof(argv[i]);
		}
		else if (s == "-c")	{
			i++;
			s = argv[i];
			conscutoff = atof(argv[i]);
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
			verbose = 1;
		}
		else if (s == "-o")	{
			i++;
			OutFile = argv[i];
		}
		else	{
			filenames[P] = argv[i];
			P++;
		}
		i++;
	}

	TaxaParametersBis* taxaparameters;

	if (! ifstream(filenames[0].c_str()))	{
		taxaparameters = new TaxaParametersBis(filenames[0] + ".treelist");
	}
	else{
		taxaparameters = new TaxaParametersBis(filenames[0]);
	}

	for(int i=0; i<P; i++){
	    	std::vector <char> buff( 1024*1024 );
	    	ifstream ifs( filenames[i].c_str() );
	    	int n = 0, c = 1;
		while(c){
			ifs.read( &buff[0], 1024*1024 );
			c=ifs.gcount();
			const char* p = &buff[0];
			for ( int j = 0; j < c; j++ ) {
				if ( p[j] == '\n' ) {
			    		n++;
				}
			}
		}
		if( (until>n) or (until == -1) ){until=n;}
		ifs.close();
	}
	cout << "The smaller chain have " << until << " points\n";

	BipartitionList** list = new BipartitionList*[P];
	for(int i=0; i<P; i++){
		list[i] = new BipartitionList(filenames[i], taxaparameters, burnin, every, until, verbose, mixstatdown, mixstatup);
	}
	BipartitionCompare* bpc = new BipartitionCompare(list, P, 0.05, OutFile);


}

