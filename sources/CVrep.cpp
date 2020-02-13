#include "Random.h"
#include "SequenceAlignment.h"
#include "Parallel.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])	{

	if (argc != 6)	{
		cerr << "cvrep datafile N1 N2 Nrep outfile\n";
		exit(1);
	}

	// Random::Random();
	string datafile = argv[1];
	int N1 = atoi(argv[2]);
	int N2 = atoi(argv[3]);
	int Nrep = atoi(argv[4]);
	string outfile = argv[5];
	
	SequenceAlignment* data = new FileSequenceAlignment(datafile,1,0);
	
	if ((N1+N2) > data->GetNsite())	{
		cerr << "error : target size should be < Nsite\n";
		exit(1);
	}

	int* mask = new int[data->GetNsite()];
	int* index = new int[N1+N2];
	
	for (int rep=0; rep<Nrep; rep++)	{
		for (int i=0; i<data->GetNsite(); i++)	{
			mask[i] = 0;
		}
		rnd::GetRandom().DrawFromUrn(index,N1+N2,data->GetNsite());
		for (int i=0; i<N1; i++)	{
			mask[index[i]] = 1;
		}
		ostringstream s;
		s << outfile << rep << "_learn.ali";
		ofstream os(s.str().c_str());
		SequenceAlignment* datalearn = new SequenceAlignment(data,mask);
		datalearn->ToStream(os);
		os.close();

		for (int i=0; i<data->GetNsite(); i++)	{
			mask[i] = 0;
		}
		for (int i=N1; i<N1+N2; i++)	{
			mask[index[i]] = 1;
		}
		ostringstream s2;
		s2 << outfile << rep << "_test.ali";
		ofstream os2(s2.str().c_str());
		SequenceAlignment* datatrain= new SequenceAlignment(data,mask);
		datatrain->ToStream(os2);
		os2.close();
	}

}

