
#include "Random.h"
#include "SequenceAlignment.h"
#include "Parallel.h"

int main(int argc, char* argv[])	{


	string datafile = argv[1];
	string partfile = argv[2];
	string contfile = argv[3];
	double minfrac = atof(argv[4]);
	string outfile = argv[5];

	SequenceAlignment* ali = new FileSequenceAlignment(datafile,0,0);

	ifstream is(argv[2]);

	int Ngene;
	is >> Ngene;
	int* genesize = new int[Ngene];
	int totsize = 0;
	for (int i=0; i<Ngene; i++)	{
		int tmp;
		int begin, end;
		is >> tmp >> begin >> end;
		genesize[i] = end - begin + 1;	
		totsize += genesize[i];
	}

	if (totsize != ali->GetNsite())	{
		cerr << "error: non matching total sequence length\n";
		cerr << totsize << '\t' << ali->GetNsite() << '\n';
		exit(1);
	}

	int* exclude = new int[Ngene];
	for (int i=0; i<Ngene; i++)	{
		exclude[i] = 0;
	}

	ifstream eis(argv[3]);

	int nex;
	eis >> nex;
	for (int i=0; i<nex; i++)	{
		int tmp;
		eis >> tmp;
		exclude[tmp-1] = 1;
	}

	double* frac = new double[Ngene];

	ali->MissingFractionPerGene(Ngene,genesize,exclude,frac);

	for (int i=0; i<Ngene; i++)	{
		cout << genesize[i] << '\t' << frac[i] << '\n';
	}

	SequenceAlignment* subali = new SequenceAlignment(ali,Ngene,genesize,exclude,frac,minfrac);
	ofstream os(outfile.c_str());
	subali->ToStream(os);
	
}


