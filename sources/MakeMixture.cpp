
#include "phylo.h"

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	ifstream pis(argv[2]);
	double epsilon = atof(argv[3]);
	ofstream os(argv[4]);

	string alphabet;
	is >> alphabet;
	int Nstate = alphabet.size();
	cerr << "nstate : " << Nstate << '\n';

	double stat[Nstate];
	for (int i=0; i<Nstate; i++)	{
		is >> stat[i];
	}

	int Ncat;
	pis >> Ncat;
	os << Nstate;
	for (int i=0; i<Nstate; i++)	{
		os << alphabet[i] << '\t';
	}
	os << Ncat +1 << '\n';
	for (int k=0; k<Ncat; k++)	{

		string pattern;
		pis >> pattern;
		int n = pattern.size();
		int pat[Nstate];
		for (int i=0; i<Nstate; i++)	{
			pat[i] = 0;
		}
		for (int i=0; i<n; i++)	{
			char c = pattern[i];
			for (int j=0; j<Nstate; j++)	{
				if (c == alphabet[j])	{
					pat[j] = 1;
				}
			}
		}

		double catstat[Nstate];
		double tot = 0;
		for (int i=0; i<Nstate; i++)	{
			if (pat[i])	{
				catstat[i] = stat[i];
			}
			else	{
				catstat[i] = stat[i] * epsilon * n / Nstate;
			}
			tot += catstat[i];
		}
		os << 1.0;
		for (int i=0; i<Nstate; i++)	{
			catstat[i] /= tot;
			os << '\t' << catstat[i];
		}
		os << '\n';
	}
	os << 1.0;
	for (int i=0; i<Nstate; i++)	{
		os << '\t' << stat[i];
	}
	os << '\n';
}

