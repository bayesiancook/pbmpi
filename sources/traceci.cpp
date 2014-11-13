
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "Random.h"
#include "TexTab.h"

int main(int argc, char* argv[])	{

	int burnin = atoi(argv[1]);
	ifstream is(argv[2]);
	int K = atoi(argv[3]);

	list<double>* listarray = new list<double>[K];

	for (int k=0; k<K; k++)	{
		string temp;
		is >> temp;
		cerr << temp << '\n';
	}
	int count = 0;
	while (!is.eof())	{
		for (int k=0; k<K; k++)	{
			double tmp;
			is >> tmp;
			if (! is.eof())	{
			listarray[k].push_back(tmp);
			if (count >= 22000)	{
				cerr << "TMP:" << tmp << "::" << '\n';
				cerr << count << '\n';
				exit(1);
			}
			}
		}
		if (! is.eof())	{
		count++;
		}
	}
	cerr << "total : " << count << '\n';

	for (int k=3; k<K; k++)	{
		cout << textabentry(listarray[k],true,false,true) << '\n';
	}
}

