#include "Random.h"
#include "SequenceAlignment.h"
#include "Parallel.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])	{

	string name = argv[1];
	int nstep = atoi(argv[2]);
	int burnin = atoi(argv[3]);
	int size = atoi(argv[4]);


    ifstream is((name + ".stepping").c_str());
    ofstream os((name + ".poststep").c_str());

    double ess[nstep];
    double logl[nstep];
    int dn[nstep];
    double lnl[size];
    double f[nstep];

    double totlogl = 0;

    for (int step=0; step<nstep; step++)    {
        double frac, dlnl, dlnlpersite;
        int dnsite;
        for (int i=0; i<burnin; i++)    {
            is >> frac >> dnsite >> dlnl >> dlnlpersite;
        }
        double max = 0;
        for (int i=0; i<size; i++)  {
            is >> frac >> dnsite >> dlnl >> dlnlpersite;
            lnl[i] = dlnl;
            if ((!i) || (max < dlnl))   {
                max = dlnl;
            }
            if (!i) {
                dn[step] = dnsite;
                f[step] = frac;
            }
            else    {
                if (dn[step] != dnsite) {
                    cerr << "error: non matching number of sites\n";
                    exit(1);
                }
                if (f[step] != frac)  {
                    cerr << "error: non matching frac\n";
                    exit(1);
                }
            }
        }
        double tot = 0;
        for (int i=0; i<size; i++)  {
            tot += exp(lnl[i] - max);
        }
        double mean = tot/size;
        double meanlogl = log(mean) + max;
        logl[step] = meanlogl;
        totlogl += meanlogl;
        double e = 0;
        for (int i=0; i<size; i++)  {
            double w = exp(lnl[i] - max) / tot;
            e += w*w;
        }
        ess[step] = 1.0 / e;
        os << f[step] << '\t' << totlogl << '\t' << logl[step] << '\t' << logl[step] / dn[step] << '\t' << ess[step] << '\n';
    }
    ofstream sos((name + ".sumstep").c_str());
    sos << totlogl << '\n';
}

