#include "Random.h"
#include "SequenceAlignment.h"
#include "Parallel.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])	{

	string name = argv[1];
	int step = atoi(argv[2]);
	int burnin = atoi(argv[3]);
	int size = atoi(argv[4]);
    int nsite = atoi(argv[5]);
    double kfold = atof(argv[6]);
    int nstep = nsite / step;
    if (nsite % step)   {
        nstep++;
    }

    ifstream is((name + ".stepping").c_str());
    ofstream os((name + ".poststep").c_str());

    double ess[nstep];
    double logl[nstep];
    double meanlog[nstep];
    double varlog[nstep];
    int dn[nstep];
    double lnl[size];
    int totnsite[nstep];

    double totlogl = 0;

    for (int step=0; step<nstep; step++)    {
        double dlnl, dlnlpersite;
        int n, dnsite;
        for (int i=0; i<burnin; i++)    {
            is >> n >> dnsite >> dlnl >> dlnlpersite;
        }
        double max = 0;
        for (int i=0; i<size; i++)  {
            is >> n >> dnsite >> dlnl >> dlnlpersite;
            lnl[i] = dlnl;
            if ((!i) || (max < dlnl))   {
                max = dlnl;
            }
            dn[step] = dnsite;
            totnsite[step] = n;
            /*
            if (!i) {
                dn[step] = dnsite;
                totnsite[step] = n;
            }
            else    {
                if (dn[step] != dnsite) {
                    cerr << "error: non matching number of sites\n";
                    exit(1);
                }
                if (totnsite[step] != n)  {
                    cerr << "error: non matching step\n";
                    exit(1);
                }
            }
            */
        }
        meanlog[step] = 0;
        varlog[step] = 0;
        for (int i=0; i<size; i++)  {
            meanlog[step] += lnl[i];
            varlog[step] += lnl[i] * lnl[i];
        }
        meanlog[step] /= size;
        varlog[step] /= size;
        varlog[step] -= meanlog[step] * meanlog[step];

        double tot = 0;
        for (int i=0; i<size; i++)  {
            tot += exp(lnl[i] - max);
        }
        double mean = tot/size;
        double meanlogl = log(mean) + max;
        logl[step] = meanlogl;
        if (step >= (1-kfold)*nstep) {
            totlogl += meanlogl;
        }
        double e = 0;
        for (int i=0; i<size; i++)  {
            double w = exp(lnl[i] - max) / tot;
            e += w*w;
        }
        ess[step] = 1.0 / e;
        os << totnsite[step] << '\t' << totlogl << '\t' << logl[step] << '\t' << logl[step] / dn[step] << '\t' << meanlog[step] << '\t' << varlog[step] << '\t' << ess[step] << '\n';
    }
    ofstream sos((name + ".sumstep").c_str());
    sos << totlogl << '\n';
}

