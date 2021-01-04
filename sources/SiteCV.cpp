#include "Random.h"
#include "SequenceAlignment.h"
#include "Parallel.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])	{

	string name = argv[1];
	int testnsite = atoi(argv[2]);
	int samplesize = atoi(argv[3]);

    cerr << name << '\t' << testnsite << '\t' << samplesize << '\n';

    vector<vector<double>> logl(testnsite, vector<double>(samplesize,0));
    ifstream is((name + ".cvsitelogl").c_str());
    for (int j=0; j<samplesize; j++)    {
        for (int i=0; i<testnsite; i++) {
            is >> logl[i][j];
        }
    }

    /*
    vector<double> tmpscore(samplesize, 0);
    vector<double> x(samplesize, 0);

    ofstream kos((name + ".ksitecv").c_str());

    // cerr << "#k\tcvscore\tbias\tmeanlogl\tvar\tstdev\tstdev2\tstdev3\tess\tminess\n";
    // kos << "#k\tcvscore\tbias\tmeanlogl\tvar\tstdev\tstdev2\tstdev3\tess\tminess\n";

    for (int k=1; k<=testnsite; k+=10) {

        double max = 0;
        double meanlog = 0;
        double varlog = 0;
        for (int j=0; j<samplesize; j++)	{
            tmpscore[j] = 0;
            for (int i=0; i<k; i++) {
                tmpscore[j] += logl[i][j];
            }
            if ((!j) || (max < tmpscore[j]))    {
                max = tmpscore[j];
            }
            meanlog += tmpscore[j];
            varlog += tmpscore[j]*tmpscore[j];
        }
        meanlog /= samplesize;
        varlog /= samplesize;
        varlog -= meanlog*meanlog;

        double tot = 0;
        for (int j=0; j<samplesize; j++)	{
            x[j] = exp(tmpscore[j] - max);
            tot += x[j];
        }

        double invess = 0;
        for (int j=0; j<samplesize; j++)	{
            double w = x[j] / tot;
            invess += w*w;
        }
        double ess = 1.0/invess;

        tot /= samplesize;

        double cvscore = log(tot) + max;
        
        double v = 0;
        for (int j=0; j<samplesize; j++)	{
            double z = x[j] / tot;
            v += z*z;
        }
        v /= samplesize;
        v -= 1.0;
        double bias = v / 2 / samplesize;
        double var1 = v / samplesize;
        // double var2 = exp(2*(meanlog-cvscore)) * varlog / samplesize;
        double var2 = varlog / samplesize;

        cerr << k << '\t' << cvscore << '\t' << bias << '\t' << sqrt(var1) << '\t' << sqrt(var2) << '\t' << meanlog << '\t' << varlog << '\t' << ess << '\n';
        kos << k << '\t' << cvscore << '\t' << bias << '\t' << sqrt(var1) << '\t' << sqrt(var2) << '\t' << meanlog << '\t' << varlog << '\t' << ess << '\n';
    }
    */

    int K = testnsite;
    int nrepmax = 1000*testnsite;
    vector<double> tmpscore(samplesize, 0);
    vector<double> x(samplesize, 0);

    ofstream kos((name + ".ksitecv").c_str());

    // cerr << "#k\tcvscore\tbias\tmeanlogl\tvar\tstdev\tstdev2\tstdev3\tess\tminess\n";
    // kos << "#k\tcvscore\tbias\tmeanlogl\tvar\tstdev\tstdev2\tstdev3\tess\tminess\n";

    int step = 1;
    for (int k=1; k<=K; k+=step) {

        int sites[k];
        double kmeancvscore = 0;
        double kmeanbias = 0;
        double kmeanvar1 = 0;
        double kmeanvar2 = 0;
        double kmeaness = 0;
        double kminess = 0;
        double kmeanlog = 0;
        double kvarlog = 0;

        int nrep = nrepmax / k;
        for (int rep=0; rep<nrep; rep++)    {
            rnd::GetRandom().DrawFromUrn(sites,k,testnsite);
            double max = 0;
            double meanlog = 0;
            double varlog = 0;
            for (int j=0; j<samplesize; j++)	{
                tmpscore[j] = 0;
                for (int i=0; i<k; i++) {
                    tmpscore[j] += logl[sites[i]][j];
                }
                if ((!j) || (max < tmpscore[j]))    {
                    max = tmpscore[j];
                }
                meanlog += tmpscore[j];
                varlog += tmpscore[j]*tmpscore[j];
            }
            meanlog /= samplesize;
            varlog /= samplesize;
            varlog -= meanlog*meanlog;
            // cerr << meanlog << " +/- " << varlog << '\t';
            kmeanlog += meanlog;
            kvarlog += varlog;

            double tot = 0;
            for (int j=0; j<samplesize; j++)	{
                x[j] = exp(tmpscore[j] - max);
                tot += x[j];
            }

            double invess = 0;
            for (int j=0; j<samplesize; j++)	{
                double w = x[j] / tot;
                invess += w*w;
            }
            double ess = 1.0/invess;
            kmeaness += ess;
            if ((!rep) || (kminess > ess))  {
                kminess = ess;
            }
            tot /= samplesize;
            
            double v = 0;
            for (int j=0; j<samplesize; j++)	{
                double z = x[j] / tot;
                v += z*z;
            }
            v /= samplesize;
            v -= 1.0;

            double bias = v / 2 / samplesize;
            kmeanbias += bias;
            double var1 = v / samplesize;
            kmeanvar1 += var1;
            double var2 = (exp(varlog) - 1) / samplesize;
            // cerr << varlog << '\t' << var2 << '\n';
            // double var2 = 1.0 / varlog / samplesize * k * k;
            // double var2 = 1.0 / (varlog-2) / samplesize;
            kmeanvar2 += var2;

            double cvscore = (log(tot) + max)/k;
            kmeancvscore += cvscore;
        }
        kmeancvscore /= nrep;
        kmeaness /= nrep;
        kmeanbias /= nrep;
        kmeanlog /= nrep;
        kvarlog /= nrep;
        kmeanvar1 /= nrep;
        kmeanvar2 /= nrep;

        cerr << k << '\t' << kmeancvscore << '\t' << kmeanbias/k << '\t' << sqrt(kmeanvar1)/k << '\t' << sqrt(kmeanvar1) << '\t' << kmeanlog << '\t' << kvarlog << '\t' << kmeaness << '\t' << kminess << '\n';
        if (k > 100)    {
            step = 100;
        }
        else if (k > 10)    {
            step = 10;
        }
    }
}
