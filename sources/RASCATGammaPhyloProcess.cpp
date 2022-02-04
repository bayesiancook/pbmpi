
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"
#include <cassert>
#include "RASCATGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void RASCATGammaPhyloProcess::GlobalUpdateParameters()	{
	// MPI2
	// should send the slaves the relevant information
	// about model parameters

	// for this model, should broadcast
	// double alpha
	// int Ncomponent
	// int* alloc
	// double* rr
	// double** profile
	// double* brancharray
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	// SetAlpha(inalpha)

	assert(myid == 0);

	// ResampleWeights();
	RenormalizeProfiles();

	int i,j,nbranch = GetNbranch(),ni,nd,L1,L2;
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 3 + nbranch + L1*L2 + GetDim() + 1;
	ni = 1 + GetNsite();
	int ivector[ni];
	double dvector[nd];
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// First we assemble the vector of doubles for distribution
	int index = 0;
	dvector[index] = GetAlpha();
	index++;
	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	index++;
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}
	
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}
	dvector[index] = kappa;
	index++;

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<GetNsite(); ++i) {
		ivector[1+i] = DPProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

    case SITELOGLCUTOFF:
        SlaveSetSiteLogLCutoff();
        break;
	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
    /*
    case ISSITELOGL:
        SlaveComputeISSiteLogL();
        break;
    */
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


void RASCATGammaPhyloProcess::SlaveUpdateParameters()	{
	int i,j,L1,L2,ni,nd,nbranch = GetNbranch();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 3 + nbranch + L1*L2 + GetDim() + 1;
	ni = 1 + GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	SetAlpha(dvector[index]);
	index++;
	branchalpha = dvector[index];
	index++;
	branchbeta = dvector[index];
	index++;
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}
	kappa = dvector[index];
	index++;

	Ncomponent = ivector[0];
	for(i=0; i<GetNsite(); ++i) {
		DPProfileProcess::alloc[i] = ivector[1+i];
	}
	delete[] dvector;
	delete[] ivector;

	UpdateZip();
}

void RASCATGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	int ss = 0;
	int map = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic
	int cv = 0;
	int sitelogl = 0;
    int verbose = 0;
	int rates = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 1;
	int savetrees = 0;
    int clusters = 0;

	int ancstatepostprobs = 0;

    siteloglcutoff = 0;
    int posthyper = 0;
    int siteprofilesuffstat = 0;
    // importance sampling estimation of site logl
    /*
    int issitelogl = 0;
    int isnrep = 0;
    string empname = "None";
    */

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-div")	{
				ppred = 2;
			}
			else if (s == "-comp")	{
				ppred = 3;
			}
            else if (s == "-ncomp") {
                clusters = 1;
            }
            else if (s == "-siteconvprob")   {
                ppred = 4;
            }
            else if (s == "-sitecomp")   {
                ppred = 5;
            }
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-allppred")	{
				ppred = -1;
			}
			else if (s == "-ppredrate")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rateprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rateprior = 0;
				}
				else	{
					cerr << "error after ppredrate: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredprofile")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					profileprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					profileprior = 0;
				}
				else	{
					cerr << "error after ppredprofile: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredroot")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rootprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rootprior = 0;
				}
				else	{
					cerr << "error after ppredroot: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-savetrees")	{
				savetrees = 1;
			}

			else if (s == "-jointcv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}

			else if (s == "-sitecv")	{
				cv = 2;
				i++;
				testdatafile = argv[i];
			}

			else if (s == "-anc")	{
				ancstatepostprobs = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}

			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
            else if (s == "-v") {
                verbose =1 ;
            }
            else if (s == "-cutoff")    {
                i++;
                siteloglcutoff = atof(argv[i]);
            }
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-ss")	{
				ss = 1;
			}
            else if (s == "-posthyper") {
                posthyper = 1;
            }
            else if (s == "-siteprofilesuffstat")  {
                siteprofilesuffstat = 1;
            }
            /*
            else if (s == "-issitelogl")    {
                issitelogl = 1;
                i++;
                isnrep = atoi(argv[i]);
                i++;
                empname = argv[i];
            }
            */
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					string tmp = argv[i];
					if (IsInt(tmp))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else {
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}
	if (burnin == -1)	{
		burnin = GetSize() / 5;
	}

	if ((GetNprocs() == 1) && (ppred || cv || sitelogl))	{
		cerr << "error : should run readpb_mpi in mpi mode, with at least 2 processes\n";
		MPI_Finalize();
		exit(1);
	}

	if (ss)	{
		ReadSiteProfiles(name,burnin,every,until);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else if (cv == 1)	{
        GlobalSetSiteLogLCutoff();
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (cv == 2)	{
        GlobalSetSiteLogLCutoff();
		ReadSiteCV(testdatafile,name,burnin,every,until);
	}
	else if (ancstatepostprobs)	{
		ReadAncestral(name,burnin,every,until);
	}
    else if (clusters)  {
        ReadClusters(name, burnin, every, until);
    }
	else if (sitelogl)	{
        GlobalSetSiteLogLCutoff();
		ReadSiteLogL(name,burnin,every,until,verbose);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
    else if (posthyper) {
		ReadPostHyper(name,burnin,every,until);
    }
    else if (siteprofilesuffstat)  {
        ReadSiteProfileSuffStat(name, burnin, every, until);
    }
    /*
    else if (issitelogl)    {
        ReadISSiteLogL(name, empname, burnin, isnrep);
    }
    */
	else if (ppred == -1)	{
		AllPostPred(name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior,savetrees);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void RASCATGammaPhyloProcess::ReadPostHyper(string name, int burnin, int every, int until)	{

    double meanratealpha = 0;
    double varratealpha = 0;
    vector<double> meandirweight(GetDim(), 0);
    vector<double> vardirweight(GetDim(), 0);
    vector<double> meanbl(GetNbranch(), 0);
    vector<double> varbl(GetNbranch(), 0);
    double meankappaalpha = 0;
    double varkappaalpha = 0;

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

        meanratealpha += alpha;
        varratealpha += alpha*alpha;
        for (int k=0; k<GetDim(); k++)  {
            meandirweight[k] += dirweight[k];
            vardirweight[k] += dirweight[k]*dirweight[k];
        }
        for (int j=1; j<GetNbranch(); j++)  {
            meanbl[j] += blarray[j];
            varbl[j] += blarray[j]*blarray[j];
        }
        meankappaalpha += kappa;
        varkappaalpha += kappa*kappa;

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".posthyper").c_str());
    meanratealpha /= samplesize;
    varratealpha /= samplesize;
    varratealpha -= meanratealpha*meanratealpha;
    os << meanratealpha*meanratealpha / varratealpha << '\t';
    os << meanratealpha / varratealpha << '\n';
    os << '\n';
    if (! dirweightprior)   {
        for (int k=0; k<GetDim(); k++)  {
            meandirweight[k] /= samplesize;
            vardirweight[k] /= samplesize;
            vardirweight[k] -= meandirweight[k]*meandirweight[k];
            os << meandirweight[k]*meandirweight[k]/vardirweight[k] << '\t';
            os << meandirweight[k]/vardirweight[k] << '\n';
        }
        os << '\n';
    }
    for (int j=1; j<GetNbranch(); j++)  {
        meanbl[j] /= samplesize;
        varbl[j] /= samplesize;
        varbl[j] -= meanbl[j]*meanbl[j];
        os << meanbl[j]*meanbl[j]/varbl[j] << '\t';
        os << meanbl[j]/varbl[j] << '\n';
    }
    os << '\n';
    meankappaalpha /= samplesize;
    varkappaalpha /= samplesize;
    varkappaalpha -= meankappaalpha * meankappaalpha;
    os << meankappaalpha*meankappaalpha/varkappaalpha << '\t';
    os << meankappaalpha / varkappaalpha << '\n';
	cerr << "posterior mean shape and scale params in " << name << ".posthyper\n";
	cerr << '\n';
}

void RASCATGammaPhyloProcess::GlobalSetEmpiricalPrior(istream& is)    {
    // read from stream
    is >> empalpha >> empbeta;
    if (!dirweightprior)    {
        for (int k=0; k<GetDim(); k++)  {
            is >> empdirweightalpha[k] >> empdirweightbeta[k];
        }
    }
    for (int j=1; j<GetNbranch(); j++)  {
        is >> branchempalpha[j] >> branchempbeta[j];
    }
    is >> empkappaalpha >> empkappabeta;

	MESSAGE signal = EMPIRICALPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&empalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empbeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (!dirweightprior)    {
        MPI_Bcast(empdirweightalpha,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(empdirweightbeta,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
	MPI_Bcast(branchempalpha,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(branchempbeta,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empkappaalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empkappabeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*
    empcount = new double[GetNsite()*GetDim()];
    for (int k=0; k<GetNsite()*GetDim(); k++)   {
        is >> empcount[k];
    }
    MPI_Bcast(empcount,GetNsite()*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    */
}

void RASCATGammaPhyloProcess::SlaveSetEmpiricalPrior()    {

	MPI_Bcast(&empalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empbeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (!dirweightprior)    {
        MPI_Bcast(empdirweightalpha,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(empdirweightbeta,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
	MPI_Bcast(branchempalpha,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(branchempbeta,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empkappaalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empkappabeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*
    empcount = new double[GetNsite()*GetDim()];
    MPI_Bcast(empcount,GetNsite()*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    */
}

void RASCATGammaPhyloProcess::GlobalSetSiteLogLCutoff()  {

	MESSAGE signal = SITELOGLCUTOFF;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&siteloglcutoff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATGammaPhyloProcess::SlaveSetSiteLogLCutoff()  {
	MPI_Bcast(&siteloglcutoff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATGammaPhyloProcess::ReadSiteProfileSuffStat(string name, int burnin, int every, int until){

    double* meancount = new double[GetNsite()*GetDim()];
    for (int k=0; k<GetNsite()*GetDim(); k++)   {
        meancount[k] = 0;
    }

  	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}
	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";

		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalCollapse();

        GlobalUpdateSiteProfileSuffStat();
        for (int k=0; k<GetNsite()*GetDim(); k++)   {
            meancount[k] += allocsiteprofilesuffstatcount[k];
        }

		GlobalUnfold();

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
    for (int k=0; k<GetNsite()*GetDim(); k++)   {
        meancount[k] /= samplesize;
    }
    ofstream os((name + ".siteprofilesuffstat").c_str());
    for (int i=0; i<GetNsite(); i++)    {
        for (int k=0; k<GetDim(); k++)  {
            os << meancount[i*GetDim() + k] << '\t';
        }
        os << '\n';
    }
    cerr << "posterior mean site profile suffstats in " << name << ".siteprofilesuffstat\n";
    cerr << '\n';
}

/*
void RASCATGammaPhyloProcess::ReadISSiteLogL(string name, string empname, int burnin, int nrep)    {

    double* empcount = new double[GetNsite()*GetDim()];
    ifstream eis(empname.c_str());
    for (int k=0; k<GetNsite()*GetDim(); k++)   {
        eis >> empcount[k];
    }

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	int maxwidth = 0;
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
		if (maxwidth < (smax[i] - smin[i]))	{
			maxwidth = smax[i] - smin[i];
		}
	}

	int i=0;
    cerr << "burnin\n";
	while (i < burnin)	{
        cerr << '.';
		FromStream(is);
		i++;
	}
    cerr << '\n';

    cerr << "computing IS site logls\n";
    QuickUpdate();
    MPI_Status stat;
    MESSAGE signal = ISSITELOGL;
    MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(empcount,GetNsite()*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	double* meansitelogl = new double[GetNsite()];
	double* varsitelogl = new double[GetNsite()];

    for(int i=1; i<GetNprocs(); ++i) {
        MPI_Recv(meansitelogl,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        MPI_Recv(varsitelogl,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
    }

    cerr << "ok\n";

    ofstream os((name + ".priorsitelogl").c_str());
    for (int i=0; i<GetNsite(); i++)    {
        os << meansitelogl[i] << '\t' << varsitelogl[i] << '\n';
    }
    cerr << "prior site logls in " << name << ".priorsitelogl\n";
    cerr << '\n';

    for(int i=1; i<GetNprocs(); ++i) {
        MPI_Recv(meansitelogl,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        MPI_Recv(varsitelogl,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
    }

    cerr << "ok\n";

    ofstream pos((name + ".postsitelogl").c_str());
    for (int i=0; i<GetNsite(); i++)    {
        pos << meansitelogl[i] << '\t' << varsitelogl[i] << '\n';
    }
    cerr << "post site logls in " << name << ".postsitelogl\n";
    cerr << '\n';

    delete[] meansitelogl;
    delete[] varsitelogl;
}

void RASCATGammaPhyloProcess::SlaveComputeISSiteLogL()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

    if (GetNmodeMax() < GetNsite()) {
        cerr << "error: nmodemax should be equal to nsite\n";
        exit(1);
    }

    for (int i=sitemin; i<sitemax; i++)	{
        PoissonDPProfileProcess::alloc[i] = i;
    }
	
    int nrep;
    MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

    double* empcount = new double[GetNsite()*GetDim()];
    MPI_Bcast(empcount,GetNsite()*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            sitelogl[i] = new double[nrep];
        }
	}

	for (int k=0; k<nrep; k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                SampleStat(profile[i]);
                UpdateZip(i);
            }
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                sitelogl[i][k] = sitelogL[i];
            }
		}
	}

	double* priormeansitelogl = new double[GetNsite()];
	double* priorvarsitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		priormeansitelogl[i] = 0;
		priorvarsitelogl[i] = 0;
	}

	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            double max = 0;
            for (int k=0; k<nrep; k++) {
                if ((!k) || (max < sitelogl[i][k]))	{
                    max = sitelogl[i][k];
                }
            }
            double tot = 0;
            for (int k=0; k<nrep; k++)	{
                tot += exp(sitelogl[i][k] - max);
            }
            tot /= nrep;
            priormeansitelogl[i] = log(tot) + max;

            double relvar = 0;
            for (int k=0; k<nrep; k++)	{
                double tmp = exp(sitelogl[i][k] - max) / tot;
                relvar += tmp*tmp;
            }
            relvar /= nrep;
            relvar -= 1;
            priorvarsitelogl[i] = relvar;
        }
    }

	MPI_Send(priormeansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(priorvarsitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	delete[] priormeansitelogl;
	delete[] priorvarsitelogl;

    // Empirical Posterior

	for (int k=0; k<nrep; k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                SampleEmpiricalStat(profile[i], empcount + i*GetDim());
                UpdateZip(i);
                sitelogl[i][k] = LogStatPrior(profile[i]) - EmpiricalLogStatPrior(profile[i], empcount + i*GetDim());
            }
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                sitelogl[i][k] += sitelogL[i];
            }
		}
	}

	double* postmeansitelogl = new double[GetNsite()];
	double* postvarsitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		postmeansitelogl[i] = 0;
		postvarsitelogl[i] = 0;
	}

	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            double max = 0;
            for (int k=0; k<nrep; k++) {
                if ((!k) || (max < sitelogl[i][k]))	{
                    max = sitelogl[i][k];
                }
            }
            double tot = 0;
            for (int k=0; k<nrep; k++)	{
                tot += exp(sitelogl[i][k] - max);
            }
            tot /= nrep;
            postmeansitelogl[i] = log(tot) + max;

            double relvar = 0;
            for (int k=0; k<nrep; k++)	{
                double tmp = exp(sitelogl[i][k] - max) / tot;
                relvar += tmp*tmp;
            }
            relvar /= nrep;
            relvar -= 1;
            postvarsitelogl[i] = relvar;
        }
    }

	MPI_Send(postmeansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(postvarsitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	delete[] postmeansitelogl;
	delete[] postvarsitelogl;

}
*/

void RASCATGammaPhyloProcess::ReadSiteProfiles(string name, int burnin, int every, int until)	{

	double** sitestat = new double*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		sitestat[i] = new double[GetDim()];
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] = 0;
		}
	}
	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int i=0; i<GetNsite(); i++)	{
			double* p = GetProfile(i);
			for (int k=0; k<GetDim(); k++)	{
				sitestat[i][k] += p[k];
			}
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".siteprofiles").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetNsite(); i++)	{
		os << i + 1;
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] /= samplesize;
			os << '\t' << sitestat[i][k];
		}
		os << '\n';
	}

	cerr << "mean site-specific profiles in " << name << ".siteprofiles\n";
	cerr << '\n';
}

void RASCATGammaPhyloProcess::ReadClusters(string name, int burnin, int every, int until)	{

    double meanncomp = 0;
    double varncomp = 0;
    double meaneffncomp = 0;
    double vareffncomp = 0;

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

    ofstream kos((name + ".kappa").c_str());

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

        kos << kappa << '\n';
        UpdateOccupancyNumbers();
        double k = GetNOccupiedComponent();
        double keff = GetEffectiveComponentNumber();
        meanncomp += k;
        varncomp += k*k;
        meaneffncomp += keff;
        vareffncomp += keff*keff;
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
    meanncomp /= samplesize;
    varncomp /= samplesize;
    varncomp -= meanncomp*meanncomp;
    meaneffncomp /= samplesize;
    vareffncomp /= samplesize;
    vareffncomp -= meaneffncomp*meaneffncomp;

	ofstream os((name + ".effncomp").c_str());
    os << "Ncomp             : " << meanncomp << " +/- " << sqrt(varncomp) << '\n';
    os << "effective Ncomp   : " << meaneffncomp << " +/- " << sqrt(vareffncomp) << '\n';
    cerr << "Ncomp             : " << meanncomp << " +/- " << sqrt(varncomp) << '\n';
    cerr << "effective Ncomp   : " << meaneffncomp << " +/- " << sqrt(vareffncomp) << '\n';
}
