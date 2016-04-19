
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include <cassert>
#include "PartitionedRASGTRGammaPhyloProcess.h"
#include "Parallel.h"
#include "StringStreamUtils.h"
#include <string>
#include <list>

void PartitionedRASGTRGammaPhyloProcess::GlobalUpdateParameters()	{
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

	int i,j,nrr,nalpha,nbranch = GetNbranch(),nd,ndim,nrrpart,nstatpart;
	nalpha = PartitionedDGamRateProcess::GetNpart();
	nrr = GetNrr();
	nrrpart = PartitionedGTRProfileProcess::GetNpart();
	ndim = GetDim();
	nstatpart = PartitionedProfileProcess::GetNpart();
	nd = 2*nalpha + 2 + nbranch + nfreerr*nrr + nfreestat*ndim + ndim;
	double* dvector = new double[nd];
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// GlobalBroadcastTree();

	// First we assemble the vector of doubles for distribution
	int index = 0;
	for(i = 0; i < nalpha; i++)
	{
		dvector[index] = GetAlpha(i);
		index++;
		dvector[index] = GetRateMultiplier(i);
		index++;
	}

	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	index++;
	
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}
	
	for(int p=0; p<nrrpart; ++p)
	{
		if(!fixrr[p])
		{
			for(i=0; i<nrr ; ++i) {
				dvector[index] = rr[p][i];
				index++;
			}
		}
	}

	for(i=0; i<nstatpart; ++i) {
		if(!fixstat[i])
		{
			for(j=0; j<ndim; ++j) {
				dvector[index] = profile[i][j];
				index++;
			}
		}
	}
	for (int i=0; i<ndim; i++)	{
		dvector[index] = dirweight[i];
		index++;
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	delete[] dvector;
}


void PartitionedRASGTRGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void PartitionedRASGTRGammaPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int i,j,nrr,nalpha,nbranch = GetNbranch(),nd,ndim,nrrpart,nstatpart;
	nalpha = PartitionedDGamRateProcess::GetNpart();
	nrr = GetNrr();
	nrrpart = PartitionedGTRProfileProcess::GetNpart();
	ndim = GetDim();
	nstatpart = PartitionedProfileProcess::GetNpart();
	nd = 2*nalpha + 2 + nbranch + nfreerr*nrr + nfreestat*ndim + ndim;

	double* dvector = new double[nd];

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	for(i=0; i<nalpha; ++i)
	{
		SetAlpha(i, dvector[index]);
		index++;
		ratemult[i] = dvector[index];
		index++;
	}

	branchalpha = dvector[index];
	index++;
	branchbeta = dvector[index];
	index++;

	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}

	for(int p=0; p<nrrpart; ++p)
	{
		if(!fixrr[p])
		{
			for(i=0; i<nrr; ++i) {
				rr[p][i] = dvector[index];
				index++;
			}
		}
	}
	for(i=0; i<nstatpart; ++i) {
		if(!fixstat[i])
		{
			for(j=0; j<ndim; ++j) {
				profile[i][j] = dvector[index];
				index++;
			}
		}
	}
	for (int i=0; i<ndim; i++)	{
		dirweight[i] = dvector[index];
		index++;
	}
	delete[] dvector;

	UpdateMatrices();
}


void PartitionedRASGTRGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	int cv = 0;
	int sitelogl = 0;
	int rr = 0;
	int rates = 0;
	int map = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;

	double tuning = 1;

	cvschemefile = "None";

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
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-p")	{
				i++;
				cvschemefile = argv[i];
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
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if (s == "-rr")	{
				rr = 1;
			}
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
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

	if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (rr)	{
		ReadRelRates(name,burnin,every,until);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void PartitionedRASGTRGammaPhyloProcess::ReadRelRates(string name, int burnin, int every, int until)	{

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

	double* meanrr = new double[GetNrr()];
	for (int k=0; k<GetNrr(); k++)	{
		meanrr[k] = 0;
	}

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int k=0; k<GetNrr(); k++)
		{
			double meank = 0.0;
			for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)
			{
				meank += rr[p][k] * PartitionedGTRProfileProcess::GetPartNsite(p);
			}

			meanrr[k] += meank / GetNsite();
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	for (int k=0; k<GetNrr(); k++)	{
		meanrr[k] /= samplesize;
	}
	ofstream os((name + ".meanrr").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetDim(); i++)	{
		for (int j=i+1; j<GetDim(); j++)	{
			os << GetStateSpace()->GetState(i) << '\t' << GetStateSpace()->GetState(j) << '\t' << meanrr[rrindex(i,j,GetDim())] << '\n';
		}
	}
	cerr << "mean relative exchangeabilities in " << name << ".meanrr\n";

}


void PartitionedRASGTRGammaPhyloProcess::GlobalSetTestData()	{
	testnsite = testdata->GetNsite();
	int* tmp = new int[testnsite * GetNtaxa()];
	testdata->GetDataVector(tmp);

	MESSAGE signal = SETTESTDATA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);

	if (cvschemefile == "None")	{
		cerr << "error : partition scheme file must be specified for cv test dataset\n";
		MPI_Finalize();
		exit(1);
	}

	PartitionScheme datascheme,testscheme;

	vector<PartitionScheme> cvschemes = PartitionedDGamRateProcess::ReadSchemes(cvschemefile, testdata->GetNsite(), myid, linkgam, unlinkgtr, rrtype);

	if(!linkgam)
	{
		datascheme = PartitionedDGamRateProcess::scheme;
		testscheme = cvschemes[2];
	}
	else
	{
		datascheme = PartitionedProfileProcess::scheme;
		testscheme = cvschemes[1];
	}

	if(datascheme.Npart != testscheme.Npart)
	{
		cerr << "error : incorrect number of partitions in cv test scheme\n";
		MPI_Finalize();
		exit(1);
	}

	MPI_Bcast(&(datascheme.Npart),1,MPI_INT,0,MPI_COMM_WORLD);
	for(int p = 0; p < datascheme.Npart; p++)
	{
		int ndatasites = datascheme.partSites[p].size();
		MPI_Bcast(&ndatasites,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&(datascheme.partSites[p][0]),ndatasites,MPI_INT,0,MPI_COMM_WORLD);

		int ntestsites = testscheme.partSites[p].size();
		MPI_Bcast(&ntestsites,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&(testscheme.partSites[p][0]),ntestsites,MPI_INT,0,MPI_COMM_WORLD);
	}

	delete[] tmp;
}

void PartitionedRASGTRGammaPhyloProcess::SlaveSetTestData()	{

	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	int* tmp = new int[testnsite * GetNtaxa()];
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);

	PartitionScheme testscheme,datascheme;
	MPI_Bcast(&(testscheme.Npart),1,MPI_INT,0,MPI_COMM_WORLD);
	datascheme.Npart = testscheme.Npart;

	datascheme.partSites.resize(datascheme.Npart);
	testscheme.partSites.resize(testscheme.Npart);

	datascheme.sitePart.resize(GetNsite());
	testscheme.sitePart.resize(testnsite);

	SetTestSiteMinAndMax();
	for(int p = 0; p < datascheme.Npart; p++)
	{
		int n;

		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		datascheme.partSites[p].resize(n);
		MPI_Bcast(&(datascheme.partSites[p][0]),n,MPI_INT,0,MPI_COMM_WORLD);

		for(int i = 0; i < datascheme.partSites[p].size(); i++)
		{
			int site = datascheme.partSites[p][i];
			datascheme.sitePart[site] = p;
		}

		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		testscheme.partSites[p].resize(n);
		MPI_Bcast(&(testscheme.partSites[p][0]),n,MPI_INT,0,MPI_COMM_WORLD);

		for(int i = 0; i < testscheme.partSites[p].size(); i++)
		{
			int site = testscheme.partSites[p][i];
			testscheme.sitePart[site] = p;
		}
	}

	int testwidth = testnsite/(nprocs-1);
	int width = GetNsite()/(nprocs-1);

	vector<int> partcounts(datascheme.Npart, 0);

	int min = 0;
	for(int proc = 0; proc < nprocs - 1; proc++)
	{
		if(min == sitemin)
			break;

		for(int site = min; site < min + width; site++)
		{
			int datasitepart = datascheme.sitePart[site];

			partcounts[datasitepart]++;
		}

		min += width;
	}

	vector<int> testcounts(datascheme.Npart, 0);

	for(int site = sitemin; site < sitemax; site++)
	{
		int datasitepart = datascheme.sitePart[site];

		testcounts[datasitepart]++;
	}

	for(int p = 0; p < partcounts.size(); p++)
	{
		int testmin = (partcounts[p] * testwidth) / width;
		int testmax = testmin + (testcounts[p] * testwidth) / width;

		int i = 0;
		for(int site = sitemin; site < sitemax; site++)
		{
			int datasitepart = datascheme.sitePart[site];

			if(datasitepart == p)
			{
				if(i < testmax - testmin)
				{
					int testsite = testscheme.partSites[p][testmin + i];

					data->SetTestData(testnsite,site,testsite,testsite+1,tmp);
				}
				else
				{
					sitemask[site-sitemin] = true;
				}

				i++;
			}
		}
	}

	delete[] tmp;
}

void PartitionedRASGTRGammaPhyloProcess::SlaveComputeCVScore()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	//sitemax = sitemin + testsitemax - testsitemin;
	
	UpdateConditionalLikelihoods();

	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
			total += sitelogL[i];
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	//sitemax = bksitemax;

}

void PartitionedRASGTRGammaPhyloProcess::SlaveComputeSiteLogL()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}


	UpdateConditionalLikelihoods();

	double* meansitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meansitelogl[i] = 0;
	}

	for (int i=sitemin; i<sitemax; i++)	{
		meansitelogl[i] = sitelogL[i];
	}

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	delete[] meansitelogl;

}
