
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
#include "PartitionedRASCATGTRFiniteGammaPhyloProcess.h"
#include "Parallel.h"
#include "StringStreamUtils.h"
#include <string>
#include <list>
#include "TexTab.h"

void PartitionedRASCATGTRFiniteGammaPhyloProcess::Create(Tree* intree, SequenceAlignment* indata, int nratecat,int ncat,PartitionScheme rrscheme, PartitionScheme dgamscheme, int infixncomp, int inempmix, string inmixtype, int insitemin,int insitemax)	{
	PartitionedExpoConjugateGTRGammaPhyloProcess::Create(intree,indata,indata->GetNstate(),rrscheme,nratecat,dgamscheme,insitemin,insitemax);
	PartitionedExpoConjugateGTRFiniteProfileProcess::Create(indata->GetNstate(), rrscheme, ncat, infixncomp, inempmix, inmixtype);
	GammaBranchProcess::Create(intree);

	if (! partoccupancy)	{
		partoccupancy = new int**[PartitionedDGamRateProcess::GetNpart()];
		for (int d=0; d<PartitionedDGamRateProcess::GetNpart(); d++)	{
			partoccupancy[d] = new int*[PartitionedGTRProfileProcess::GetNpart()];
			for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
				partoccupancy[d][p] = new int[GetNmodeMax()];
			}
		}
	}
}

void PartitionedRASCATGTRFiniteGammaPhyloProcess::Delete() {
	if (partoccupancy)	{
		for (int d=0; d<PartitionedDGamRateProcess::GetNpart(); d++)	{
			for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)	{
				delete [] partoccupancy[d][p];
			}
			delete [] partoccupancy[d];
		}
		delete [] partoccupancy;
		partoccupancy = 0;
	}

	GammaBranchProcess::Delete();
	PartitionedExpoConjugateGTRFiniteProfileProcess::Delete();
	PartitionedExpoConjugateGTRGammaPhyloProcess::Delete();
}

void PartitionedRASCATGTRFiniteGammaPhyloProcess::UpdatePartOccupancyNumbers()
{
	for(int d = 0; d < PartitionedDGamRateProcess::GetNpart(); d++)
	{
		for (int p=0; p<PartitionedGTRProfileProcess::GetNpart(); p++)
		{
			for (int i=0; i<GetNcomponent(); i++)
				partoccupancy[d][p][i] = 0;
		}

		vector<int> partsites = PartitionedDGamRateProcess::GetPartSites(d);

		for (int i=0; i < partsites.size(); i++)	{
			partoccupancy[d][PartitionedGTRProfileProcess::GetSitePart(partsites[i])][MixtureProfileProcess::alloc[partsites[i]]]++;
		}
	}
}

// Importantly, this assumes that DGam partitions are always sub-partitions of GTR partitions
double PartitionedRASCATGTRFiniteGammaPhyloProcess::GetNormalizationFactor()
{
	if(occupancyNeedsUpdating)
	{
		UpdatePartOccupancyNumbers();

		double total = 0;
		for (int dgampart=0; dgampart<PartitionedDGamRateProcess::GetNpart(); dgampart++)
		{
			vector<int> partsites = PartitionedDGamRateProcess::GetPartSites(dgampart);

			size_t rrpart = PartitionedGTRProfileProcess::GetSitePart(partsites.front());

			total += GetNormPartRate(dgampart, rrpart) * PartitionedDGamRateProcess::GetRateMultiplier(dgampart) * partsites.size();
		}
		normFactor = total / GetNsite();

		occupancyNeedsUpdating = false;
	}

	return normFactor;
}

double PartitionedRASCATGTRFiniteGammaPhyloProcess::GetNormPartRate(int d, int p)
{
		double norm = 0;
		int tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (partoccupancy[d][p][k])	{
				norm += (partoccupancy[d][p][k] + 1) * GetNormRate(p, k);
				tot += partoccupancy[d][p][k] + 1;
			}
		}

		norm /= tot;
		return norm;
	}

void PartitionedRASCATGTRFiniteGammaPhyloProcess::GlobalUpdateParameters()	{
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

	int i,j,nrr,nalpha,nbranch = GetNbranch(),nd,ndim,nrrpart,ncomponent,ni;
	nalpha = PartitionedDGamRateProcess::GetNpart();
	nrr = GetNrr();
	nrrpart = PartitionedGTRProfileProcess::GetNpart();
	ndim = GetDim();
	ncomponent = GetNmodeMax();
	nd = 2*nalpha + 2 + nbranch + nfreerr*nrr + ncomponent*(ndim + 1) + ndim;
	ni = 1 + GetNsite();

	int* ivector = new int[ni];
	double* dvector = new double[nd];

	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// GlobalBroadcastTree();

	// First we assemble the vector of doubles for distribution
	int index = 0;
	for(i=0; i<nalpha; ++i)
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
			for(i=0; i<nrr; ++i) {
				dvector[index] = rr[p][i];
				index++;
			}
		}
	}
	for(i=0; i<ncomponent; ++i) {
		for(j=0; j<ndim; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
		dvector[index] = weight[i];
		index++;
	}
	for (int i=0; i<ndim; i++)	{
		dvector[index] = dirweight[i];
		index++;
	}


	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<GetNsite(); ++i) {
		ivector[1+i] = FiniteProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	delete [] ivector;
	delete [] dvector;
}


void PartitionedRASCATGTRFiniteGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	case REALLOC_MOVE:
		SlaveIncrementalFiniteMove();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void PartitionedRASCATGTRFiniteGammaPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int i,j,nrr,nalpha,nbranch = GetNbranch(),nd,ndim,nrrpart,ncomponent,ni;
	nalpha = PartitionedDGamRateProcess::GetNpart();
	nrr = GetNrr();
	nrrpart = PartitionedGTRProfileProcess::GetNpart();
	ndim = GetDim();
	ncomponent = GetNmodeMax();
	nd = 2*nalpha + 2 + nbranch + nfreerr*nrr + ncomponent*(ndim + 1) + ndim;
	ni = 1 + GetNsite();

	int* ivector = new int[ni];
	double* dvector = new double[nd];

	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	for(i = 0; i < nalpha; i++)
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
			for(i=0; i<nrr ; ++i) {
				rr[p][i] = dvector[index];
				index++;
			}
		}
	}

	for(i=0; i<ncomponent; ++i) {
		for(j=0; j<ndim; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
		weight[i] = dvector[index];
		index++;
	}

	for (int i=0; i<ndim; i++)	{
		dirweight[i] = dvector[index];
		index++;
	}

	Ncomponent = ivector[0];
	for(i=0; i<GetNsite(); ++i) {
		FiniteProfileProcess::alloc[i] = ivector[1+i];
	}
	delete[] dvector;
	delete[] ivector;
	// some upate here ?
}


void PartitionedRASCATGTRFiniteGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

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
	int map = 0;
	string testdatafile = "";
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
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if (s == "-p")	{
				i++;
				cvschemefile = argv[i];
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
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
}


void PartitionedRASCATGTRFiniteGammaPhyloProcess::GlobalSetTestData()	{
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
		datascheme = PartitionedGTRProfileProcess::scheme;
		testscheme = cvschemes[0];
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

void PartitionedRASCATGTRFiniteGammaPhyloProcess::SlaveSetTestData()	{

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
		int testmin = partcounts[p] * testwidth / width;
		int testmax = testmin + testcounts[p] * testwidth / width;

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

void PartitionedRASCATGTRFiniteGammaPhyloProcess::SlaveComputeCVScore()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	//sitemax = sitemin + testsitemax - testsitemin;
	double** sitelogl = new double*[GetNsite()];
	for(int i = sitemin; i < sitemax; i++)
	{
		if(!sitemask[i-sitemin])
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for(int i = sitemin; i < sitemax; i++)	{
			PartitionedExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
		}
		UpdateConditionalLikelihoods();
		for(int i = sitemin; i < sitemax; i++)	{
			if(!sitemask[i-sitemin])
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double total = 0;
	for(int i = sitemin; i < sitemax; i++)	{
		if(!sitemask[i-sitemin])
		{
			double max = 0;
			for (int k=0; k<GetNcomponent(); k++)	{
				if ((!k) || (max < sitelogl[i][k]))	{
					max = sitelogl[i][k];
				}
			}
			double tot = 0;
			double totweight = 0;
			for (int k=0; k<GetNcomponent(); k++)	{
				tot += weight[k] * exp(sitelogl[i][k] - max);
				totweight += weight[k];
			}
			total += log(tot) + max;
		}
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for(int i = sitemin; i < sitemax; i++){
		if(!sitemask[i-sitemin])
		delete[] sitelogl[i];
	}
	delete[] sitelogl;

	//sitemax = bksitemax;

}

void PartitionedRASCATGTRFiniteGammaPhyloProcess::SlaveComputeSiteLogL()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			PartitionedExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double* meansitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meansitelogl[i] = 0;
	}
	for (int i=sitemin; i<sitemax; i++)	{
		double max = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if ((!k) || (max < sitelogl[i][k]))	{
				max = sitelogl[i][k];
			}
		}
		double tot = 0;
		double totweight = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			tot += weight[k] * exp(sitelogl[i][k] - max);
			totweight += weight[k];
		}
		meansitelogl[i] = log(tot) + max;
	}

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;
	delete[] meansitelogl;
}
