
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

void PartitionedRASCATGTRFiniteGammaPhyloProcess::Create(Tree* intree, SequenceAlignment* indata, int nratecat,bool inlinkmult, int ncat,PartitionScheme rrscheme, PartitionScheme dgamscheme, int infixncomp, int inempmix, string inmixtype, int insitemin,int insitemax)	{
	PartitionedExpoConjugateGTRGammaPhyloProcess::Create(intree,indata,indata->GetNstate(),rrscheme,nratecat,inlinkmult,dgamscheme,insitemin,insitemax);
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

	UpdateMatrices();
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
	int rates = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;

	int rr = 0;
	int siteprofile = 0;
	
	bool ss = false;

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
			else if (s == "-var") {
                ppred = 4;
                i++;
                testdatafile = argv[i];
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
			else if (s == "-siteprofile")	{
				siteprofile = 1;
			}
			else if (s == "-ss")	{
				ss = true;
			}
			else if (s == "-rr")	{
				rr = 1;
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
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

	if (until == -1 && !ss)	{
		until = GetSize();
	}
	if (burnin == -1 && !ss)	{
		burnin = GetSize() / 5;
	}

	if ((GetNprocs() == 1) && (ppred || cv || sitelogl))	{
		cerr << "error : should run readpb_mpi in mpi mode, with at least 2 processes\n";
		MPI_Finalize();
		exit(1);
	}

	if (cv)	{
		cerr << "error : cannot use cross-validation with partitioned model\n";
		exit(1);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (siteprofile)	{
		ReadSiteProfiles(name,burnin,every,until);
	}
	else if (rr)	{
		ReadRelRates(name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior,testdatafile);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else if (ss)	{
		ReadSteppingStone(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
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

void PartitionedRASCATGTRFiniteGammaPhyloProcess::ReadRelRates(string name, int burnin, int every, int until)	{

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

void PartitionedRASCATGTRFiniteGammaPhyloProcess::ReadSiteProfiles(string name, int burnin, int every, int until)	{

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
		os << i+1;
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] /= samplesize;
			os << '\t' << sitestat[i][k];
		}
		os << '\n';
	}
	cerr << "mean site-specific profiles in " << name << ".siteprofiles\n";
	cerr << '\n';
}
