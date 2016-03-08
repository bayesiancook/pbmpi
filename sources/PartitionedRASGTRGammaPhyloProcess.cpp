
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

	int i,j,nrr,nalpha,nbranch = GetNbranch(),ni,nd,ndim,nrrpart,nstatpart;
	nalpha = PartitionedDGamRateProcess::GetNpart();
	nrr = GetNrr();
	nrrpart = PartitionedGTRProfileProcess::GetNpart();
	ndim = GetDim();
	nstatpart = PartitionedProfileProcess::GetNpart();
	nd = 2 + 2*nalpha + nbranch + nrrpart*nrr + nstatpart*ndim + ndim + 1;
	ni = 1 + GetNsite();
	int ivector[ni];
	double dvector[nd]; 
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
	
	for(int p=0; p<nrrpart; ++i)
	{
		for(i=0; i<nrr ; ++i) {
			dvector[index] = rr[p][i];
			index++;
		}
	}

	for(i=0; i<nstatpart; ++i) {
		for(j=0; j<ndim; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
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

	int i,j,nrr,nalpha,nbranch = GetNbranch(),ni,nd,ndim,nrrpart,nstatpart;
	nalpha = PartitionedDGamRateProcess::GetNpart();
	nrr = GetNrr();
	nrrpart = PartitionedGTRProfileProcess::GetNpart();
	ndim = GetDim();
	nstatpart = PartitionedProfileProcess::GetNpart();
	nd = 2 + 2*nalpha + nbranch + nrrpart*nrr + nstatpart*ndim + ndim + 1;
	ni = 1 + GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int index = 0;
	for(i=0; i<nalpha; ++i)
	{
		SetAlpha(i, dvector[index]);
		index++;
		ratemult[i] = dvector[index];
		index++;
	}
	index++;
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
		for(i=0; i<nrr; ++i) {
			rr[p][i] = dvector[index];
			index++;
		}
	}
	for(i=0; i<nstatpart; ++i) {
		for(j=0; j<ndim; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}
	delete[] dvector;
	delete[] ivector;
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
	int map = 0;
	string testdatafile = "";

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
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
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

void PartitionedRASGTRGammaPhyloProcess::SlaveComputeCVScore()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	sitemax = sitemin + testsitemax - testsitemin;
	
	UpdateConditionalLikelihoods();

	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
		total += sitelogL[i];
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	sitemax = bksitemax;

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

vector<PartitionScheme> PartitionedRASGTRGammaPhyloProcess::ReadSchemes(string schemefile)
{
	string error = "Error: improperly formatted scheme file\n";

	ifstream theStream((Path + schemefile).c_str());

	map<string, size_t> fixrrparts;
	map<string, size_t> fixstatparts;

	PartitionScheme rrscheme(GetNsite());
	PartitionScheme statscheme(GetNsite());
	PartitionScheme dgamscheme(GetNsite());

	vector<size_t> partsites;

	// get the unlinked partition info
	string line;
	while(getline(theStream, line))
	{
		stringstream ss(line);

		string item;
		vector<string> elems;
		while (getline(ss, item, ',')) {
			elems.push_back(item);
		}

		if(elems.size() < 2)
		{
			cerr << error;
			exit(1);
		}

		string type;
		bool fixprof;

		for(vector<string>::iterator it = elems.begin(); it != elems.end(); it++)
		{
			if(it == elems.begin())
			{
				ss.clear();
				ss.str(*it);

				ss >> type;

				if(type.empty())
				{
						cerr << error;
						exit(1);
				}

				std::transform(type.begin(), type.end(), type.begin(), ::tolower);

				if(type == "gtr")
				{
					type = "None";
					fixprof = false;
				}
				else
				{
					fixprof = ((*type.rbegin()) != 'f' || type == "dayhoff" );

					if(!fixprof)
					{
						// remove 'f' character
						if (type.size () > 0)  type.resize (type.size () - 1);
					}
					else if(fixstatparts.find(type) == fixstatparts.end())
					{
						statscheme.partType.push_back(type);
						fixstatparts[type] = statscheme.Npart++;
					}
				}

				if(fixrrparts.find(type) == fixrrparts.end())
				{
					rrscheme.partType.push_back(type);
					fixrrparts[type] = rrscheme.Npart++;
				}
			}
			else
			{
				ss.clear();
				ss.str(*it);

				while(ss >> item) {}

				ss.clear();
				ss.str(item);

				vector<int> range;
				while (getline(ss, item, '-'))
				{
					if(IsInt(item))
					{
						range.push_back(atoi(item.c_str()) - 1);
					}
					else
					{
						cerr << error;
						exit(1);
					}
				}

				if(range.empty() || range.size() > 2)
				{
					cerr << error;
					exit(1);
				}

				if(range.size() == 1)
				{
					partsites.push_back(range[0]);

					dgamscheme.sitePart[range[0]] = dgamscheme.Npart;

					rrscheme.sitePart[range[0]] = fixrrparts[type];

					if(fixprof)
					{
						statscheme.sitePart[range[0]] = fixstatparts[type];
					}
					else
					{
						statscheme.sitePart[range[0]] = statscheme.Npart;
					}
				}
				else
				{
					for(size_t i = range[0]; i <= range[1]; i++)
					{
						partsites.push_back(i);

						dgamscheme.sitePart[i] = dgamscheme.Npart;

						rrscheme.sitePart[i] = fixrrparts[type];

						if(fixprof)
						{
							statscheme.sitePart[i] = fixstatparts[type];
						}
						else
						{
							statscheme.sitePart[i] = statscheme.Npart;
						}
					}
				}
			}
		}

		dgamscheme.Npart++;

		if(!fixprof)
		{
			statscheme.partType.push_back("None");
			statscheme.Npart++;
		}

	}

	if(partsites.size() != GetNsite())
	{
		size_t rrpart;
		if(fixrrparts.find("None") != fixrrparts.end())
		{
			rrpart = fixrrparts["None"];
		}
		else
		{
			rrscheme.partType.push_back("None");
			rrpart = rrscheme.Npart++;
		}

		statscheme.partType.push_back("None");

		std::sort(partsites.begin(), partsites.end());

		vector<size_t>::iterator it = partsites.begin();
		for(size_t site = 0; site < GetNsite(); site++)
		{
			if(*it != site)
			{
				rrscheme.sitePart[site] = rrpart;
				statscheme.sitePart[site] = statscheme.Npart;
				dgamscheme.sitePart[site] = dgamscheme.Npart;
			}
			else
			{
				it++;
			}
		}

		statscheme.Npart++;
		dgamscheme.Npart++;
	}

	rrscheme.update();
	statscheme.update();
	dgamscheme.update();

	vector<PartitionScheme> schemes;

	schemes.push_back(rrscheme);
	schemes.push_back(statscheme);
	schemes.push_back(dgamscheme);

	return schemes;
}
