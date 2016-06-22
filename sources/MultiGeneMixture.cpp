
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
#include "Parallel.h"
#include <string.h>

#include "MultiGeneMixture.h"

void MultiGeneMixture::MakeFiles()	{
	for (int gene=0; gene<Ngene; gene++)	{
		ofstream os((name + genename[gene] + ".treelist").c_str());
	}
}

void MultiGeneMixture::AllocateAlignments(string datafile, string treefile, int dc)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
	ifstream* tis = 0;
	if (treefile != "None")	{
		tis = new ifstream(treefile.c_str());
		int n;
		(*tis) >> n;
		if (n != Ngene)	{
			cerr << "error: non matching number of alignments and trees\n";
			exit(1);
		}
	}
	genename = new string[Ngene];
	treename = new string[Ngene];
	genesize = new int[Ngene];
	genealloc = new int[Ngene];
	int* geneweight = new int[Ngene];
	// genedata = new SequenceAlignment*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		if (tis)	{
			(*tis) >> treename[gene];
		}
		else	{
			treename[gene] = "None";
		}
		SequenceAlignment* data = new FileSequenceAlignment(genename[gene],0,myid);
		if (dc)	{
			data->DeleteConstantSites();
		}
		int nstate = data->GetNstate();
		if (! gene)	{
			Nstate = nstate;
			statespace = data->GetStateSpace();
		}
		else	{
			if (Nstate != nstate)	{
				cerr << "error: all data files do not have the same alphabet\n";
				exit(1);
			}
		}

		genesize[gene] = data->GetNsite();
		// geneweight[gene] = data->GetNsite();
		geneweight[gene] = data->GetNsite() * data->GetNtaxa();
		delete data;
	}
	delete tis;
	tis = 0;

	// sort alignments by decreasing size
	int permut[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		permut[gene] = gene;
	}
	for (int i=0; i<Ngene; i++)	{
		for (int j=Ngene-1; j>i; j--)	{
			if (geneweight[permut[i]] < geneweight[permut[j]])	{
			// if (genesize[permut[i]] < genesize[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}

	int totsize[nprocs];
	for (int i=0; i<nprocs; i++)	{
		totsize[i] = 0;
	}

	for (int i=0; i<Ngene; i++)	{
		int gene = permut[i];
		int size = geneweight[gene];
		// int size = genesize[gene];

		int min = 0;
		int jmin = 0;
		for (int j=1; j<nprocs; j++)	{
			if ((j==1) || (min > totsize[j]))	{
				min = totsize[j];
				jmin = j;
			}
		}
		genealloc[gene] = jmin;
		totsize[jmin] += size;
	}

	if (totsize[0])	{
		cerr << "error in alloc\n";
		exit(1);
	}
	int total = 0;
	for (int i=1; i<nprocs; i++)	{
		int tot = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				tot += geneweight[gene];
				// tot += genesize[gene];
				total++;
			}
		}
		if (tot != totsize[i])	{
			cerr << "error in allocation\n";
			cerr << tot << '\t' << totsize[i] << '\n';
			exit(1);
		}
	}
	if (total != Ngene)	{
		cerr << "error in total allocation\n";
		exit(1);
	}

	genelnL = new double[Ngene];
	tmpgenelnL = new double[Ngene];

	globalnsite = new int[nprocs];
	for (int i=0; i<nprocs; i++)	{
		globalnsite[i] = 0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if ((genealloc[gene] < 0) || (genealloc[gene] >= nprocs))	{
			cerr << "alloc : " << genealloc[gene] << '\t' << gene << '\n';
			exit(1);
		}
		globalnsite[0] += genesize[gene];
		globalnsite[genealloc[gene]] += genesize[gene];
	}
	if (! myid)	{
		cerr << '\n';
		for (int i=1; i<nprocs; i++)	{
			cerr << i << '\t' << globalnsite[i] << '\n';
		}
		cerr << '\n';
		cerr << "total: " << GetGlobalNsite() << '\n';
	}
	
	// check total size
	if (! myid)	{
		int tot = 0;
		for (int i=1; i<nprocs; i++)	{
			tot += globalnsite[i];
		}
		if (tot != globalnsite[0])	{
			cerr << "error in total size during gene allocation\n";
			exit(1);
		}
		int tot2 = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			tot2 += genesize[gene];
		}
		if (tot2 != tot)	{
			cerr << "error during alloc: total size does not match\n";
			exit(1);
		}
	}

	genealpha = new double[Ngene];
	genelength = new double[Ngene];
	delete[] geneweight;
}

double MultiGeneMixture::GetLogLikelihood()	{
	GlobalCollectGeneLikelihoods();
	lnL = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		lnL += genelnL[gene];
	}
	return lnL;
}

void MultiGeneMixture::GlobalCollectGeneLikelihoods()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = LIKELIHOOD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genelnL[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGeneMixture::SlaveSendGeneLikelihoods()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetLogLikelihood();
		}
	}
	MPI_Send(genelnL,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneMixture::GetMeanLength()	{
	GlobalCollectGeneLengths();
	double total = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		total += genelength[gene];
	}
	return total / Ngene;
}

void MultiGeneMixture::GlobalCollectGeneLengths()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = LENGTH;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genelength[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGeneMixture::SlaveSendGeneLengths()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelength[gene] = process[gene]->GetTotalLength();
		}
	}
	MPI_Send(genelength,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneMixture::GetMeanAlpha()	{
	GlobalCollectGeneAlphas();
	double total = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		total += genealpha[gene];
	}
	return total / Ngene;
}

void MultiGeneMixture::GlobalCollectGeneAlphas()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = ALPHA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genealpha[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGeneMixture::SlaveSendGeneAlphas()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genealpha[gene] = process[gene]->GetAlpha();
		}
	}
	MPI_Send(genealpha,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalGeneMove()	{
	assert(myid == 0);
	MESSAGE signal = GENE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalSample()	{
	assert(myid == 0);
	MESSAGE signal = SAMPLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalUnfold()	{
	assert(myid == 0);
	MESSAGE signal = UNFOLD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalCollapse()	{
	assert(myid == 0);
	MESSAGE signal = COLLAPSE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::SaveTrees()	{
	assert(myid == 0);
	MESSAGE signal = SAVETREES;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalFromStream(istream& is)	{
	assert(myid == 0);
	MESSAGE signal = FROMSTREAM;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for (int j=1; j<nprocs; j++)	{
		unsigned int len = 0;
		is >> len;
		char* bvector = new char[len];
		// unsigned char* bvector = new unsigned char[len];
		for (unsigned int i=0; i<len; i++)	{
			// is >> bvector[i];
			is.get(bvector[i]);
		}
		MPI_Send(&len,1,MPI_INT,j,TAG1,MPI_COMM_WORLD);
		MPI_Send(bvector,len,MPI_UNSIGNED_CHAR,j,TAG1,MPI_COMM_WORLD);
		delete[] bvector;
	}
}

void MultiGeneMixture::SlaveFromStream()	{

	MPI_Status stat;
	unsigned int len = 0;
	MPI_Recv(&len, 1, MPI_INT,0,TAG1,MPI_COMM_WORLD,&stat);
	char* bvector = new char[len];
	MPI_Recv(bvector,len,MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD,&stat);
	ostringstream os;
	for (unsigned int i=0; i<len; i++)	{
		os << bvector[i];
	}
	delete[] bvector;
	string s = os.str();
	istringstream is(s);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->FromStream(is);
		}
	}
}

void MultiGeneMixture::GlobalToStream(ostream& os)	{

	assert(myid == 0);
	MESSAGE signal = TOSTREAM;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		unsigned int len = 0;
		MPI_Recv(&len, 1, MPI_INT,j,TAG1,MPI_COMM_WORLD,&stat);
		os << len << '\n';
		char* bvector = new char[len];
		MPI_Recv(bvector,len,MPI_UNSIGNED_CHAR,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int i=0; i<len; i++)	{
			os << bvector[i];
		}
		delete[] bvector;
	}
}

void MultiGeneMixture::SlaveToStream()	{
	ostringstream os;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ToStream(os);
		}
	}
	string s = os.str();
	unsigned int len = s.length();
	char* bvector = new char[len];
	for (unsigned int i=0; i<len; i++)	{
		bvector[i] = s[i];
	}
	MPI_Send(&len,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(bvector,len,MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;
}

void MultiGeneMixture::SlaveGeneMove()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Move();
		}
	}
}

void MultiGeneMixture::SlaveSample()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Sample();
		}
	}
}

void MultiGeneMixture::SlaveUnfold()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Unfold();
		}
	}
}

void MultiGeneMixture::SlaveCollapse()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Collapse();
		}
	}
}

void MultiGeneMixture::SlaveSaveTrees()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SaveTree(name + genename[gene]);
		}
	}
}

void MultiGeneMixture::WaitLoop()	{
	MESSAGE signal;
	do {
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		if (signal == KILL) break;
		SlaveExecute(signal);
	} while(true);
}

void MultiGeneMixture::SlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case FROMSTREAM:
		SlaveFromStream();
		break;
	case TOSTREAM:
		SlaveToStream();
		break;
	case SAVETREES:
		SlaveSaveTrees();
		break;
	case LIKELIHOOD:
		SlaveSendGeneLikelihoods();
		break;
	case LENGTH:
		SlaveSendGeneLengths();
		break;
	case ALPHA:
		SlaveSendGeneAlphas();
		break;
	case SAMPLE:
		SlaveSample();
		break;
	case UNFOLD:
		SlaveUnfold();
		break;
	case COLLAPSE:
		SlaveCollapse();
		break;
	case UPDATE_SPROFILE:
		SlaveUpdateSiteProfileSuffStat();
		break;
	case PARAMETER_DIFFUSION:
		SlaveUpdateParameters();
		break;
	case GENE_MOVE:
		SlaveGeneMove();
		break;
	case MIX_MOVE:
		SlaveMixMove();
		break;
	
	default:
		MultiGeneMixture::SlaveExecute(signal);
	}
}

