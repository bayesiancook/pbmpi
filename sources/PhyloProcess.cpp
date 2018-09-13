
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

#include "Random.h"
#include "PhyloProcess.h"
#include <string>

#include "Parallel.h"
extern MPI_Datatype Propagate_arg;

#include "TexTab.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//	* New
//-------------------------------------------------------------------------

void PhyloProcess::New(int unfold)	{

	ReadData(datafile);
	SetProfileDim();
	CreateMPI(GetData()->GetNsite());
	if (reshuffle)	{
		if ((sis == 1) || (topobf == 1))	{
			GlobalReshuffleSites();
		}
	}
	
	SetTree(treefile);
	if (partitionfile != "None")	{
		SetPartition(partitionfile);
	}

	if (! GetMyid())	{
        if (roottax1 != "None")	{
            cerr << "reroot tree\n";
            Link* newroot = GetTree()->GetLCA(roottax1,roottax2);
            if (!newroot)	{
                cerr << "error when rerooting\n";
                exit(1);
            }
            GlobalRootAt(newroot);
        }
	}

	Create();

	if (! GetMyid())	{

		if (topobf == 1)	{
			bffrac = 0;
		}
		if (topobf == 2)	{
			bffrac = -bfnfrac;
		}
		/*
		if (topobf == 3)	{
			cerr << "fix topo bf\n";
			cerr << "bffrac : " << bffrac << '\n';
		}
		if (topobf == 4)	{
			cerr << "fix topo bl\n";
			cerr << "bffrac : " << bffrac << '\n';
		}
		*/
		if (sis == 1)	{
			sisfrac = 0;
			if (! siscutoff)	{
				sislevel = 1;
			}
			SetSIS();
		}
		if (sis == 2)	{
			/*
			cerr << "fix BF\n";
			cerr << "sisfrac : " << sisfrac << '\n';
			*/
			if (! siscutoff)	{
				sislevel = 1;
			}
			SetSIS();
		}
		if (sis == 1)	{
			PriorSample();
		}
		else	{
			Sample();
		}

		if (unfold)	{
			if (topobf)	{
				GlobalSetTopoBF();
			}
			if (sis)	{

				if (sumovercomponents > 0)	{
					GlobalActivateSumOverComponents();
				}
			}
			GlobalUnfold();
		}
	}

	if (BPP)	{
		BPP->RegisterWithTaxonSet(GetData()->GetTaxonSet());
	}
}

//-------------------------------------------------------------------------
//	* Open
//-------------------------------------------------------------------------

void PhyloProcess::Open(istream& is, int unfold)	{

	ReadData(datafile);
	SetProfileDim();
	CreateMPI(GetData()->GetNsite());

	SetTreeFromString(treestring);

	if (partitionfile != "None")	{
		SetPartition(partitionfile);
	}

	Create();


	if (unfold)	{

		if (! GetMyid()) {

			FromStream(is);

			if (sis)	{
				SetSIS();
			}

			if (topobf)	{
				GlobalSetTopoBF();
			}
			GlobalReadSiteRankFromStream(is);
			GlobalUnfold();
		}
	}
}

void PhyloProcess::MakeObservedArray()	{

	observedarray = new int*[data->GetNsite()];
	for (int i=0; i<data->GetNsite(); i++)	{
		observedarray[i] = new int[data->GetNstate()];
		for (int k=0; k<data->GetNstate(); k++)	{
			observedarray[i][k] = 0;
		}
	}
	for (int i=0; i<data->GetNsite(); i++)	{
		for (int j=0; j<GetNtaxa(); j++)	{
			if (data->GetState(j,i) != -1)	{
				observedarray[i][data->GetState(j,i)] = 1;
			}
		}
	}	
}

void PhyloProcess::Monitor(ostream& os)  {
	os << "matrix uni" << '\t' << SubMatrix::GetUniSubCount() << '\n';
	os << "inf prob  " << '\t' << GetInfProbCount() << '\n';
	os << "stat inf  " << '\t' << GetStatInfCount() << '\n';
	if (sprtry)	{
		os << "spr " << '\t' << (100 * spracc) / sprtry << '\n';
	}
	if (mhsprtry)	{
		os << "mhspr " << '\t' << (100 * mhspracc) / mhsprtry << '\n';
	}
	if (tsprtry)	{
		os << "tspr " << '\t' << (100 * tspracc) / tsprtry << '\n';
		os << "fraction of temperedmoves: " << (100 * tsprtmp) / tsprtot << '\n';
		if (tsprtmp)	{
			os << "10" << '\t' << (100 * tsprtmpacc10) / tsprtmp << '\n';
			os << "01" << '\t' << (100 * tsprtmpacc01) / tsprtmp << '\n';
			os << "11" << '\t' << (100 * tsprtmpacc11) / tsprtmp << '\n';
			os << "00" << '\t' << (100 * tsprtmpacc00) / tsprtmp << '\n';
		}
	}
	if (nnitry)	{
		os << "nni " << '\t' << (100 * nniacc) / nnitry << '\n';
	}
	if (bppsprtry)	{
		os << "bppspr " << '\t' << (100 * bppspracc) / bppsprtry << '\n';
	}
	if (tbppsprtry)	{
		os << "tbppspr " << '\t' << (100 * tbppspracc) / tbppsprtry << '\n';
	}
	if (proftry)	{
		os << "profile moves " << '\t' << (100 * profacc) / proftry << '\n';
	}
	if (rrtry)	{
		os << "rr moves " << '\t' << (100 * rracc) / rrtry << '\n';
	}
	if (ziptopotry)	{
		os << "zip topo : " << '\t' << (100 * ziptopoacc) / ziptopotry << '\n';
	}
	if (fasttopotry)	{
		os << "fast topo moves: \n";
		os << "topo changed : " << '\t' << (100 * fasttopochange) / fasttopotry << '\n';
		if (fasttopochange)	{
			os << "accepted     : " << '\t' << (100 * fasttopoacc) / fasttopochange << '\n';
		}
		if (anntot)	{
			os << "tempered fraction: " << '\t' << (100 * anntmp) / anntot << '\n';
			if (anntmp)	{
				os << "10" << '\t' << (100 * anntmpacc10) / anntmp << '\n';
				os << "01" << '\t' << (100 * anntmpacc01) / anntmp << '\n';
				os << "11" << '\t' << (100 * anntmpacc11) / anntmp << '\n';
				os << "00" << '\t' << (100 * anntmpacc00) / anntmp << '\n';
			}
		}
	}
	if (! fixtopo)	{
		double totaltime = nnichrono.GetTime() + sprchrono.GetTime() + tsprchrono.GetTime();
		os << "nni  time : " << nnichrono.GetTime() / totaltime << '\n';
		os << "spr  time : " << sprchrono.GetTime() / totaltime << '\n';
		os << "tspr time : " << tsprchrono.GetTime() / totaltime << '\n';
	}
}

void PhyloProcess::SetParameters(string indatafile, string intreefile, string inpartitionfile, int iniscodon, GeneticCodeType incodetype, int insis, double insisfrac, int insisnfrac, int insisnrep, double insiscutoff, int infixtopo, int infixroot, string inroottax1, string inroottax2, int intopoburnin, int intopobf, int inbfburnin, double inbffrac, int inbfnfrac, int inbfnrep, double inblfactor, string inblfile, int inNSPR, int inNMHSPR, int inNTSPR, int intemperedbl, int intemperedgene, int intemperedrate, double intopolambda, double intopomu, int intoponstep, int inNNNI, int innspec, int inntspec, string intaxon1, string intaxon2, string intaxon3, string intaxon4, int inbpp, int innbpp, int inntbpp, int inbppnstep, string inbppname, double inbppcutoff, double inbppbeta, int inprofilepriortype, int indc, int infixbl, int insumovercomponents, int inproposemode, int inallocmode, int infasttopo, double infasttopofracmin, int infasttoponstep, int infastcondrate, int indirpriortype, int innstatcomp, int inpriorempmix, string inpriormixtype, int infixstatweight, int infixstatalpha, int infixstatcenter, int inreshuffle)	{

	reshuffle = inreshuffle;

	datafile = indatafile;
	treefile = intreefile;
	partitionfile = inpartitionfile;
	iscodon = iniscodon;
	codetype = incodetype;

	sis = insis;
	sisfrac = insisfrac;
	sisnfrac = insisnfrac;
	sisnrep = insisnrep;
	siscutoff = insiscutoff;

	fixtopo = infixtopo;
	fixroot = infixroot;
	roottax1 = inroottax1;
	roottax2 = inroottax2;
	topoburnin = intopoburnin;

	topobf = intopobf;
	bfburnin = inbfburnin;
	bffrac = inbffrac;
	bfnfrac = inbfnfrac;
	bfnrep = inbfnrep;
	blfactor = inblfactor;
	blfile = inblfile;

	NSPR = inNSPR;
	NMHSPR = inNMHSPR;
	NTSPR = inNTSPR;
	temperedbl = intemperedbl;
	temperedgene = intemperedgene;
	temperedrate = intemperedrate;
	topolambda = intopolambda;
	topomu = intopomu;
	toponstep = intoponstep;
	NNNI = inNNNI;
	nspec = innspec;
	ntspec = inntspec;
	SetSpecialSPR(intaxon1,intaxon2,intaxon3,intaxon4);
	bpp = inbpp;
	nbpp = innbpp;
	ntbpp = inntbpp;
	bppnstep = inbppnstep;
	bppname = inbppname;
	bppcutoff = inbppcutoff;
	bppbeta = inbppbeta;
	
	BPP = 0;
	if (bpp == 1)	{
		cerr << "make new BPP\n";
		BPP = new UnrootedBPP(bppname,bppcutoff,bppbeta);
	}
	else if (bpp == 2)	{
		cerr << "make new CCP\n";
		BPP = new UnrootedCCP(bppname,bppcutoff,bppbeta);
	}
	else if (bpp == 3)	{
		cerr << "make new CCP\n";
		BPP = new UnrootedCCP(bppname,bppcutoff,bppbeta);
	}

	profilepriortype = inprofilepriortype;
	dc = indc;
	fixbl = infixbl;
	sumovercomponents = insumovercomponents;
	if (sis && (! sumovercomponents))	{
		sumovercomponents = 50;
	}
	proposemode = inproposemode;
	allocmode = inallocmode;
	// sumratealloc = insumratealloc;
	fasttopo = infasttopo;
	fasttopofracmin = infasttopofracmin;
	fasttoponstep = infasttoponstep;
	fastcondrate = infastcondrate;

	dirpriortype = indirpriortype;
	Nstatcomp = innstatcomp;
	priorempmix = inpriorempmix;
	priormixtype = inpriormixtype;
	fixstatweight = infixstatweight;
	fixstatalpha = infixstatalpha;
	fixstatcenter = infixstatcenter;
}

void PhyloProcess::SetBranchesToCollapse(string blfile)	{

	cerr << "in set branches to collapse\n";
	ResetBranchAlloc();
	if (blfile != "None")	{
		ifstream is(blfile.c_str());
		int n;
		is >> n;
		for (int i=0; i<n; i++)	{
			string tax1, tax2;
			is >> tax1 >> tax2;
			Link* down = GetTree()->GetLCA(tax1,tax2);
			if (! down)	{
				cerr << "error in set topobf: did not find MRCA of " << tax1 << " and " << tax2 << '\n';
				exit(1);
			}
			Link* up = GetTree()->GetAncestor(down);
			SetBranchAlloc(up->GetBranch()->GetIndex(),1);
			cerr << tax1 << '\t' << tax2 << '\t' << up->GetBranch()->GetIndex() << '\t' << branchalloc[up->GetBranch()->GetIndex()] << '\n';
		}
	}
	else	{
		Link* down = GetTree()->GetLCA(taxon1,taxon2);
		if (! down)	{
			cerr << "error in PhyloProcess::SetTopoBF: did not find LCA of " << taxon1 << " and " << taxon2 << '\n';
			exit(1);
		}
		Link* up = GetTree()->GetAncestor(down);
		SetBranchAlloc(up->GetBranch()->GetIndex(),1);
		cerr << taxon1 << '\t' << taxon2 << '\t' << up->GetBranch()->GetIndex() << '\t' << branchalloc[up->GetBranch()->GetIndex()] << '\n';
	}
}

void PhyloProcess::GlobalSetTopoBF()	{

	if (sumovercomponents > 0)	{
		GlobalActivateSumOverComponents();
	}

	GlobalSetBFFrac();
	if ((topobf == 2) || (topobf == 4))	{
		double scale = 0.1;
		if (bffrac < -1)	{
			scale *= exp((bfnfrac + bffrac)*log(blfactor));
		}
		else if (bffrac == 0)	{
			scale *= exp((bfnfrac - 1)*log(blfactor));
		}
		else	{
			scale *= exp((bfnfrac - bffrac -1)*log(blfactor));
		}
		SetBranchScaling(scale,1);	
		cerr << "set branch scaling to : " << scale << '\n';
		// cerr << bffrac << '\t' << bfnfrac << '\n';
	}
	GlobalBackupTree();
	Link* down = GetTree()->GetLCA(taxon1,taxon2);

	if (! down)	{
		cerr << "error in PhyloProcess::SetTopoBF: did not find LCA of " << taxon1 << " and " << taxon2 << '\n';
		exit(1);
	}
	Link* up = GetTree()->GetAncestor(down);

	if ((topobf == 2) || (topobf == 4))	{
		SetBranchesToCollapse(blfile);
	}
	Link* fromdown = GlobalDetach(down,up);
	Link* todown = GetTree()->GetLCA(taxon3,taxon4);
	if (! todown)	{
		cerr << "error in PhyloProcess::SetTopoBF: did not find LCA of " << taxon3 << " and " << taxon4 << '\n';
		exit(1);
	}
	Link* toup = GetTree()->GetAncestor(todown);
	GlobalAttach(down,up,todown,toup);

	ofstream os("topobftree2");
	GetTree()->ToStream(os);
	GlobalSwapTree();
	ofstream os1("topobftree1");
	GetTree()->ToStream(os1);

	if (((topobf == 2) || (topobf == 4)) && (bffrac >= 0))	{
		GlobalSwapTree();
	}
}

void PhyloProcess::SetSIS()	{

	GlobalSetMinMax(0,sisfrac);

}


void PhyloProcess::IncSize()	{
	size++;

	if (topobf == 1)	{
		if (size > bfburnin)	{

			double delta = 1.0 / bfnfrac;
			double deltalogp = GlobalComputeTopoBFLogLikelihoodRatio(bffrac,bffrac+delta);

			ofstream os((name + ".bf").c_str(),ios_base::app);
			os << bffrac << '\t' << deltalogp << '\n';
			os.close();

			int c = (size - bfburnin) % bfnrep;
			if (! c)	{
				bffrac += delta;
				GlobalSetBFFrac();
				GlobalUpdateConditionalLikelihoods();
			}
		}
	}

	if (topobf == 3)	{

		if (size > bfburnin)	{

			double delta = 1.0 / bfnfrac;
			double deltalogp1 = 0;
			double deltalogp2 = 0;
			if (bffrac + delta -1 < 1e-6)	{
				deltalogp1 = GlobalComputeTopoBFLogLikelihoodRatio(bffrac,bffrac+delta);
			}
			/*
			if (bffrac - delta > -1e-6)	{
				deltalogp2 = GlobalComputeTopoBFLogLikelihoodRatio(bffrac-delta,bffrac);
			}
			*/

			ofstream os((name + ".bf").c_str(),ios_base::app);
			os << bffrac << '\t' << deltalogp1 << '\n';
			// os << bffrac << '\t' << deltalogp1 << '\t' << deltalogp2 << '\n';
			os.close();
		}

	}

	if (topobf == 2)	{

		if (size > bfburnin)	{

			double deltalogp = 0;
			if (bffrac == -1)	{
				deltalogp = GlobalComputeTopoBFLogLikelihoodRatio(0,1);
			}
			else	{
				deltalogp -= LogLengthPrior();
				if (bffrac < -1)	{
					RescaleBranchPrior(blfactor,1);
				}
				else	{
					RescaleBranchPrior(1.0/blfactor,1);
				}
				deltalogp += LogLengthPrior();
				if (bffrac < -1)	{
					RescaleBranchPrior(1.0/blfactor,1);
				}
				else	{
					RescaleBranchPrior(blfactor,1);
				}
			}

			ofstream os((name + ".bf").c_str(),ios_base::app);
			os << bffrac << '\t' << deltalogp << '\t' << GetTotalNsite() * GetAllocTotalLength(1) << '\n';
			os.close();

			int c = (size - bfburnin) % bfnrep;
			if (! c)	{

				bffrac += 1;
				if (bffrac == 0)	{
					GlobalSwapTree();
					GlobalUpdateConditionalLikelihoods();
				}
				else if (bffrac < 0)	{
					RescaleBranchPrior(blfactor,1);
				}
				else	{
					RescaleBranchPrior(1.0/blfactor,1);
				}
			}

		}

	}
	if (topobf == 4)	{

		if (size > bfburnin)	{

			double deltalogp = 0;
			if (bffrac == -1)	{
				deltalogp = GlobalComputeTopoBFLogLikelihoodRatio(0,1);
			}
			else	{
				deltalogp -= LogLengthPrior();
				if (bffrac < -1)	{
					RescaleBranchPrior(blfactor,1);
				}
				else	{
					RescaleBranchPrior(1.0/blfactor,1);
				}
				deltalogp += LogLengthPrior();
				if (bffrac < -1)	{
					RescaleBranchPrior(1.0/blfactor,1);
				}
				else	{
					RescaleBranchPrior(blfactor,1);
				}
			}

			ofstream os((name + ".bf").c_str(),ios_base::app);
			os << bffrac << '\t' << deltalogp << '\t' << GetTotalNsite() * GetAllocTotalLength(1) << '\n';
			os.close();
		}

	}

	if (sis == 1)	{

		int fast = 0;
		if (sisfrac < 1.0)	{

			int c = size % sisnrep;
			if (c)	{
				if (! fast)	{
					reverseafterfull = 1;
					double delta = 1.0 / sisnfrac / sislevel;
					GlobalSetMinMax(sisfrac,sisfrac+delta);
					if (sumovercomponents > 0)	{
						GlobalChooseMultipleTryAlloc();
					}
					double deltalogp = GlobalGetFullLogLikelihood();
					reverseafterfull = 0;

					ofstream os((name + ".sis").c_str(),ios_base::app);
					os << sisfrac << '\t' << deltalogp << '\t' << logZ << '\n';
					os.close();
				}

			}
			else	{
				double delta = 1.0 / sisnfrac / sislevel;
				GlobalSetMinMax(sisfrac,sisfrac+delta);
				if (sumovercomponents > 0)	{
					GlobalChooseMultipleTryAlloc();
				}
				double deltalogp = GlobalGetFullLogLikelihood();
				logZ += deltalogp;

				ofstream os((name + ".sis").c_str(),ios_base::app);
				os << sisfrac << '\t' << deltalogp << '\t' << logZ << '\n';
				os.close();

				sisfrac += delta;
				if (sisfrac > 1.0)	{
					sisfrac = 1.0;
				}
				if ((sislevel == 10) && (sisfrac >= siscutoff))	{
					sislevel = 1;
				}
			}
		}
		GlobalSetMinMax(0,sisfrac);
		GlobalUpdateConditionalLikelihoods();
	}

	if (sis == 2)	{

		reverseafterfull = 1;
		double delta = 1.0 / sisnfrac / sislevel;
		GlobalSetMinMax(sisfrac,sisfrac+delta);
		if (sumovercomponents > 0)	{
			GlobalChooseMultipleTryAlloc();
		}
		double deltalogp = GlobalGetFullLogLikelihood();
		reverseafterfull = 0;

		ofstream os((name + ".sis").c_str(),ios_base::app);
		os << sisfrac << '\t' << deltalogp << '\n';
		os.close();

		GlobalSetMinMax(0,sisfrac);
		GlobalUpdateConditionalLikelihoods();
	}
}

void PhyloProcess::QuickUpdate()	{

	MPI_Status stat;
	MESSAGE signal = BCAST_TREE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	GlobalBroadcastTree();
	GlobalCollapse();
	GlobalUnfold();
}


/*
double PhyloProcess::ComputeBLLogLikelihoodRatio(double bffrac)	{

	double deltalogp = 0;
	if (bffrac == -1)	{
		deltalogp = GlobalComputeTopoBFLogLikelihoodRatio(0,1);
	}
	else	{
		deltalogp -= LogLengthPrior();
		if (bffrac < -1)	{
			RescaleBranchPrior(blfactor,1);
		}
		else	{
			RescaleBranchPrior(1.0/blfactor,1);
		}
		deltalogp += LogLengthPrior();
		if (bffrac < -1)	{
			RescaleBranchPrior(1.0/blfactor,1);
		}
		else	{
			RescaleBranchPrior(blfactor,1);
		}
	}
	return deltalogp;
}
*/

void PhyloProcess::ToStreamHeader(ostream& os)	{
	os << version << '\n';
	propchrono.ToStream(os);
	chronototal.ToStream(os);
	os << size << '\n';
	os << datafile << '\n';
	os << iscodon << '\n';
	os << codetype << '\n';
	os << fixtopo << '\n';
	os << fixroot << '\n';
	os << topoburnin << '\n';
	os << NSPR << '\t' << NMHSPR << '\t' << NTSPR << '\n';
	os << temperedbl << '\t' << temperedgene << '\t' << temperedrate << '\n';
	os << topolambda << '\t' << topomu << '\t' << toponstep << '\n';
	os << NNNI << '\n';
	os << nspec << '\t' << ntspec << '\n';
	os << taxon1 << '\t' << taxon2 << '\t' << taxon3 << '\t' << taxon4 << '\n';
	os << bpp << '\t' << nbpp << '\t' << ntbpp << '\t' << bppnstep << '\t' << bppname << '\t' << bppcutoff << '\t' << bppbeta << '\n';
	os << dc << '\n';
	os << fixbl << '\n';
	os << proposemode << '\n';
	os << allocmode << '\n';
	// os << sumratealloc << '\n';
	os << fasttopo << '\t' << fasttopofracmin << '\t' << fasttoponstep << '\n';
	os << fastcondrate << '\n';
	os << sumovercomponents << '\n';
	if (((topobf == 2) || (topobf == 4)) && (bffrac >= 0))	{
		GlobalSwapTree();
	}
	SetNamesFromLengths();
	GetTree()->ToStream(os);
	if (((topobf == 2) || (topobf == 4)) && (bffrac >= 0))	{
		GlobalSwapTree();
	}
	os << topobf << '\t' << bfburnin << '\t' << bfnfrac << '\t' << bfnrep << '\t' << bffrac << '\n';
	os << treefile << '\n';
	if ((topobf == 2) || (topobf == 4))	{
		os << blfactor << '\n';
		os << blfile << '\n';
	}
	os << sis << '\t' << sisnfrac << '\t' << sisnrep << '\t' << sisfrac << '\t' << siscutoff << '\t' << logZ << '\n';
	os << 0 << '\n';
	os << dirpriortype << '\t' << Nstatcomp << '\t' << priorempmix << '\t' << priormixtype << '\t' << fixstatweight << '\t' << fixstatalpha << '\t' << fixstatcenter << '\n';
	os << 0 << '\n';
	os << partitionfile << '\n';
	os << 1 << '\n';
}

void PhyloProcess::FromStreamHeader(istream& is)	{
	is >> version;
	if (atof(version.substr(0,3).c_str()) < 1.2)	{
		cerr << "error: version is too old : " << version << '\n';
		exit(1);
	}
	propchrono.FromStream(is);
	chronototal.FromStream(is);
	string indatafile;
	is >> size;
	is >> datafile;
	is >> iscodon;
	is >> codetype;
	is >> fixtopo;
	is >> fixroot;
	is >> topoburnin;
	is >> NSPR >> NMHSPR >> NTSPR;
	is >> temperedbl >> temperedgene >> temperedrate;
	is >> topolambda >> topomu >> toponstep;
	is >> NNNI;
	is >> nspec >> ntspec;
	is >> taxon1 >> taxon2 >> taxon3 >> taxon4;
	is >> bpp >> nbpp >> ntbpp >> bppnstep >> bppname >> bppcutoff >> bppbeta;
	is >> dc;
	is >> fixbl;
	is >> proposemode;
	is >> allocmode;
	// is >> sumratealloc;
	is >> fasttopo >> fasttopofracmin >> fasttoponstep;
	is >> fastcondrate;
	is >> sumovercomponents;
	is >> treestring;
	is >> topobf >> bfburnin >> bfnfrac >> bfnrep >> bffrac;
	is >> treefile ;
	if ((topobf == 2) || (topobf == 4))	{
		is >> blfactor;
		is >> blfile;
	}
	is >> sis >> sisnfrac >> sisnrep >> sisfrac >> siscutoff >> logZ;
	int check;
	is >> check;
	if (!check)	{
		is >> dirpriortype >> Nstatcomp >> priorempmix >> priormixtype >> fixstatweight >> fixstatalpha >> fixstatcenter;
		is >> check;
		partitionfile = "None";
		if (! check)	{
			is >> partitionfile;
			is >> check;
			if (! check)	{
				cerr << "error when reading stream header \n";
				exit(1);
			}
		}
	}
}

//-------------------------------------------------------------------------
//	* Create / Delete
//-------------------------------------------------------------------------

void PhyloProcess::Create()	{

	if (! empfreq)	{
		if (! data)	{
			cerr << "error in PhyloProcess::Create: data have not been specified\n";
			exit(1);
		}
		RateProcess::Create();
		ProfileProcess::Create();
		BranchProcess::Create();

		empfreq = new double[GetData()->GetNstate()];
		GetData()->GetEmpiricalFreq(empfreq);

		loglarray = new double[GetNbranch()];
		// MPI : slaves only
		// for each slave, should specify the range of sites (sitemin <= i < sitemax)
		if ((GetMyid() > 0) || (GetNprocs() == 1)) {
			SubstitutionProcess::Create();

			nodestate = new int*[GetNnode()];
			condlmap = new double***[GetNlink()];
			for (int j=0; j<GetNlink(); j++)	{
				condlmap[j] = 0;
			}

            if (withoutwardcondlmap)    {
                condlmap2 = new double***[GetNlink()];
                for (int j=0; j<GetNlink(); j++)	{
                    condlmap2[j] = 0;
                }
            }

			CreateCondSiteLogL();
			ActivateSumOverRateAllocations();
			CreateConditionalLikelihoods();
			CreateNodeStates();
			CreateMappings();
			CreateSuffStat();
			CreateMissingMap();
			activesuffstat = false;
		}
		else {
			nodestate = new int*[GetNnode()];
			CreateCondSiteLogL();
			ActivateSumOverRateAllocations();
			CreateNodeStates();
			CreateSuffStat();
			CreateMissingMap();
			activesuffstat = false;
		}
		CreateSiteConditionalLikelihoods();
	}
}

void PhyloProcess::CreateMissingMap()	{

	missingmap = new int*[GetNbranch()];
	for (int j=0; j<GetNbranch(); j++)	{
		missingmap[j] = new int[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			missingmap[j][i] = -1;
		}
	}
}

void PhyloProcess::DeleteMissingMap()	{

	for (int j=0; j<GetNbranch(); j++)	{
		delete[] missingmap[j];
	}
	delete[] missingmap;
}

void PhyloProcess::FillMissingMap(int reduced)	{
    if (! reduced)   {
        for (int i=0; i<GetNsite(); i++)	{
            missingmap[0][i] = 2;
        }
        for (int j=1; j<GetNbranch(); j++)	{
            for (int i=0; i<GetNsite(); i++)	{
                missingmap[j][i] = 1;
            }
        }
    }
    else    {
        BackwardFillMissingMap(GetRoot());
        ForwardFillMissingMap(GetRoot(),GetRoot());
    }
}

void PhyloProcess::BackwardFillMissingMap(const Link* from)	{

	int index = GetBranchIndex(from->GetBranch());
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			missingmap[index][i] = 0;
		}
	}
	if (from->isLeaf())	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				int state = GetData(from)[i];
				if (state != -1)	{
					missingmap[index][i] = 1;
				}
			}
		}
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			BackwardFillMissingMap(link->Out());
			int j = GetBranchIndex(link->Out()->GetBranch());
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				if (ActiveSite(i))	{
					if (missingmap[j][i])	{
						missingmap[index][i] ++;
					}
				}
			}
		}
	}
}

void PhyloProcess::ForwardFillMissingMap(const Link* from, const Link* up)	{

	int index = GetBranchIndex(from->GetBranch());
	int upindex = GetBranchIndex(up->GetBranch());
	if (from->isRoot())	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				if (missingmap[index][i] <= 1)	{
					missingmap[index][i] = 0;
				}
				else	{
					missingmap[index][i] = 2;
				}
			}
		}
	}
	else	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				if (missingmap[index][i] > 0)	{
					if (missingmap[upindex][i])	{
						missingmap[index][i] = 1;
					}
					else	{
						if (from->isLeaf() || (missingmap[index][i] > 1))	{
							missingmap[index][i] = 2;
						}
						else	{
							missingmap[index][i] = 0;
						}
					}
				}
			}
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		ForwardFillMissingMap(link->Out(),from);
	}
}

void PhyloProcess::Delete() {

	if (data)	{
		// MPI slaves only
		if ((GetMyid() > 0) || (GetNprocs() == 1)) {
			DeleteConditionalLikelihoods();
			DeleteCondSiteLogL();
			DeleteNodeStates();
			FullDeleteMappings();
			DeleteSuffStat();
			DeleteMissingMap();
			delete[] nodestate;
			delete[] condlmap;
            delete[] condlmap2;
			SubstitutionProcess::Delete();
		}
		else	{
			DeleteNodeStates();
			delete[] nodestate;
			DeleteCondSiteLogL();
			DeleteSuffStat();
			DeleteMissingMap();
		}
		// MPI master and slaves
		BranchProcess::Delete();
		ProfileProcess::Delete();
		RateProcess::Delete();

		delete[] empfreq;
		DeleteSiteConditionalLikelihoods();
	}
}

//-------------------------------------------------------------------------
//	* Create / Delete Mappings
//-------------------------------------------------------------------------

void PhyloProcess::CreateMappings()	{

	submap = new BranchSitePath**[GetNbranch()];
	for (int j=0; j<GetNbranch(); j++)	{
		submap[j] = new BranchSitePath*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			submap[j][i] = 0;
		}
	}
}

void PhyloProcess::DeleteMappings()	{

	for (int j=0; j<GetNbranch(); j++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				delete submap[j][i];
				submap[j][i] = 0;
			}
		}
	}
}

void PhyloProcess::FullDeleteMappings()	{

	DeleteMappings();
	for (int j=0; j<GetNbranch(); j++)	{
		delete[] submap[j];
	}
	delete[] submap;
}

//-------------------------------------------------------------------------
//	* Create / Delete SuffStat
//-------------------------------------------------------------------------

void PhyloProcess::CreateSuffStat()	{

	if (! siteratesuffstatcount)	{
		siteratesuffstatcount = new double[GetNsite()];
		siteratesuffstatbeta = new double[GetNsite()];
	}
	if (! branchlengthsuffstatcount)	{
		branchlengthsuffstatcount = new double[GetNbranch()];
		branchlengthsuffstatbeta = new double[GetNbranch()];
	}
}

void PhyloProcess::DeleteSuffStat()	{

	delete[] siteratesuffstatcount;
	delete[] siteratesuffstatbeta;
	siteratesuffstatcount = 0;
	siteratesuffstatbeta = 0;
	delete[] branchlengthsuffstatcount;
	delete[] branchlengthsuffstatbeta;
	branchlengthsuffstatcount = 0;
	branchlengthsuffstatbeta = 0;
}

//-------------------------------------------------------------------------
//	* Create / Delete NodeStates
//-------------------------------------------------------------------------

void PhyloProcess::CreateNodeStates()	{

	for (int j=0; j<GetNnode(); j++)	{
		nodestate[j] = new int[GetNsite()];
	}
}

void PhyloProcess::DeleteNodeStates()	{

	for (int j=0; j<GetNnode(); j++)	{
		delete[] nodestate[j];
	}
}

//-------------------------------------------------------------------------
//	* Create / Delete Conditional Likelihoods
//-------------------------------------------------------------------------

void PhyloProcess::CreateConditionalLikelihoods()	{

	// do not create for leaves
	for (int j=0; j<GetNlink(); j++)	{
		if (condlmap[j])	{
			cerr << "error: condl already exists\n";
			exit(1);
		}
		condlmap[j] =  CreateConditionalLikelihoodVector();
        if (withoutwardcondlmap)    {
            if (condlmap2[j])   {
                cerr << "error: condl2 already exists\n";
                exit(1);
            }
            condlmap2[j] =  CreateConditionalLikelihoodVector();
        }
	}
}

void PhyloProcess::DeleteConditionalLikelihoods()	{

	for (int j=0; j<GetNlink(); j++)	{
		DeleteConditionalLikelihoodVector(condlmap[j]);
		condlmap[j] = 0;
        if (withoutwardcondlmap)    {
            DeleteConditionalLikelihoodVector(condlmap2[j]);
            condlmap2[j] = 0;
        }
	}
}

//-------------------------------------------------------------------------
//	* Broadcast Tree
//-------------------------------------------------------------------------

void PhyloProcess::GlobalBroadcastTree()	{

	if (GetNprocs() > 1)	{
		ostringstream os;
		tree->ToStream(os);
		string s = os.str();
		unsigned int len = s.length();
		unsigned char* bvector = new unsigned char[len];
		for (unsigned int i=0; i<len; i++)	{
			bvector[i] = s[i];
		}
		MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
		delete[] bvector;
	}
}

void PhyloProcess::SlaveBroadcastTree()	{

	int len;
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	unsigned char* bvector = new unsigned char[len];
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	ostringstream os;
	for (int i=0; i<len; i++)	{
		os << bvector[i];
	}
	istringstream is(os.str());
	tree->ReadFromStream(is);
	delete[] bvector;
}

//-------------------------------------------------------------------------
//	* Unfold
//-------------------------------------------------------------------------

void PhyloProcess::Unfold()	{

	DeleteMappings();
	ActivateSumOverRateAllocations();
	/*
	if (!sumratealloc)	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations();
	}
	*/
	activesuffstat = false;
}

void PhyloProcess::GlobalUnfold()	{

	if (GetNprocs() > 1)	{
		GlobalUpdateParameters();

		MESSAGE signal = UNFOLD;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		GlobalUpdateConditionalLikelihoods();
	}
	else	{
		Unfold();
		UpdateConditionalLikelihoods();
        /*
        cerr << "check likelihoods\n";
        CheckLikelihood();
        exit(1);
        */
	}
}

void PhyloProcess::SlaveUnfold()	{
	Unfold();
}

//-------------------------------------------------------------------------
//	* Collapse
//-------------------------------------------------------------------------

void PhyloProcess::Collapse()	{

	// if (sumratealloc)	{
	DrawAllocations(0);
	InactivateSumOverRateAllocations();
	// }
	if ((topobf == 1) || (topobf == 3))	{
		SetMinMax(bffrac,1);
		SampleNodeStates();
		FillMissingMap();
		SampleSubstitutionMappings(GetRoot());
		SlaveSwapTree();
		SetMinMax(0,bffrac);
		SampleNodeStates();
		FillMissingMap();
		SampleSubstitutionMappings(GetRoot());
		SlaveSwapTree();
		SetMinMax(0,1);
	}
	else	{
		SampleNodeStates();
		FillMissingMap();
		SampleSubstitutionMappings(GetRoot());
	}
	activesuffstat = true;
}

void PhyloProcess::GlobalCollapse()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = COLLAPSE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		Collapse();
	}
}

void PhyloProcess::SlaveCollapse()	{
	Collapse();
}

//-------------------------------------------------------------------------
//	* Activate / Inactivate sum over rate categories
//-------------------------------------------------------------------------
// local functions are defined in RateProcess

void PhyloProcess::GlobalActivateSumOverRateAllocations()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = SUMRATE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		ActivateSumOverRateAllocations();
	}
}

void PhyloProcess::SlaveActivateSumOverRateAllocations()	{
	ActivateSumOverRateAllocations();
	// sumratealloc = 1;
}

void PhyloProcess::GlobalInactivateSumOverRateAllocations()	{

	// assumes cond likelihoods already updated
	if (GetNprocs() > 1)	{
		MESSAGE signal = CONDRATE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations();
	}
}

void PhyloProcess::SlaveInactivateSumOverRateAllocations()	{
	DrawAllocations(0);
	InactivateSumOverRateAllocations();
	// sumratealloc = 0;
}

//-------------------------------------------------------------------------
//	* Activate / Inactivate Zip
//-------------------------------------------------------------------------

void PhyloProcess::GlobalActivateZip()	{

	MESSAGE signal = ACTIVATEZIP;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	ActivateZip();
}

void PhyloProcess::GlobalInactivateZip()	{

	MESSAGE signal = INACTIVATEZIP;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	InactivateZip();
}

void PhyloProcess::SlaveActivateZip()	{
	ActivateZip();
}

void PhyloProcess::SlaveInactivateZip()	{
	InactivateZip();
}

//-------------------------------------------------------------------------
//	* Likelihood computation
//-------------------------------------------------------------------------

double PhyloProcess::ComputeTopoBFLogLikelihoodRatio(double fracmin, double fracmax)	{

	double delta = 0;

	if (sumovercomponents)	{

		if (sumovercomponents != -1)	{
			ChooseMultipleTryAlloc();
		}

		SetMinMax(fracmin,fracmax);

		delta -= GetFullLogLikelihood();
		SlaveSwapTree();
		delta += GetFullLogLikelihood();
		SlaveSwapTree();

		SetMinMax(0,1);
	}
	else	{

		SetMinMax(fracmin,fracmax);

		PostOrderPruning(GetRoot(),condlmap[0]);
		MultiplyByStationaries(condlmap[0],condflag);
		ComputeLikelihood(condlmap[0],condflag);
		PreOrderPruning(GetRoot(),condlmap[0]);

		delta -= SumLogLikelihood();

		SlaveSwapTree();

		PostOrderPruning(GetRoot(),condlmap[0]);
		MultiplyByStationaries(condlmap[0],condflag);
		ComputeLikelihood(condlmap[0],condflag);
		PreOrderPruning(GetRoot(),condlmap[0]);

		delta += SumLogLikelihood();

		SlaveSwapTree();

		PostOrderPruning(GetRoot(),condlmap[0]);
		MultiplyByStationaries(condlmap[0],condflag);
		ComputeLikelihood(condlmap[0],condflag);
		PreOrderPruning(GetRoot(),condlmap[0]);

		SetMinMax(0,1);
	}

	return delta;
}

double PhyloProcess::GlobalComputeTopoBFLogLikelihoodRatio(double fracmin, double fracmax)	{

	double ret = 0;

	if (GetNprocs() > 1)	{

		double frac[2];
		frac[0] = fracmin;
		frac[1] = fracmax;
		MESSAGE signal = TOPOLRT;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(frac,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

		MPI_Status stat;
		double delta;
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(&delta,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			ret += delta;
		}
	}
	else	{
		ret = ComputeTopoBFLogLikelihoodRatio(fracmin,fracmax);
	}
	return ret;
}

void PhyloProcess::SlaveComputeTopoBFLogLikelihoodRatio()	{

	double frac[2];
	MPI_Bcast(frac,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	double delta = ComputeTopoBFLogLikelihoodRatio(frac[0],frac[1]);
	MPI_Send(&delta,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}


void PhyloProcess::GlobalSetBFFrac()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = SETBFFRAC;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&bffrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
}

void PhyloProcess::SlaveSetBFFrac()	{

	MPI_Bcast(&bffrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalSetSISFrac()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = SETSISFRAC;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sisfrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
}

void PhyloProcess::SlaveSetSISFrac()	{

	MPI_Bcast(&sisfrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void PhyloProcess::UpdateConditionalLikelihoods()	{

	if ((topobf == 1) || (topobf == 3))	{

		SetMinMax(bffrac,1);

		PostOrderPruning(GetRoot(),condlmap[0]);
		MultiplyByStationaries(condlmap[0],condflag);
		ComputeLikelihood(condlmap[0],condflag);
		PreOrderPruning(GetRoot(),condlmap[0]);

		SlaveSwapTree();
		SetMinMax(0,bffrac);

		PostOrderPruning(GetRoot(),condlmap[0]);
		MultiplyByStationaries(condlmap[0],condflag);
		ComputeLikelihood(condlmap[0],condflag);
		PreOrderPruning(GetRoot(),condlmap[0]);

		SlaveSwapTree();
		SetMinMax(0,1);
		SumLogLikelihood();
	}
	else	{
		PostOrderPruning(GetRoot(),condlmap[0]);
		MultiplyByStationaries(condlmap[0],condflag);
		ComputeLikelihood(condlmap[0],condflag);
		PreOrderPruning(GetRoot(),condlmap[0]);
	}
}

void PhyloProcess::CheckLikelihood()	{

	vector<double> nodelogl;
	vector<double> branchlogl;
	RecursiveComputeLikelihood(GetRoot(),-1,nodelogl,branchlogl);
	double max = 0;
	for (unsigned int i=0; i<nodelogl.size(); i++)	{
		double tmp = fabs(nodelogl[i] - nodelogl[0]);
		if (max < tmp)	{
			max = tmp;
		}
	}
	for (unsigned int i=0; i<branchlogl.size(); i++)	{
		double tmp = fabs(branchlogl[i] - branchlogl[0]);
		if (max < tmp)	{
			max = tmp;
		}
	}
	if (max > 1e-10)	{
		cerr << "error in check likelihoods\n";
		cerr << max << '\n';
		cerr.precision(25);
        cerr << "node log likelihoods\n";
		for (unsigned int i=0; i<nodelogl.size(); i++)	{
			cerr << nodelogl[i] << '\n';
		}
        cerr << "branch log likelihoods\n";
		for (unsigned int i=0; i<branchlogl.size(); i++)	{
			cerr << branchlogl[i] << '\n';
		}
		exit(1);
	}
}

void PhyloProcess::RecursiveComputeLikelihood(const Link* from, int auxindex, vector<double>& nodelogl, vector<double>& branchlogl)	{

    // WARNING: preorder pruning does not update leaf condtional likelihood vectors (not necessary in the present case)
    // so the following will issue an error message if tried on leaf
    if (! from->isLeaf())   {
        double nodelnL = ComputeNodeLikelihood(from,auxindex);
        nodelogl.push_back(nodelnL);
    }

    if (! from->isRoot())   {
        double branchlnL = ComputeBranchLikelihood(from,auxindex);
        branchlogl.push_back(branchlnL);
    }

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
        RecursiveComputeLikelihood(link->Out(),auxindex,nodelogl,branchlogl);
	}
}

void PhyloProcess::GlobalUpdateConditionalLikelihoods()	{

	if (GetNprocs() > 1)	{
		// MPI
		// just send Updateconlikelihood message to all slaves
		MESSAGE signal = UPDATE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalCollectLogLikelihood();
	}
	else	{
		UpdateConditionalLikelihoods();
		// SumLogLikelihood();
		// ComputeNodeLikelihood(GetRoot(),0);
	}
}

void PhyloProcess::SlaveUpdateConditionalLikelihoods()	{
	UpdateConditionalLikelihoods();
}

double PhyloProcess::ComputeNodeLikelihood(const Link* from, int auxindex)	{

	double*** aux = 0;
	bool localaux = false;
	if (auxindex != -1)	{
		aux = condlmap[auxindex];
	}
	else	{
		localaux = true;
		aux = CreateConditionalLikelihoodVector();
	}

	if (from->isLeaf())	{
		Initialize(aux,GetData(from),condflag);
	}
	else	{
		Reset(aux,condflag);
	}

	if (! from->isRoot())	{
		Multiply(GetConditionalLikelihoodVector(from),aux,condflag);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! link->isRoot())	{
			Multiply(GetConditionalLikelihoodVector(link),aux,condflag);
		}
	}
	MultiplyByStationaries(aux,condflag);
	double lnL = ComputeLikelihood(aux,condflag);
	if (localaux)	{
		DeleteConditionalLikelihoodVector(aux);
	}
	return lnL;
}

double PhyloProcess::ComputeBranchLikelihood(const Link* from, int auxindex)	{

    if (from->isRoot()) {
        cerr << "error: compute branch likelihood called on root\n";
        exit(1);
    }

	double*** aux = 0;
	bool localaux = false;
	if (auxindex != -1)	{
		aux = condlmap[auxindex];
	}
	else	{
		localaux = true;
		aux = CreateConditionalLikelihoodVector();
	}

    Propagate(GetOutwardConditionalLikelihoodVector(from->Out()),aux,GetLength(from->GetBranch()),condflag);
    Multiply(GetOutwardConditionalLikelihoodVector(from),aux,condflag);
	MultiplyByStationaries(aux,condflag);
	double lnL = ComputeLikelihood(aux,condflag);
	if (localaux)	{
		DeleteConditionalLikelihoodVector(aux);
	}
	return lnL;
}

double PhyloProcess::GlobalCollectLogLikelihood()	{

	if (GetNprocs() > 1)	{

		MESSAGE signal = COLLECTLIKELIHOOD;
		MPI_Status stat;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		logL = 0.0;
		double sum;
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(&sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			logL += sum;
		}
	}
	else	{
		// logL = SumLogLikelihood();
	}
	if (isnan(logL))	{
		cerr << "in PhyloProcess::GlobalCollectLogLikelihood: logL is nan\n";
		exit(1);
	}
	return logL;
}

void PhyloProcess::SlaveCollectLogLikelihood()	{
	MPI_Send(&logL,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	/*
	double ret = SumLogLikelihood();
	MPI_Send(&ret,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	*/
}

double PhyloProcess::SumLogLikelihood()	{

	logL = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			logL += sitelogL[i];
		}
	}
	return logL;
}

double PhyloProcess::GlobalComputeNodeLikelihood(const Link* from, int auxindex)	{ 

	if (GetNprocs() > 1)	{

		MESSAGE signal = LIKELIHOOD;
		MPI_Status stat;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		int args[] = {GetLinkIndex(from),auxindex};
		MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
		// master : sums up all values sent by slaves
		// store this sum into member variable logL
		// and return it

		logL = 0.0;
		double sum;
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(&sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			logL += sum;
		}
	}
	else	{
		logL = ComputeNodeLikelihood(from,auxindex);
	}
	if (isnan(logL))	{
		cerr << "in PhyloProcess::GlobalComputeNodeLikelihood: logL is nan\n";
		exit(1);
	}
	return logL;
}

void PhyloProcess::SlaveComputeNodeLikelihood(int fromindex,int auxindex) {
	double ret = LocalComputeNodeLikelihood(fromindex,auxindex);
	MPI_Send(&ret,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double PhyloProcess::LocalComputeNodeLikelihood(int fromindex,int auxindex) {
	return ComputeNodeLikelihood(GetLinkForGibbs(fromindex),auxindex);
}

double PhyloProcess::GlobalGetFullLogLikelihood()	{

	double totlogL = 0;

	if (! sumovercomponents)	{
		GlobalUpdateConditionalLikelihoods();
		totlogL = logL;
	}

	else	{

		logL = 0;

		if (GetNprocs() > 1)	{

			MESSAGE signal = FULLLIKELIHOOD;
			MPI_Status stat;
			MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&reverseafterfull,1,MPI_INT,0,MPI_COMM_WORLD);

			double sum[2];
			for(int i=1; i<GetNprocs(); i++) {
				MPI_Recv(&sum,2,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
				totlogL += sum[0];
				logL += sum[1];
			}
		}
		else	{
			totlogL = GetFullLogLikelihood();
		}
	}
	if (isnan(logL))	{
		cerr << "in PhyloProcess::GlobalGetFullLogLikelihood: logL is nan\n";
		exit(1);
	}
	if (isnan(totlogL))	{
		cerr << "in PhyloProcess::GlobalGetFullLogLikelihood: totlogL is nan\n";
		exit(1);
	}
	return totlogL;
}

void PhyloProcess::SlaveGetFullLogLikelihood()	{

	MPI_Bcast(&reverseafterfull,1,MPI_INT,0,MPI_COMM_WORLD);
	double totlogL = GetFullLogLikelihood();
	double sum[2];
	sum[0] = totlogL;
	sum[1] = logL;
	MPI_Send(&totlogL,2,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::PostOrderPruning(const Link* from, double*** aux)	{

	if (from->isLeaf())	{
		Initialize(aux,GetData(from),condflag);
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			PostOrderPruning(link->Out(),aux);
			Propagate(aux,GetConditionalLikelihoodVector(link),GetLength(link->GetBranch()),condflag);
		}
		Reset(aux,condflag);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			Multiply(GetConditionalLikelihoodVector(link),aux,condflag);
		}
		Offset(aux,condflag);
	}
    if (withoutwardcondlmap)    {
        Reset(GetOutwardConditionalLikelihoodVector(from),condflag);
        Multiply(aux,GetOutwardConditionalLikelihoodVector(from),condflag);
    }
}

void PhyloProcess::PreOrderPruning(const Link* from, double*** aux)	{

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{

		Reset(aux,condflag);
		for (const Link* link2=link->Next(); link2!=link; link2=link2->Next())	{
			if (! link2->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link2),aux,condflag);
			}
		}
		// Here, in principle
		// should be done even if link->Out()->isLeaf()
		// in order for all the conditional likelihood vectors, including those at the leaves, to be updated
		// but in practice, the leaf likelihood vectors are not used anyway (and they represent half of the whole set of likelihood vectors)
		// so not computing them saves 50% CPU time

        if (withoutwardcondlmap)    {
            Reset(GetOutwardConditionalLikelihoodVector(link),condflag);
            Multiply(aux,GetOutwardConditionalLikelihoodVector(link),condflag);
        }

		if (! link->Out()->isLeaf())	{
			Propagate(aux,GetConditionalLikelihoodVector(link->Out()),GetLength(link->GetBranch()),condflag);
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! link->Out()->isLeaf())	{
			PreOrderPruning(link->Out(),aux);
		}
	}
}

void PhyloProcess::GlobalReset(const Link* link, bool condalloc)	{

	if (GetNprocs() > 1)	{
		// MPI
		// send a Reset message with GetLinkIndex(link) as argument
		// slaves: upon receiving message
		// call the Reset function with link corresponding to index received as argument of the message
		MESSAGE signal = RESET;
		int args[2];
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		args[0] = GetLinkIndex(link);
		args[1] = (condalloc) ? 1 : 0;
		MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		Reset(condlmap[GetLinkIndex(link)],condalloc);
	}
}

void PhyloProcess::SlaveReset(int n,bool v) {
	Reset(condlmap[n],v);	
}

void PhyloProcess::GlobalMultiply(const Link* from, const Link* to, bool condalloc)	{

	if (GetNprocs() > 1)	{
		// MPI
		// send a Multiply message with GetLinkIndex(from) and GetLinkIndex(to) as argument
		// slaves: upon receiving message
		// call the Multiply function with links corresponding to the two indices received as argument
		MESSAGE signal = MULTIPLY;
		int args[3];
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		args[0] = GetLinkIndex(from);
		args[1] = GetLinkIndex(to);
		args[2] = (condalloc) ? 1 : 0;
		MPI_Bcast(args,3,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		Multiply(condlmap[GetLinkIndex(from)],condlmap[GetLinkIndex(to)],condalloc);
	}
}

void PhyloProcess::SlaveMultiply(int n,int m,bool v) {
	Multiply(condlmap[n],condlmap[m],v);
}

void PhyloProcess::GlobalMultiplyByStationaries(const Link* from, bool condalloc)	{

	if (GetNprocs() > 1)	{
		// MPI
		MESSAGE signal = SMULTIPLY;
		int args[2];
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		args[0] = GetLinkIndex(from);
		args[1] = (condalloc) ? 1 : 0;
		MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		MultiplyByStationaries(condlmap[GetLinkIndex(from)],condalloc);
	}
}

void PhyloProcess::SlaveMultiplyByStationaries(int n,bool v) {
	MultiplyByStationaries(condlmap[n],v);
}

void PhyloProcess::GlobalInitialize(const Link* from, const Link* link, bool condalloc)	{

	if (GetNprocs() > 1)	{
		// MPI
		MESSAGE signal = INITIALIZE;
		int args[3];
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		args[0] = GetLinkIndex(from);
		args[1] = GetLinkIndex(link);
		args[2] = (condalloc) ? 1 : 0;
		MPI_Bcast(args,3,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		Initialize(condlmap[GetLinkIndex(from)],GetData(link),condalloc);
	}
}

void PhyloProcess::SlaveInitialize(int n,int m,bool v) {
	const Link* link = GetLink(m);
	Initialize(condlmap[n],GetData(link),v);
}


void PhyloProcess::GlobalPropagate(const Link* from, const Link* to, double time, bool condalloc)	{

	if (GetNprocs() > 1)	{
		// MPI
		MESSAGE signal = PROPAGATE;
		prop_arg args;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		args.from = GetLinkIndex(from);
		args.to = GetLinkIndex(to);
		args.condalloc = (condalloc) ? 1 : 0;
		args.time = time;
		MPI_Bcast(&args,1,Propagate_arg,0,MPI_COMM_WORLD);
	}
	else	{
		Propagate(condlmap[GetLinkIndex(from)],condlmap[GetLinkIndex(to)],time,condalloc);
	}
}

void PhyloProcess::SlavePropagate(int n,int m,bool v,double t) {
	Propagate(condlmap[n],condlmap[m],t,v);
	Offset(condlmap[m],v);
}

//-------------------------------------------------------------------------
//	* Sample Node States
//-------------------------------------------------------------------------

void PhyloProcess::SampleNodeStates()	{
	SampleNodeStates(GetRoot(),condlmap[0]);
}


void PhyloProcess::SampleNodeStates(const Link* from, double*** aux)	{
	
	if (from->isLeaf())	{
		Initialize(aux,GetData(from),condflag);
	}
	else	{
		Reset(aux,true);
	}
	// make product of conditional likelihoods around node
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Multiply(GetConditionalLikelihoodVector(link),aux,true);
	}
	if (!from->isRoot())	{
		Multiply(GetConditionalLikelihoodVector(from),aux,true);
	}
	MultiplyByStationaries(aux,true);
	// let substitution process choose states based on this vector
	// this should collapse the vector into 1s and 0s
	ChooseStates(aux,GetStates(from->GetNode()));

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// propagate forward
		Propagate(aux,GetConditionalLikelihoodVector(link->Out()),GetLength(link->GetBranch()),true);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleNodeStates(link->Out(),aux);
	}
}

//-------------------------------------------------------------------------
//	* Sample Substitution Mappings
//-------------------------------------------------------------------------

void PhyloProcess::SampleSiteSubstitutionMapping(int site, const Link* from)	{

	if (from->isRoot())	{
		submap[0][site] = SampleRootSitePath(site,GetStates(from->GetNode())[site]);
	}
	else	{
		submap[GetBranchIndex(from->GetBranch())][site] = SampleSitePath(site,GetStates(from->Out()->GetNode())[site], GetStates(from->GetNode())[site], GetLength(from->GetBranch()));
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleSiteSubstitutionMapping(site,link->Out());
	}
}

void PhyloProcess::SampleSubstitutionMappings(const Link* from)	{

	if (from->isRoot())	{
		SampleRootPaths(submap[0],GetStates(from->GetNode()));
	}
	else	{
		SamplePaths(submap[GetBranchIndex(from->GetBranch())],GetStates(from->Out()->GetNode()), GetStates(from->GetNode()), GetLength(from->GetBranch()));
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleSubstitutionMappings(link->Out());
	}
}

//-------------------------------------------------------------------------
//	* BranchLength Moves 
//-------------------------------------------------------------------------

// MPI master functions
double PhyloProcess::BranchLengthMove(double tuning)	{

	if (GetNprocs() == 1)	{
		return NonMPIBranchLengthMove(tuning);
	}
	// uses condlmap[0] as auxiliary variable
	int n = 0;
	double total = RecursiveBranchLengthMove(GetRoot(),tuning,n);
	return total / n;
}

// assumes aux contains the product of incoming likelihoods
double PhyloProcess::RecursiveBranchLengthMove(const Link* from, double tuning, int& n)	{

	// uses condlmap[0] as auxiliary variable
	double total = 0;

	if (! from->isRoot())	{
		total += LocalBranchLengthMove(from,tuning);
		n++;
	}
	
	vector<const Link*> v;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		v.push_back(link);
	}

	double x = rnd::GetRandom().Uniform();

	for (unsigned int i=0; i<v.size(); i++)	{
		const Link* link = 0;
		if (x < 0.5)	{
			link = v[i];
		}
		else	{
			link = v[v.size()-1-i];
		}
	// for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		GlobalReset(0);
		// Reset(condlmap[0]);
		for (const Link* link2=link->Next(); link2!=link; link2=link2->Next())	{
			if (! link2->isRoot())	{
				GlobalMultiply(link2,0);
				// Multiply(GetConditionalLikelihoodVector(link2),condlmap[0]);
			}
		}
		total += RecursiveBranchLengthMove(link->Out(),tuning,n);
	}

	if (from->isLeaf())	{
		GlobalInitialize(0,from);
		// Initialize(condlmap[0],GetData(from));
	}
	else	{
		GlobalReset(0);
		// Reset(condlmap[0]);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (! link->isRoot())	{
				GlobalMultiply(link,0);
				// Multiply(GetConditionalLikelihoodVector(link),condlmap[0]);
			}
		}
	}
	
	if (! from->isRoot())	{
		total += LocalBranchLengthMove(from->Out(),tuning);
		n++;
	}

	return total;
};

double PhyloProcess::LocalBranchLengthMove(const Link* from, double tuning)	{

	// uses condlmap[0] as auxiliary variable

	double currentloglikelihood = logL;
	double currentlogprior = LogBranchLengthPrior(from->GetBranch());
	double loghastings = GlobalProposeMove(from->GetBranch(),tuning);
	// double loghastings = ProposeMove(from->GetBranch(),tuning);

	GlobalPropagate(0,from,GetLength(from->GetBranch()));
	// Propagate(aux,GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()));

	double newloglikelihood = GlobalComputeNodeLikelihood(from);
	// double newloglikelihood = ComputeNodeLikelihood(from);
	double newlogprior = LogBranchLengthPrior(from->GetBranch());
	double delta = newlogprior + newloglikelihood - currentlogprior - currentloglikelihood + loghastings;
	
	int accepted = (log(rnd::GetRandom().Uniform()) < delta);
	if (!accepted)	{
		GlobalRestoreBranch(from->GetBranch());
		GlobalPropagate(0,from,GetLength(from->GetBranch()));
		GlobalComputeNodeLikelihood(from);
		// not useful: done by ComputeNodeLikelihood(from) just above
		// logL = currentloglikelihood;
	}

	return (double) accepted;
}

// MPI master functions
double PhyloProcess::NonMPIBranchLengthMove(double tuning)	{

	UpdateConditionalLikelihoods();
	// uses condlmap[0] as auxiliary variable
	int n = 0;
	double total = RecursiveNonMPIBranchLengthMove(GetRoot(),tuning,n);
	return total / n;
}

// assumes aux contains the product of incoming likelihoods
double PhyloProcess::RecursiveNonMPIBranchLengthMove(const Link* from, double tuning, int& n)	{

	// uses condlmap[0] as auxiliary variable
	double total = 0;

	if (! from->isRoot())	{
		total += LocalNonMPIBranchLengthMove(from,tuning);
		n++;
	}
	//This update should be in the previous "if" ?
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Reset(condlmap[0],condflag);
		for (const Link* link2=link->Next(); link2!=link; link2=link2->Next())	{
			if (! link2->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link2),condlmap[0],condflag);
			}
		}
		total += RecursiveNonMPIBranchLengthMove(link->Out(),tuning,n);
	}

	if (from->isLeaf())	{
		Initialize(condlmap[0],GetData(from),condflag);
	}
	else	{
		Reset(condlmap[0],condflag);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (! link->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link),condlmap[0],condflag);
			}
		}
	}
	
	if (! from->isRoot())	{
		total += LocalNonMPIBranchLengthMove(from->Out(),tuning);
		n++;
	}

	return total;
};

double PhyloProcess::LocalNonMPIBranchLengthMove(const Link* from, double tuning)	{

	// uses condlmap[0] as auxiliary variable

	double currentloglikelihood = logL;
	double currentlogprior = LogBranchLengthPrior(from->GetBranch());
	double loghastings = ProposeMove(from->GetBranch(),tuning);

	Propagate(condlmap[0],GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()),condflag);

	double newloglikelihood = ComputeNodeLikelihood(from,-1);
	double newlogprior = LogBranchLengthPrior(from->GetBranch());
	double delta = newlogprior + newloglikelihood - currentlogprior - currentloglikelihood + loghastings;
	
	int accepted = (log(rnd::GetRandom().Uniform()) < delta);
	if (!accepted)	{
		RestoreBranch(from->GetBranch());
		Propagate(condlmap[0],GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()),condflag);
		ComputeNodeLikelihood(from,-1);
		// not useful: done by ComputeNodeLikelihood(from) just above
		// logL = currentloglikelihood;
	}
	return (double) accepted;
}

double PhyloProcess::GlobalProposeMove(const Branch* branch, double tuning)	{

	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);

	if (GetNprocs() > 1)	{
		// MPI
		// master and all slaves should all call MoveBranch(branch,m)
		// should send a message with arguments: GetBranchIndex(branch), m
		// slaves should interpret the message, and apply on branch with index received as message argument
		MESSAGE signal = PROPOSE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		prop_arg args;
		args.time = m;
		args.condalloc = branch->GetIndex();
		MPI_Bcast(&args,1,Propagate_arg,0,MPI_COMM_WORLD);
	}

	MoveBranch(branch,m);
	return m;
}

void PhyloProcess::SlaveProposeMove(int n,double x) {
	const Branch* br = GetBranch(n);
	MoveBranch(br,x);	
}

void PhyloProcess::GlobalRestoreBranch(const Branch* branch)	{

	if (GetNprocs() > 1)	{
		// MPI
		// master and all slaves should all call RestoreBranch(branch)
		MESSAGE signal = RESTOREBRANCH;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		int n = branch->GetIndex();
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	RestoreBranch(branch);
}

void PhyloProcess::SlaveRestoreBranch(int n) {
	const Branch* br = GetBranch(n);
	RestoreBranch(br);
}

//-------------------------------------------------------------------------
//	* TopoMoves 
//-------------------------------------------------------------------------

double PhyloProcess::SimpleTopoMoveCycle(int nrep, double tuning)	{

	for (int rep=0; rep<nrep; rep++)	{

		/*
		double x = rnd::GetRandom().Uniform();
		if (x < 0.1)	{
			BranchLengthMove(tuning);
		}
		else if (x < 0.2)	{
			BranchLengthMove(0.1 * tuning);
		}
		*/

		if (! fixtopo)	{
			MoveTopo();
		}

		/*
		double y = rnd::GetRandom().Uniform();
		if (y < 0.1)	{
			BranchLengthMove(tuning);
		}
		else if (y < 0.2)	{
			BranchLengthMove(0.1 * tuning);
		}
		*/
	}
}

double PhyloProcess::SPRMove(int nrep)	{
	sprchrono.Start();
	double tmp = GibbsSPR(nrep,0);
	spracc += tmp;
	sprtry ++;
	sprchrono.Stop();
}

double PhyloProcess::NNIMove(int nrep, double tuning)	{

	nnichrono.Start();
	for(int i=0; i<nrep; i++){
		double tmp = GibbsNNI(tuning,1);
		nniacc += tmp;
		nnitry ++;
	}
	nnichrono.Stop();
}

double PhyloProcess::MoveTopo()	{

	if (size < topoburnin)	{
		return 0;
	}

	// all moves reroot the tree and make their own likelihood updates before starting
	// but they don't leave with likelihoods updated
	double success = 0;

	if (NNNI)	{
		nnichrono.Start();
		if (GetNprocs() == 1)	{
			cerr << "error in PhyloProcess::MoveTopo: NNI not implemented in non-mpi mode\n";
			exit(1);
		}
		for(int i=0; i<NNNI; i++){
			double tmp = GibbsNNI(0.1,1);
			success += tmp;
			nniacc += tmp;
			nnitry ++;
		}
		nnichrono.Stop();
	}

	if (NSPR)	{
		sprchrono.Start();
		double tmp = GibbsSPR(NSPR,0);
		sprchrono.Stop();
		success += tmp;
		spracc += tmp;
		sprtry ++;
	}
	if (NMHSPR)	{
		sprchrono.Start();
		double tmp = GibbsMHSPR(topolambda,NMHSPR,0);
		sprchrono.Stop();
		success += tmp;
		mhspracc += tmp;
		mhsprtry ++;
	}
	if (NTSPR)	{
        for (int rep=0; rep<NTSPR; rep++)   {
            tsprchrono.Start();
            double logbf = 0;
            double deltalogp = 0;
            double tmp = TemperedGibbsSPR(topolambda,topomu,toponstep,1,0,deltalogp,logbf);
            // double tmp = TemperedGibbsSPR(topolambda,topomu,toponstep,1,2,deltalogp,logbf);
            tsprchrono.Stop();
            success += tmp;
            tspracc += tmp;
            tsprtry ++;
        }
	}
	if (nspec)	{
		sprchrono.Start();
		double tmp = GibbsSPR(nspec,1);
		//double tmp = GibbsMHSPR(topolambda,nspec,1);
		sprchrono.Stop();
		success += tmp;
		spracc += tmp;
		sprtry ++;
	}

	if (ntspec)	{
		tsprchrono.Start();
		double logbf;
		double deltalogp = 0;
		double tmp = TemperedGibbsSPR(0,topomu,toponstep,ntspec,1,deltalogp,logbf);
		tsprchrono.Stop();
		success += tmp;
		tspracc += tmp;
		tsprtry ++;
	}

	if (nbpp)	{
		sprchrono.Start();
		double tmp = 0;
		if (bpp == 3)	{
			tmp = GibbsMHSPR(bppbeta,nbpp,0);
		}
		else	{
			tmp = BPPSPR(nbpp);
		}
		sprchrono.Stop();
		success += tmp;
		bppspracc += tmp;
		bppsprtry ++;
	}

	if (ntbpp)	{
		tsprchrono.Start();
		double tmp = 0;
		if (bpp == 3)	{
			double logbf = 0;
			double deltalogp = 0;
			tmp = TemperedGibbsSPR(bppbeta,topomu,toponstep,ntbpp,0,deltalogp,logbf);
		}
		else	{
			tmp = TemperedBPPSPR(ntbpp,bppnstep);
		}
		tsprchrono.Stop();
		success += tmp;
		tbppspracc += tmp;
		tbppsprtry ++;
	}

	GlobalUpdateConditionalLikelihoods();

	return success;
}

void PhyloProcess::GlobalRootAt(Link* newroot)	{

	int n = GetLinkIndex(newroot);

	if (GetNprocs() > 1)	{
		// MPI
		// call slaves, send a reroot message with argument newroot
		MESSAGE signal = ROOTAT;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	}

	GetTree()->RootAt(newroot);
	// this may not always be necessary
	// should perhaps delegate this update to calling functions
	GlobalUpdateConditionalLikelihoods();	
}


void PhyloProcess::SlaveRootAt(int n) {
	Link* newroot = GetLink(n);
	GetTree()->RootAt(newroot);
}

void PhyloProcess::GlobalRootAtRandom()	{

	int n = GetTree()->CountInternalNodes(GetRoot());
	int choose = (int) (n * rnd::GetRandom().Uniform());

	if (GetNprocs() > 1)	{
		// MPI
		// call slaves, send a reroot message with argument newroot
		MESSAGE signal = ROOT;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&choose,1,MPI_INT,0,MPI_COMM_WORLD);
	}

	Link* tmp = 0;
	Link* newroot = GetTree()->ChooseInternalNode(GetRoot(),tmp,choose);
	if (newroot->isLeaf())	{
		cerr << "error : root at leaf\n";
		exit(1);
	}
	GetTree()->RootAt(newroot);
	// this may not always be necessary
	// should perhaps delegate this update to calling functions
	GlobalUpdateConditionalLikelihoods();	
}

void PhyloProcess::SlaveRoot(int n) {
	Link* tmp = 0;
	Link* newroot = GetTree()->ChooseInternalNode(GetRoot(),tmp,n);
	GetTree()->RootAt(newroot);
}

/*
void PhyloProcess::GlobalResetTree()	{
	
	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = RESETTREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	ResetTree();
}

void PhyloProcess::SlaveResetTree()	{

	ResetTree();
}

void PhyloProcess::ResetTree()	{

	if (tree)	{
		delete tree;
		tree = 0;
	}
	tree = new Tree(treefile);
	tree->RegisterWith(GetData()->GetTaxonSet());
}
*/

void PhyloProcess::GlobalBackupTree()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = BACKUPTREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	GetTree()->Backup();
	Backup();
}

void PhyloProcess::SlaveBackupTree()	{
	GetTree()->Backup();
	Backup();
}

void PhyloProcess::GlobalRestoreTree()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = RESTORETREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	GetTree()->Restore();
	Restore();
	GlobalUpdateParameters();
}

void PhyloProcess::SlaveRestoreTree()	{
	GetTree()->Restore();
	Restore();
}

void PhyloProcess::GlobalSwapTree()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = SWAPTREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	GetTree()->Swap();
}

void PhyloProcess::SlaveSwapTree()	{

	GetTree()->Swap();
}

void PhyloProcess::SetTestSiteMinAndMax()	{

	if (GetMyid() > 0) {
		int width = testnsite/(GetNprocs()-1);
		testsitemin = (GetMyid()-1)*width;
		testsitemax = 0;
		if (GetMyid() == (GetNprocs()-1)) {
			testsitemax = testnsite;
		}
		else {
			testsitemax = GetMyid()*width;
		} 
	}
}

void PhyloProcess::WaitLoop()	{
	MESSAGE signal;
	do {
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		if (signal == KILL) break;
		SlaveExecute(signal);
	} while(true);
}

void PhyloProcess::SlaveExecute(MESSAGE signal)	{
	int n,arg[4];
	prop_arg alpha;
	bool tvalue;

	switch(signal) {

	case SETRATEPRIOR:
		SlaveSetRatePrior();
		break;
	case SETPROFILEPRIOR:
		SlaveSetProfilePrior();
		break;
	case SETROOTPRIOR:
		SlaveSetRootPrior();
		break;
	case RESHUFFLE:
		SlaveReshuffleSites();
		break;
	case SUMRATE:
		SlaveActivateSumOverRateAllocations();
		break;
	case CONDRATE:
		SlaveInactivateSumOverRateAllocations();
		break;
	case SMCADDSITES:
		SMCAddSites();
		break;
	case RESETNSITE:
		ResetNsite();
		break;
	case SETNSITE:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		IncrementNsite(n);
		break;
	case ROOT:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveRoot(n);
		break;
	case ROOTAT:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveRootAt(n);
		break;
	case TOPOLRT:
		SlaveComputeTopoBFLogLikelihoodRatio();
		break;
	case SETBFFRAC:
		SlaveSetBFFrac();
		break;
	/*
	case RESCALEBRANCH:
		SlaveRescaleBranchPrior();
		break;
	case SETBRANCHALLOC:
		SlaveSetBranchAlloc();
		break;
	*/
	case SETSISFRAC:
		SlaveSetSISFrac();
		break;
	case COLLECTLIKELIHOOD:
		SlaveCollectLogLikelihood();
		break;
	case LIKELIHOOD:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveComputeNodeLikelihood(arg[0],arg[1]);
		break;
	case FULLLIKELIHOOD:
		SlaveGetFullLogLikelihood();
		break;
	case SCAN:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveGibbsSPRScan(arg[0],arg[1]);
		break;
	case PROPOSE:
		MPI_Bcast(&alpha,1,Propagate_arg,0,MPI_COMM_WORLD);
		SlaveProposeMove(alpha.condalloc,alpha.time);
		break;
	/*
	case RESETTREE:
		SlaveResetTree();
		break;
	*/
	case BACKUPTREE:
		SlaveBackupTree();
		break;
	case RESTORETREE:
		SlaveRestoreTree();
		break;
	case SWAPTREE:
		SlaveSwapTree();
		break;
	case RESTOREBRANCH:
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveRestoreBranch(n);
		break;
	case RESET:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[1] == 1) ? true : false;
		SlaveReset(arg[0],tvalue);
		break;
	case MULTIPLY:
		MPI_Bcast(arg,3,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[2] == 1) ? true : false;
		SlaveMultiply(arg[0],arg[1],tvalue);
		break;
	case SMULTIPLY:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[1] == 1) ? true : false;
		SlaveMultiplyByStationaries(arg[0],tvalue);
		break;
	case INITIALIZE:
		MPI_Bcast(arg,3,MPI_INT,0,MPI_COMM_WORLD);
		tvalue = (arg[2] == 1) ? true : false;
		SlaveInitialize(arg[0],arg[1],tvalue);
		break;
	case PROPAGATE:
		MPI_Bcast(&alpha,1,Propagate_arg,0,MPI_COMM_WORLD);
		tvalue = (alpha.condalloc == 1) ? true : false;
		SlavePropagate(alpha.from,alpha.to,tvalue,alpha.time);
		break;
	case ATTACH:
		MPI_Bcast(arg,4,MPI_INT,0,MPI_COMM_WORLD);
		SlaveAttach(arg[0],arg[1],arg[2],arg[3]);
		break;
	case DETACH:
		MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);
		SlaveDetach(arg[0],arg[1]);
		break;
	case MINMAX:
		SlaveSetMinMax();
		break;
	case NNI:
		SlaveNNI();
		break;
	case KNIT:
		SlaveKnit();
		break;
	case BRANCHPROPAGATE:
		MPI_Bcast(arg,1,MPI_INT,0,MPI_COMM_WORLD);
		SlavePropagateOverABranch(arg[0]);
		break;
	case UNFOLD:
		SlaveUnfold();
		break;
	case COLLAPSE:
		SlaveCollapse();
		break;
	case 	UPDATE:
		SlaveUpdateConditionalLikelihoods();
		break;
	case UPDATE_SRATE:
		SlaveUpdateSiteRateSuffStat();
		break;
	case UPDATE_SPROFILE:
		SlaveUpdateSiteProfileSuffStat();
		break;
	case UPDATE_MPROFILE:
		SlaveUpdateModeProfileSuffStat();
		break;
	case UPDATE_BLENGTH:
		SlaveUpdateBranchLengthSuffStat();
		break;
	case PARAMETER_DIFFUSION:
		SlaveUpdateParameters();
		break;
	case BCAST_TREE:
		SlaveBroadcastTree();
		break;
	case UNCLAMP:
		SlaveUnclamp();
		break;
	case RESTOREDATA:
		SlaveRestoreData();
		break;
	case SETDATA:
		SlaveSetDataFromLeaves();
		break;
	case SETNODESTATES:
		SlaveSetNodeStates();
		break;
	case GETDIV:
		SlaveGetMeanDiversity();
		break;
	case CVSCORE:
		SlaveComputeCVScore();
		break;
	case SITELOGL:
		SlaveComputeSiteLogL();
		break;
	case SITERATE:
		SlaveSendMeanSiteRate();
		break;
	case SETTESTDATA:
		SlaveSetTestData();
		break;
	case WRITE_MAPPING:
		SlaveWriteMappings();
		break;
	case COUNTMAPPING:
		SlaveCountMapping();
		break;
	case ACTIVATEZIP:
		SlaveActivateZip();
		break;
	case INACTIVATEZIP:
		SlaveInactivateZip();
		break;
	
	default:
		cerr << "slave could not process signal : " << signal << '\n';
		exit(1);
	}
}

//-----------------------------PhyloProcess--------------------------------------------
//-------------------------------------------------------------------------
//	* MPI SuffStat 
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PhyloProcess::GlobalUpdateBranchLengthSuffStat()	{

	if (GetNprocs() > 1)	{
	MPI_Status stat;
	MESSAGE signal = UPDATE_BLENGTH;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(int i=0; i<GetNbranch(); i++) {
		branchlengthsuffstatcount[i] = 0;
		branchlengthsuffstatbeta[i] = 0.0;
	}

	double ivector[GetNbranch()];
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(ivector,GetNbranch(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNbranch(); j++) {
			branchlengthsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	double dvector[GetNbranch()];
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(dvector,GetNbranch(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNbranch(); ++j) {
			branchlengthsuffstatbeta[j] += dvector[j];
		}
	}

	if (branchlengthsuffstatcount[0])	{
		cerr << "error at root\n";
		cerr << branchlengthsuffstatcount[0] << '\n';
	}
	if (branchlengthsuffstatbeta[0])	{
		cerr << "error at root\n";
		cerr << branchlengthsuffstatbeta[0] << '\n';
	}
	}
	else	{
		UpdateBranchLengthSuffStat();
	}

	// check for nan
	for(int j=0; j<GetNbranch(); ++j) {
		if (isnan(branchlengthsuffstatbeta[j]))	{
			cerr << "in PhyloProcess::GlobalUpdateBranchLengthSuffStat: nan\n";
			exit(1);
		}
	}
}

void PhyloProcess::SlaveUpdateBranchLengthSuffStat()	{

	UpdateBranchLengthSuffStat();
	if (branchlengthsuffstatcount[0])	{
		cerr << "error at root in slave " << GetMyid() << "\n";
		cerr << branchlengthsuffstatcount[0] << '\n';
	}
	if (branchlengthsuffstatbeta[0])	{
		cerr << "error at root in slave " << GetMyid() << "\n";
		cerr << branchlengthsuffstatbeta[0] << '\n';
	}
	MPI_Send(branchlengthsuffstatcount,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(branchlengthsuffstatbeta,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalUpdateSiteRateSuffStat()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = UPDATE_SRATE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateSiteRateSuffStat();
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (isnan(siteratesuffstatbeta[i]))	{
				cerr << "in PhyloProcess::GlobalUpdateSiteRateSuffStat: nan ratesuffstatbeta\n";
				exit(1);
			}
		}
	}
}

void PhyloProcess::SlaveUpdateSiteRateSuffStat()	{

	UpdateSiteRateSuffStat();
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (isnan(siteratesuffstatbeta[i]))	{
			cerr << "in PhyloProcess::GlobalUpdateSiteRateSuffStat: nan ratesuffstatbeta\n";
			exit(1);
		}
	}
}


