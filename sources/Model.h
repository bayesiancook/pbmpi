
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "CodonSequenceAlignment.h"
#include "RASCATGTRFiniteGammaPhyloProcess.h"
#include "RASCATGTRSBDPGammaPhyloProcess.h"
#include "GeneralPathSuffStatRASCATGTRDPGammaPhyloProcess.h"
#include "GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.h"
#include "RASCATGTRDPGammaPhyloProcess.h"
#include "RASCATGammaPhyloProcess.h"
#include "RASCATFiniteGammaPhyloProcess.h"
#include "RASCATSBDPGammaPhyloProcess.h"
#include "GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.h"
#include "AAMutSelFinitePhyloProcess.h"
#include "AACodonMutSelFinitePhyloProcess.h"
#include "AAMutSelSBDPPhyloProcess.h"
#include "AAMutSelSiteSpecificPhyloProcess.h"
#include "AAMutSelDPPhyloProcess.h"
#include "CodonMutSelFinitePhyloProcess.h"
#include "CodonMutSelSBDPPhyloProcess.h"
#include "AACodonMutSelSBDPPhyloProcess.h"
#include "PartitionedRASGTRGammaPhyloProcess.h"
#include "PartitionedRASCATGTRSBDPGammaPhyloProcess.h"
#include "PartitionedRASCATGTRFiniteGammaPhyloProcess.h"
#include "Parallel.h"
#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

MPI_Datatype Propagate_arg;



class Model	{

	public:

	PhyloProcess* process;
	string type;
	string name;
	int every;
	int until;
	int saveall;
	int incinit;

	int w;
	int stone_index;
	int num_stones;
	stringstream ss_ext;

	Model(string datafile, string treefile, string schemefile, bool inlinkgam, bool inunlinkgtr, bool inlinkmult,  int modeltype, int nratecat, int mixturetype, int ncat, GeneticCodeType codetype, int suffstat, int fixncomp, int empmix, string mixtype, string rrtype, int iscodon, int fixtopo, int NSPR, int NNNI, int fixcodonprofile, int fixomega, int fixbl, int omegaprior, int kappaprior, int dirweightprior, double mintotweight, int dc, int inevery, int inuntil, int insaveall, int inincinit, string inname, int myid, int nprocs)	{

		every = inevery;
		until = inuntil;
		name = inname;
		saveall = insaveall;
		incinit = inincinit;
		stone_index = -1;
		num_stones = 0;

		// 1 : CAT
		// 2 : CATGTR
		// 3 : MutSel

		// mixturetype
		// 0 : one for all
		// 1 : finite
		// 2 : dp
		// 3 : tdp : removed
		// 4 : sbdp
		
		// PARTITIONED
		if (schemefile != "None")
		{
			string inrrtype = (modeltype == 2) ? rrtype : "Part";
			if (mixturetype == 1)	{
				type = "PARTCATFINITE";
				process = new PartitionedRASCATGTRFiniteGammaPhyloProcess(datafile,treefile,schemefile,inlinkgam,inunlinkgtr,inlinkmult,inrrtype,nratecat,ncat,fixncomp,empmix,mixtype,fixtopo,NSPR,NNNI,myid,nprocs);
			}
			else if (mixturetype == 3)	{
				type = "PARTCATSBDP";
				process = new PartitionedRASCATGTRSBDPGammaPhyloProcess(datafile,treefile,schemefile,inlinkgam,inunlinkgtr,inlinkmult,inrrtype,nratecat,iscodon,codetype,fixtopo,NSPR,NNNI,kappaprior,mintotweight,myid,nprocs);
			}
			else
			{
				type = "PARTCATFIX";
				process = new PartitionedRASGTRGammaPhyloProcess(datafile,treefile,schemefile,inlinkgam,inunlinkgtr,inlinkmult,inrrtype,nratecat,iscodon,codetype,fixtopo,NSPR,NNNI,mintotweight,myid,nprocs);
			}
		}

		// CAT
		else if (modeltype == 1)	{
			if (myid == 0) {
				// cerr << "cat model\n";
			}
			if (mixturetype == 1)	{
				type = "CATFINITE";
				process = new RASCATFiniteGammaPhyloProcess(datafile,treefile,nratecat,ncat,fixncomp,empmix,mixtype,fixtopo,NSPR,NNNI,dc,myid,nprocs); 
			}
			else	{
				type = "CATSBDP";
				process = new RASCATSBDPGammaPhyloProcess(datafile,treefile,nratecat,iscodon,codetype,fixtopo,NSPR,NNNI,kappaprior,mintotweight,dc,incinit,myid,nprocs); 
			}
		}

		// CATGTR
		else if (modeltype == 2)	{
			if (myid == 0) {
				// cerr << "catgtr model\n";
			}
			if (mixturetype == 1)	{
				if (suffstat)	{
					type = "CATGTRFINITE";
					process = new RASCATGTRFiniteGammaPhyloProcess(datafile,treefile,nratecat,ncat,fixncomp,empmix,mixtype,rrtype,fixtopo,NSPR,NNNI,dc,myid,nprocs); 
				}
				else	{
					cerr << "gpss deprecated\n";
					exit(1);
				}
			}
			else if (mixturetype == 2)	{
				cerr << "simple dp deprecated\n";
				exit(1);
			}
			else if (mixturetype == 3)	{
				if (suffstat)	{
					type = "CATGTRSBDP";
					process = new RASCATGTRSBDPGammaPhyloProcess(datafile,treefile,nratecat,iscodon,codetype,rrtype,fixtopo,NSPR,NNNI,kappaprior,mintotweight,dc,incinit,myid,nprocs); 
				}
				else	{
					cerr << "gpss deprecated\n";
					exit(1);
				}
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized\n";
				exit(1);
			}
		}

		// AAMUTSEL
		else if (modeltype == 3)	{
			cerr << "deprecated.\n";
			exit(1);
			if (mixturetype == 1)	{
				type = "AAMUTSELFINITE";
				process = new AAMutSelFinitePhyloProcess(datafile,treefile,codetype,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,dc,myid,nprocs);
			}
			else if (mixturetype == 2)	{
				type = "AAMUTSELDP";
				process = new AAMutSelDPPhyloProcess(datafile,treefile,codetype,fixtopo,fixbl,kappaprior,dc,myid,nprocs);
			}
			else if (mixturetype == 3)	{
				type = "AAMUTSELSBDP";
				process = new AAMutSelSBDPPhyloProcess(datafile,treefile,codetype,fixtopo,fixbl,kappaprior,dc,myid,nprocs);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized\n";
				exit(1);
			}
		}
		

		// CodonMutSel
		else if (modeltype == 4)	{
			if (mixturetype == 1)	{
				type = "CODONMUTSELFINITE";
				process = new CodonMutSelFinitePhyloProcess(datafile,treefile,codetype,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,NSPR,NNNI,dirweightprior,dc,myid,nprocs);
			}
			else if (mixturetype == 3)	{
				type = "CODONMUTSELSBDP";
				process = new CodonMutSelSBDPPhyloProcess(datafile,treefile,codetype,fixtopo,fixbl,NSPR,NNNI,kappaprior,mintotweight,dc,myid,nprocs);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized or not yet implemented.\n";
				exit(1);
			}
		}

		// AACodonMutSel
		else	{
			if (mixturetype == 1)	{
				type = "AACODONMUTSELFINITE";
				process = new AACodonMutSelFinitePhyloProcess(datafile,treefile,codetype,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,NSPR,NNNI,fixcodonprofile,fixomega,omegaprior,dirweightprior,dc,myid,nprocs);
			}
			else if (mixturetype == 3)	{
				type = "AACODONMUTSELSBDP";
				process = new AACodonMutSelSBDPPhyloProcess(datafile,treefile,codetype,fixtopo,fixbl,NSPR,NNNI,fixcodonprofile,fixomega,omegaprior,kappaprior,dirweightprior,mintotweight,dc,myid,nprocs);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized or not yet implemented.\n";
				exit(1);
			}
			
		}

	}

	Model(string inname, int myid, int nprocs, string treefile = "None", int innum_stones = 0, int instone_size = -1) : process(NULL)	{

	    name = inname;

        ifstream is((name + ".param").c_str());
        if (! is)   {
            cerr << "error: cannot open " << name << ".param\n";
            exit(1);
        }

	    int size;

	    is >> type;
        is >> every >> until >> size;
        is >> saveall;
        is >> stone_index >> num_stones;

	    // if we are elongating multiple stones
        // begin with the first stone
        if(innum_stones > 0)
        {
            num_stones = innum_stones;
            until = instone_size;
            stone_index = 0;
            saveall = false;
            w = num_stones > 0 ? (int) log10 ((double) num_stones) + 1 : 1;

            ss_ext.str("");
            ss_ext << "_ss" << setfill('0') << setw(w) << stone_index;
        }
        // otherwise we elongate a single chain
        else
        {
            // if the chain is a steppingstone
            // then continue forever
            // note that this can be overridden by -x
            if(num_stones > 0)
                until = -1;
        }
		
        if (type == "PARTCATFINITE")    {
            process = new PartitionedRASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "PARTCATSBDP") {
            process = new PartitionedRASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "PARTCATFIX")  {
            process = new PartitionedRASGTRGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CATDP")   {
            process = new RASCATGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CATSBDP") {
            process = new RASCATSBDPGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CATFINITE")   {
            process = new RASCATFiniteGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CATGTRDP")    {
            process = new RASCATGTRDPGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CATGTRSBDP")  {
            process = new RASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "GPSSCATGTRSBDP")  {
            process = new GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CATGTRFINITE")    {
            process = new RASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "GPSSCATGTRFINITE")    {
            process = new GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs);
        }
        else if (type == "AAMUTSELFINITE")  {
            process = new AAMutSelFinitePhyloProcess(is,myid,nprocs);
        }
        else if (type == "AACODONMUTSELFINITE") {
            process = new AACodonMutSelFinitePhyloProcess(is,myid,nprocs);
        }
        else if (type == "AAMUTSELSBDP")    {
            process = new AAMutSelSBDPPhyloProcess(is,myid,nprocs);
        }
        else if (type == "AAMUTSELDP")  {
            process = new AAMutSelDPPhyloProcess(is,myid,nprocs);
        }
        else if (type == "CODONMUTSELFINITE")   {
            process = new CodonMutSelFinitePhyloProcess(is,myid,nprocs);
        }
        else if (type == "CODONMUTSELSBDP") {
            process = new CodonMutSelSBDPPhyloProcess(is,myid,nprocs);
        }
        else if (type == "AACODONMUTSELSBDP")   {
            process = new AACodonMutSelSBDPPhyloProcess(is,myid,nprocs);
        }
        else    {
            cerr << "error, does not recognize model type : " << type << '\n';
            exit(1);
        }

	process->SetSize(size);

        if(num_stones > 0)
        {
            // allow use of a fixed topology with steppingstones
            if(treefile != "None")
            {
                process->FixTopo(treefile);
            }
            
            // set stepping stone if elongating a single stone
            if(myid == 0 && innum_stones == 0)
            {
                process->GlobalSetSteppingStone(stone_index, num_stones);
                process->GlobalCollapse();
                process->GlobalUnfold();
            }
        }
	}

	void ToStream(ostream& os, bool header)	{
		if (header)	{
			os << type << '\n';
			os << every << '\t' << until << '\t' << GetSize() << '\n';
			os << saveall << '\n';
			os << stone_index << '\t' << num_stones << '\n';
			process->ToStreamHeader(os);
		}
		process->ToStream(os);
	}

	~Model()	{
		// things here !!
	}

	void WaitLoop()	{
		process->WaitLoop();
	}

	double Move(double tuning, int nrep)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += process->Move(tuning);
		}
		return total / nrep;
	}

	int RunningStatus()	{
		ifstream ris((name + ss_ext.str() + ".run").c_str());
		int i;
		ris >> i;
		ris.close();
		return i;
	}

	int GetSize()	{
		return process->GetSize();
	}

	void RunStone(int burnin)
	{
	    cerr << "run " << (name + ss_ext.str()) << " started\n";
	    if (burnin != 0)    {
            if (GetSize() < burnin) {
                process->SetBurnin(true);
            }
        }

        while (((until == -1) || (GetSize() < until)))  {
            if (GetSize() >= burnin)    {
                process->SetBurnin(false);
            }

            Move(1,every);

            process->IncSize();

            ofstream os((name + ss_ext.str() + ".treelist").c_str(), ios_base::app);
            process->SetNamesFromLengths();
            process->RenormalizeBranchLengths();
            GetTree()->ToStream(os);
            process->DenormalizeBranchLengths();
            os.close();

            ofstream tos((name + ss_ext.str() + ".trace").c_str(), ios_base::app);
            Trace(tos);
            tos.close();

            if(num_stones > 0)
            {
                ofstream tos((name + ss_ext.str() + ".ss").c_str(), ios_base::app);
                tos << GetSize() - 1 << "\t" << process->GetMaskedLogLikelihood() << "\n";
                tos.close();
            }

            ofstream pos((name + ss_ext.str() + ".param").c_str());
            pos.precision(numeric_limits<double>::digits10);
            ToStream(pos,true);
            pos.close();

            if (saveall)    {
                ofstream cos((name + ss_ext.str() + ".chain").c_str(),ios_base::app);
                cos.precision(12);
                ToStream(cos,false);
                cos.close();
            }

        }
        cerr << (name + ss_ext.str()) << ": completed after " << GetSize() << " samples.\n\n";
	}

	void Run(int burnin)	{

	    while(stone_index < num_stones)
        {
	        // if we are elongating multiple stones
            if(!ss_ext.str().empty())
            {
                // check if the next stone already exists
                ifstream is((name + ss_ext.str() + ".param").c_str());
                if (! is )
                {
                    // if not, reset the process and create new output files
                    process->SetSize(0);
                    process->GlobalSetSteppingStone(stone_index, num_stones);
                    process->GlobalCollapse();
                    process->GlobalUnfold();

                    ofstream os((name + ss_ext.str() + ".treelist").c_str());
                    ofstream tos((name + ss_ext.str() + ".trace").c_str());
                    TraceHeader(tos);
                    tos.close();
                    ofstream ssos((name + ss_ext.str() + ".ss").c_str());
                    ssos << "iter\tss\n";
                    ssos.close();
                    ofstream pos((name + ss_ext.str() + ".param").c_str());
                    pos.precision(numeric_limits<double>::digits10);
                    ToStream(pos,true);
                    pos.close();
                }
                else
                {
                    // if so, then restore the process
                    int size;
                    is >> type;
                    int unt;
                    is >> every >> unt >> size;
                    is >> saveall;
                    is >> stone_index >> num_stones;

                    process->FromStreamHeader(is);
                    process->FromStream(is);
                    MESSAGE signal = BCAST_TREE;
		    MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		    process->GlobalBroadcastTree();
		    process->SetSize(size);
                    process->GlobalSetSteppingStone(stone_index, num_stones);
                    process->GlobalCollapse();
                    process->GlobalUnfold();

                    // if the stone has enough states, then continue to the next stone
                    if(size >= until)
                    {
                        cerr << (name + ss_ext.str()) << " already complete. continuing...\n\n";

                        is.close();

                        stone_index++;

                        ss_ext.str("");
                        ss_ext << "_ss" << setfill('0') << setw(w) << stone_index;
                        continue;
                    }
                }
            }

            RunStone(burnin);

            // if we were simply elongating a single stone,
            // then we can quit here
            if(ss_ext.str().empty())
                break;

            stone_index++;

            ss_ext.str("");
            ss_ext << "_ss" << setfill('0') << setw(w) << stone_index;
        }
	}

	NewickTree* GetTree() {return process->GetLengthTree();}

	void TraceHeader(ostream& os)	{
		process->TraceHeader(os);
	}

	void Trace(ostream& os)	{
		process->Trace(os);
	}

	void ReadPB(int argc, char* argv[])	{
		process->ReadPB(argc,argv);
	}
};
