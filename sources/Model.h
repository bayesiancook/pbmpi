
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
#include "RASCATFiniteGammaPhyloProcess.h"
#include "RASCATSBDPGammaPhyloProcess.h"
#include "AACodonMutSelSBDPPhyloProcess.h"
#include "AACodonMutSelFinitePhyloProcess.h"
#include "CodonMutSelFinitePhyloProcess.h"
#include "CodonMutSelSBDPPhyloProcess.h"

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
    int steppingstep;
    int steppingtaxstep;
    int steppingburnin;
    int steppingsize;

	Model(string datafile, string treefile, int modeltype, int nratecat, int mixturetype, int ncat, int nmodemax, GeneticCodeType codetype, int suffstat, int fixncomp, int empmix, string mixtype, string rrtype, int iscodon, int fixtopo, int NSPR, int NNNI, int fixcodonprofile, int fixomega, int fixbl, int omegaprior, int kappaprior, int dirweightprior, double mintotweight, int dc, int inevery, int inuntil, int insaveall, int inincinit, int topoburnin, int insteppingstep, int insteppingtaxstep, int insteppingburnin, int insteppingsize, string inname, int myid, int nprocs)	{

		every = inevery;
		until = inuntil;
		name = inname;
		saveall = insaveall;
		incinit = inincinit;
        steppingstep = insteppingstep;
        steppingtaxstep = insteppingtaxstep;
        steppingburnin = insteppingburnin;
        steppingsize = insteppingsize;

		// 1 : CAT
		// 2 : CATGTR
		// 3 : MutSel

		// mixturetype
		// 0 : one for all
		// 1 : finite
		// 2 : dp
		// 3 : tdp : removed
		// 4 : sbdp
		// 5 : site specific
		
		// CAT
		if (modeltype == 1)	{
			if (myid == 0) {
				// cerr << "cat model\n";
			}
			if (mixturetype == 1)	{
				type = "CATFINITE";
				process = new RASCATFiniteGammaPhyloProcess(datafile,treefile,nratecat,nmodemax,ncat,fixncomp,empmix,mixtype,dirweightprior,fixtopo,NSPR,NNNI,dc,myid,nprocs); 
			}
			else	{
				type = "CATSBDP";
				process = new RASCATSBDPGammaPhyloProcess(datafile,treefile,nratecat,nmodemax,iscodon,codetype,fixtopo,NSPR,NNNI,kappaprior,dirweightprior,mintotweight,dc,incinit,myid,nprocs); 
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
					process = new RASCATGTRFiniteGammaPhyloProcess(datafile,treefile,nratecat,nmodemax,ncat,fixncomp,empmix,mixtype,rrtype,dirweightprior,fixtopo,NSPR,NNNI,dc,myid,nprocs); 
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
					process = new RASCATGTRSBDPGammaPhyloProcess(datafile,treefile,nratecat,nmodemax,iscodon,codetype,rrtype,fixtopo,NSPR,NNNI,kappaprior,dirweightprior,mintotweight,dc,incinit,myid,nprocs); 
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
		}
		

		// CodonMutSel
		else if (modeltype == 4)	{
			if (mixturetype == 1)	{
				type = "CODONMUTSELFINITE";
				process = new CodonMutSelFinitePhyloProcess(datafile,treefile,codetype,nmodemax,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,NSPR,NNNI,dirweightprior,dc,myid,nprocs);
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
				process = new AACodonMutSelFinitePhyloProcess(datafile,treefile,codetype,nmodemax,ncat,fixncomp,empmix,mixtype,fixtopo,fixbl,NSPR,NNNI,fixcodonprofile,fixomega,omegaprior,dirweightprior,dc,myid,nprocs);
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

		process->SetFixBL(fixbl);
		if (fixbl && (! myid))	{
			cerr << "set lengths from tree file\n";
			process->SetLengthsFromNames();
		}
		process->SetTopoBurnin(topoburnin);
	}

	Model(string inname, int myid, int nprocs)	{

		name = inname;

		ifstream is((name + ".param").c_str());
		if (! is)	{
			cerr << "error: cannot open " << name << ".param\n";
			exit(1);
		}

		is >> type;
        if (type == "STEPPING") {
            is >> steppingstep >> steppingtaxstep >> steppingburnin >> steppingsize;
            is >> type;
        }
		int size;
		is >> every >> until >> size;
		is >> saveall;
		
		if (type == "CATSBDP")	{
			process = new RASCATSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATFINITE")	{
			process = new RASCATFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRSBDP")	{
			process = new RASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRFINITE")	{
			process = new RASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "AACODONMUTSELFINITE")	{
			process = new AACodonMutSelFinitePhyloProcess(is,myid,nprocs);
		}
		else if (type == "CODONMUTSELFINITE")	{
			process = new CodonMutSelFinitePhyloProcess(is,myid,nprocs);
		}
		else if (type == "CODONMUTSELSBDP")	{
			process = new CodonMutSelSBDPPhyloProcess(is,myid,nprocs);
		}
		else if (type == "AACODONMUTSELSBDP")	{
			process = new AACodonMutSelSBDPPhyloProcess(is,myid,nprocs);
		}
		else	{
			cerr << "error, does not recognize model type : " << type << '\n';
			exit(1);
		}

		process->SetSize(size);
	}

	void ToStream(ostream& os, bool header)	{
		stringstream ss;
		if (header)	{
            if (steppingstep)  {
                ss << "STEPPING\n";
                ss << steppingstep << '\t' << steppingtaxstep << '\t' << steppingburnin << '\t' << steppingsize << '\n';
            }
			ss << type << '\n';
			ss << every << '\t' << until << '\t' << GetSize() << '\n';
			ss << saveall << '\n';
			process->ToStreamHeader(ss);
		}
		process->ToStream(ss);
		os << ss.str();
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
		ifstream ris((name + ".run").c_str());
		int i;
		ris >> i;
		ris.close();
		return i;
	}

	int GetSize()	{
		return process->GetSize();
	}
	
	void IncSize()	{
		process->IncSize();
	}

    void Run(int burnin)    {
        if (! steppingstep)    {
            MCMCRun(burnin);
        }
        else    {
            SteppingRun(steppingstep, steppingtaxstep, steppingburnin, steppingsize);
        }
    }

	void MCMCRun(int burnin)	{
		if (burnin != 0)	{
			if (GetSize() < burnin)	{
				process->SetBurnin(true);
			}
		}
		ofstream ros((name + ".run").c_str()); stringstream buf;
		buf << 1 << '\n';
		ros << buf.str();
		ros.close();
	
		while (RunningStatus() && ((until == -1) || (GetSize() < until)))	{
			if (GetSize() >= burnin)	{
				process->SetBurnin(false);
			}

			Move(1,every);
			
			process->IncSize();

			ofstream os((name + ".treelist").c_str(), ios_base::app);
			TreeTrace(os);
			os.close();

			ofstream tos((name + ".trace").c_str(), ios_base::app);
			Trace(tos);
			tos.close();

			ofstream mos((name + ".monitor").c_str());
			Monitor(mos);
			mos.close();

			ofstream pos((name + ".param").c_str());
			pos.precision(numeric_limits<double>::digits10);
			ToStream(pos,true);
			pos.close();

			if (saveall)	{
				ofstream cos((name + ".chain").c_str(),ios_base::app);
				cos.precision(numeric_limits<double>::digits10);
				ToStream(cos,false);
				cos.close();
			}

		}	
		cerr << name << ": stopping after " << GetSize() << " points.\n";
		cerr << '\n';
	}

	void SteppingRun(int step, int taxstep, int burnin, int stepsize)	{

		process->GlobalPrepareStepping();
	
		// total size : nstep * (burnin + stepsize)
        // number of steps:
        //

        int sitencycle = process->GetNsite() / step;
        if (process->GetNsite() % step) {
            sitencycle++;
        }
        int taxncycle = process->GetNtaxa() / taxstep;
        if (process->GetNtaxa() % taxstep)  {
            taxncycle++;
        }

        int until = sitencycle * taxncycle * (burnin + stepsize);

		ofstream ros((name + ".run").c_str()); stringstream buf;
		buf << 1 << '\n';
		ros << buf.str();
		ros.close();
	
        if (! GetSize())    {
            process->PriorSample();
            process->GlobalUpdateParameters();
        }

		while (RunningStatus() && (GetSize() < until))	{

            int cycle = int(GetSize() / (burnin + stepsize));
            int sitecycle = cycle / taxncycle;
            int taxcycle = cycle % taxncycle;
            int nsite = sitecycle * step;
            if (nsite > process->GetNsite()) {
                nsite = process->GetNsite();
            }
            int ntaxa = taxcycle * taxstep;
            if (ntaxa >= process->GetNtaxa()) {
                cerr << "error: ntaxa > Ntaxa\n";
                exit(1);
            }
            int cutoff = process->GetNtaxa()*nsite + ntaxa;

            int taxcycle2 = taxcycle;
            int sitecycle2 = sitecycle;
            taxcycle2++;
            if (taxcycle2 == taxncycle)   {
                taxcycle2 = 0;
                sitecycle2++;
            }
            int nsite2 = sitecycle2 * step;
            if (nsite2 > process->GetNsite()) {
                nsite2 = process->GetNsite();
            }
            int ntaxa2 = taxcycle2 * taxstep;
            if (ntaxa2 >= process->GetNtaxa()) {
                cerr << "error: ntaxa2 > Ntaxa\n";
                exit(1);
            }
            int cutoff2 = process->GetNtaxa()*nsite2 + ntaxa2;

            /*
            cerr << sitecycle << '\t' << nsite << '\t' << taxcycle << '\t' << ntaxa << '\t' << cutoff << '\t';
            cerr << '\n';
            cerr << sitecycle2 << '\t' << nsite2 << '\t' << taxcycle2 << '\t' << ntaxa2 << '\t' << cutoff2 << '\n';
            cerr << '\n';
            cerr << '\n';
            */

			process->GlobalSetSteppingFraction(cutoff);

			Move(1,every);
			
			process->IncSize();

			ofstream os((name + ".treelist").c_str(), ios_base::app);
			TreeTrace(os);
			os.close();

			ofstream tos((name + ".trace").c_str(), ios_base::app);
			Trace(tos);
			tos.close();

			ofstream mos((name + ".monitor").c_str());
			Monitor(mos);
			mos.close();

			ofstream pos((name + ".param").c_str());
			pos.precision(numeric_limits<double>::digits10);
			ToStream(pos,true);
			pos.close();

			if (saveall)	{
				ofstream cos((name + ".chain").c_str(),ios_base::app);
				cos.precision(numeric_limits<double>::digits10);
				ToStream(cos,false);
				cos.close();
			}

			ofstream los((name + ".stepping").c_str(), ios_base::app);
			double lnL1 = process->GlobalGetFullLogLikelihood();
			process->GlobalSetSteppingFraction(cutoff2);
			double lnL2 = process->GlobalGetFullLogLikelihood();
            double dlnL = lnL2 - lnL1;
            if (std::isnan(dlnL))   {
                cerr << "nan lnl\n";
                cerr << lnL1 << '\t' << lnL2 << '\n';
                exit(1);
            }
            int dnsite = nsite2 - nsite;
            int dntaxa = ntaxa2 - ntaxa;
			los << cutoff << '\t' << dnsite << '\t' << dlnL << '\t' << dntaxa << '\n';
			// los << cutoff << '\t' << dnsite << '\t' << dlnL << '\t' << dlnL / dnsite << '\n';
			los.close();
		}	
		cerr << name << ": stopping after " << GetSize() << " points.\n";
		cerr << '\n';
	}

	void EmpiricalSteppingRun(string empname, int step, int taxstep, int burnin, int stepsize)	{

        ifstream is(empname.c_str());
        process->GlobalSetEmpiricalPrior(is);

		process->GlobalPrepareStepping();
	
		// total size : nstep * (burnin + stepsize)
        // number of steps:
        //

        int sitencycle = process->GetNsite() / step;
        if (process->GetNsite() % step) {
            sitencycle++;
        }
        int taxncycle = process->GetNtaxa() / taxstep;
        if (process->GetNtaxa() % taxstep)  {
            taxncycle++;
        }

        int until = sitencycle * taxncycle * (burnin + stepsize);

		ofstream ros((name + ".run").c_str()); stringstream buf;
		buf << 1 << '\n';
		ros << buf.str();
		ros.close();
	
        if (! GetSize())    {
            process->PriorSample();
            process->GlobalUpdateParameters();
        }

		while (RunningStatus() && (GetSize() < until))	{

            int cycle = int(GetSize() / (burnin + stepsize));
            int sitecycle = cycle / taxncycle;
            int taxcycle = cycle % taxncycle;
            int nsite = sitecycle * step;
            if (nsite > process->GetNsite()) {
                nsite = process->GetNsite();
            }
            int ntaxa = taxcycle * taxstep;
            if (ntaxa >= process->GetNtaxa()) {
                cerr << "error: ntaxa > Ntaxa\n";
                exit(1);
            }
            int cutoff = process->GetNtaxa()*nsite + ntaxa;

            int taxcycle2 = taxcycle;
            int sitecycle2 = sitecycle;
            taxcycle2++;
            if (taxcycle2 == taxncycle)   {
                taxcycle2 = 0;
                sitecycle2++;
            }
            int nsite2 = sitecycle2 * step;
            if (nsite2 > process->GetNsite()) {
                nsite2 = process->GetNsite();
            }
            int ntaxa2 = taxcycle2 * taxstep;
            if (ntaxa2 >= process->GetNtaxa()) {
                cerr << "error: ntaxa2 > Ntaxa\n";
                exit(1);
            }
            int cutoff2 = process->GetNtaxa()*nsite2 + ntaxa2;

            /*
            cerr << sitecycle << '\t' << nsite << '\t' << taxcycle << '\t' << ntaxa << '\t' << cutoff << '\t';
            cerr << '\n';
            cerr << sitecycle2 << '\t' << nsite2 << '\t' << taxcycle2 << '\t' << ntaxa2 << '\t' << cutoff2 << '\n';
            cerr << '\n';
            cerr << '\n';
            */

			process->GlobalSetSteppingFraction(cutoff);

			Move(1,every);
			
			process->IncSize();

			ofstream os((name + ".treelist").c_str(), ios_base::app);
			TreeTrace(os);
			os.close();

			ofstream tos((name + ".trace").c_str(), ios_base::app);
			Trace(tos);
			tos.close();

			ofstream mos((name + ".monitor").c_str());
			Monitor(mos);
			mos.close();

			ofstream pos((name + ".param").c_str());
			pos.precision(numeric_limits<double>::digits10);
			ToStream(pos,true);
			pos.close();

			if (saveall)	{
				ofstream cos((name + ".chain").c_str(),ios_base::app);
				cos.precision(numeric_limits<double>::digits10);
				ToStream(cos,false);
				cos.close();
			}

			ofstream los((name + ".stepping").c_str(), ios_base::app);
			double lnL1 = process->GlobalGetFullLogLikelihood();
			process->GlobalSetSteppingFraction(cutoff2);
			double lnL2 = process->GlobalGetFullLogLikelihood();
            double dlnL = lnL2 - lnL1;
            if (std::isnan(dlnL))   {
                cerr << "nan lnl\n";
                cerr << lnL1 << '\t' << lnL2 << '\n';
                exit(1);
            }
            int dnsite = nsite2 - nsite;
            int dntaxa = ntaxa2 - ntaxa;
			los << cutoff << '\t' << dnsite << '\t' << dlnL << '\t' << dntaxa << '\n';
			// los << cutoff << '\t' << dnsite << '\t' << dlnL << '\t' << dlnL / dnsite << '\n';
			los.close();
		}	
		cerr << name << ": stopping after " << GetSize() << " points.\n";
		cerr << '\n';
	}

	NewickTree* GetTree() {return process->GetLengthTree();}

	void TraceHeader(ostream& os)	{
		stringstream ss;
		process->TraceHeader(ss);
		os << ss.str();
	}

	void Trace(ostream& os)	{
		stringstream ss;
        ss.precision(10);
		process->Trace(ss);
		os << ss.str();
	}
	
	void TreeTrace(ostream& os)	{
		stringstream ss;
		process->SetNamesFromLengths();
		process->RenormalizeBranchLengths();
		GetTree()->ToStream(ss);
		process->DenormalizeBranchLengths();
		os << ss.str();
	}
	
	void Monitor(ostream& os)	{
		stringstream ss;
		process->Monitor(ss);
		os << ss.str();
	}

	void ReadPB(int argc, char* argv[])	{
		process->ReadPB(argc,argv);
	}
};
