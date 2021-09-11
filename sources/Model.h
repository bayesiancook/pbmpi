
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
    int steppingdnsite;
    int steppingburnin;
    int steppingsize;
    double steppingmaxvar;
    int steppingmaxsize;
    string empstepping;
    double empramp;
    int steppingcycle;
    int randstepping;

	Model(string datafile, string treefile, int modeltype, int nratecat, int mixturetype, int ncat, int nmodemax, GeneticCodeType codetype, int suffstat, int fixncomp, int empmix, string mixtype, string rrtype, int iscodon, int fixtopo, int NSPR, int NNNI, int fixcodonprofile, int fixomega, int fixbl, int omegaprior, int kappaprior, int dirweightprior, double mintotweight, int dc, int inevery, int inuntil, int insaveall, int inincinit, int topoburnin, int insteppingdnsite, int insteppingburnin, int insteppingsize, double insteppingmaxvar, int insteppingmaxsize, int inrandstepping, string inempstepping, double inempramp, string inname, int myid, int nprocs)	{

		every = inevery;
		until = inuntil;
		name = inname;
		saveall = insaveall;
		incinit = inincinit;
        steppingdnsite = insteppingdnsite;
        steppingburnin = insteppingburnin;
        steppingsize = insteppingsize;
        steppingmaxvar = insteppingmaxvar;
        steppingmaxsize = insteppingmaxsize;
        randstepping = inrandstepping;
        empstepping = inempstepping;
        empramp = inempramp;
        steppingcycle = 0;

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

        steppingdnsite = 0;

		name = inname;

		ifstream is((name + ".param").c_str());
		if (! is)	{
			cerr << "error: cannot open " << name << ".param\n";
			exit(1);
		}

		is >> type;
        if (type == "STEPPING") {
            is >> steppingdnsite >> steppingburnin >> steppingsize >> steppingmaxvar >> steppingmaxsize;
            is >> empstepping;
            is >> empramp;
            is >> steppingcycle;
            is >> randstepping;
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
            if (steppingdnsite)  {
                ss << "STEPPING\n";
                ss << steppingdnsite << '\t' << steppingburnin << '\t' << steppingsize << '\t' << steppingmaxvar << '\t' << steppingmaxsize << '\n';
                ss << empstepping << '\n';
                ss << empramp << '\n';
                ss << steppingcycle << '\n';
                ss << randstepping << '\n';
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
        if (! steppingdnsite)    {
            MCMCRun(burnin);
        }
        else    {
            SteppingRun(steppingdnsite, steppingburnin, steppingsize,
                        steppingmaxvar, steppingmaxsize, randstepping, empstepping, empramp);
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

            if (! process->fixtopo) {
                ofstream os((name + ".treelist").c_str(), ios_base::app);
                TreeTrace(os);
                os.close();
            }

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

	void SteppingRun(int step, int burnin, int minnpoint, double maxvar, double maxnpoint, int rand, string empname, double empramp)    {


        // deprecated option for importance sampling estimation of mixture log likelihood
        // deactivated when nrep = 0
        int nrep = 0;

        int empiricalprior = 0;
        if (empname != "None")  {
            empiricalprior = 1;
            ifstream is(empname.c_str());
            process->GlobalSetEmpiricalPrior(is);
        }
		process->GlobalPrepareStepping(name, GetSize(), rand);
	
		ofstream ros((name + ".run").c_str()); stringstream buf;
		buf << 1 << '\n';
		ros << buf.str();
		ros.close();
	
        if (! GetSize())    {
            if (empiricalprior)  {
                process->GlobalSetEmpiricalFrac(0);
            }
            process->PriorSample();
            process->GlobalUpdateParameters();
        }

        process->GlobalResetAllConditionalLikelihoods();

        int ncycle = process->GetNsite() / step;
        if (process->GetNsite() % step) {
            ncycle++;
        }

		while (RunningStatus() && (steppingcycle < ncycle)) {

            double frac1 = ((double) steppingcycle) / ncycle;
            frac1 *= empramp;
            if (frac1 > 1.0)    {
                frac1 = 1.0;
            }
            int nsite1 = steppingcycle * step;
            if (nsite1 > process->GetNsite()) {
                nsite1 = process->GetNsite();
            }

            double frac2 = ((double) steppingcycle + 1.0) / ncycle;
            frac2 *= empramp;
            if (frac2 > 1.0)    {
                frac2 = 1.0;
            }
            int nsite2 = (steppingcycle + 1) * step;
            if (nsite2 > process->GetNsite()) {
                nsite2 = process->GetNsite();
            }

            process->GlobalSetSteppingFraction(0, nsite1);

            // not necessary: increasing series of sites
            // process->GlobalResetAllConditionalLikelihoods();

            process->GlobalUpdateConditionalLikelihoods();

            if (empiricalprior)  {
                process->GlobalSetEmpiricalFrac(frac1);
            }

            for (int i=0; i<burnin; i++)    {
                Move(1,every);
                process->IncSize();
            }

            double premaxlogp = 0;
            double pretotp1 = 0;
            double pretotp2 = 0;
            double pretotlogp1 = 0;
            double pretotlogp2 = 0;
            double pretotlogprior = 0;

            double prelogZ = 0;
            double preeffsize = 0;
            double premeanlogp = 0;
            double prevarlogp = 0;
            double premeanlogprior = 0;

            double targetnpoint = minnpoint;

            if (maxvar) {
                int npoint = 0;
                while (npoint < minnpoint)  {
                    Move(1,every);
                    process->IncSize();
                    npoint++;

                    process->GlobalSetSteppingFraction(nsite1, nsite2);
                    double delta = process->GlobalGetSteppingLogLikelihood(nrep, 1);

                    double dlogp = 0;
                    if (empiricalprior)  {
                        double lnP1 = process->GetLogPrior();
                        process->GlobalSetEmpiricalFrac(frac2);
                        double lnP2 = process->GetLogPrior();
                        delta += lnP2 - lnP1;
                        dlogp = lnP2 - lnP1;
                    }
                    if (std::isnan(delta))   {
                        cerr << "nan delta\n";
                        exit(1);
                    }

                    pretotlogprior += dlogp;
                    pretotlogp1 += delta;
                    pretotlogp2 += delta*delta;
                    if ((!premaxlogp) || (premaxlogp < delta))    {
                        pretotp1 *= exp(premaxlogp-delta);
                        pretotp1 += 1.0;
                        pretotp2 *= exp(2*(premaxlogp-delta));
                        pretotp2 += 1.0;
                        premaxlogp = delta;
                    }
                    else    {
                        pretotp1 += exp(delta - premaxlogp);
                        pretotp2 += exp(2*(delta - premaxlogp));
                    }

                    if (npoint == minnpoint)    {
                        prelogZ = log(pretotp1 / npoint) + premaxlogp;
                        preeffsize = pretotp1 * pretotp1 / pretotp2;
                        premeanlogp = pretotlogp1/npoint;
                        prevarlogp = pretotlogp2/npoint - premeanlogp*premeanlogp;
                        premeanlogprior = pretotlogprior / npoint;
                    }

                    process->GlobalSetSteppingFraction(0, nsite1);
                    process->GlobalSetEmpiricalFrac(frac1);
                }

                if (prevarlogp > maxvar)   {
                    targetnpoint *= exp(prevarlogp)/exp(maxvar);
                    if (targetnpoint > maxnpoint) {
                        targetnpoint = maxnpoint;
                    }
                }
            }

            int finalnpoint = int(targetnpoint);

            double maxlogp = 0;
            double totp1 = 0;
            double totp2 = 0;
            double totlogp1 = 0;
            double totlogp2 = 0;
            double totlogprior = 0;

            int npoint = 0;
            while (npoint < finalnpoint)    {
                Move(1,every);
                process->IncSize();
                npoint++;

                process->GlobalSetSteppingFraction(nsite1, nsite2);
                int restore = (npoint == finalnpoint) ? 0 : 1;
                double delta = process->GlobalGetSteppingLogLikelihood(nrep, restore);
                double dlogp = 0;
                if (empiricalprior)  {
                    double lnP1 = process->GetLogPrior();
                    process->GlobalSetEmpiricalFrac(frac2);
                    double lnP2 = process->GetLogPrior();
                    delta += lnP2 - lnP1;
                    dlogp = lnP2 - lnP1;
                }
                if (std::isnan(delta))   {
                    cerr << "nan delta\n";
                    exit(1);
                }

                totlogprior += dlogp;

                totlogp1 += delta;
                totlogp2 += delta*delta;
                if ((!maxlogp) || (maxlogp < delta))    {
                    totp1 *= exp(maxlogp-delta);
                    totp1 += 1.0;
                    totp2 *= exp(2*(maxlogp-delta));
                    totp2 += 1.0;
                    maxlogp = delta;
                }
                else    {
                    totp1 += exp(delta - maxlogp);
                    totp2 += exp(2*(delta - maxlogp));
                }

                if (npoint < finalnpoint)   {
                    process->GlobalSetSteppingFraction(0, nsite1);
                    process->GlobalSetEmpiricalFrac(frac1);
                }
                else    {

                    double logZ = log(totp1 / npoint) + maxlogp;
                    double effsize = totp1 * totp1 / totp2;
                    double meanlogp = totlogp1/npoint;
                    double varlogp = totlogp2/npoint - meanlogp*meanlogp;
                    double meanlogprior = totlogprior / npoint;

                    ofstream los((name + ".stepping").c_str(), ios_base::app);
                    if (maxvar) {
                        los << frac1 << '\t' << nsite1 << '\t' << logZ << '\t' << meanlogp << '\t' << meanlogprior << '\t' << varlogp << '\t' << npoint << '\t' << effsize << '\t' << prelogZ << '\t' << premeanlogp << '\t' << premeanlogprior << '\t' << prevarlogp << '\t' << minnpoint << '\t' << preeffsize <<  '\n';
                    }
                    else    {
                        los << frac1 << '\t' << nsite1 << '\t' << logZ << '\t' << meanlogp << '\t' << meanlogprior << '\t' << varlogp << '\t' << npoint << '\t' << effsize << '\n';
                    }
                    los.close();

                    steppingcycle++;
                    ofstream pos((name + ".param").c_str());
                    pos.precision(numeric_limits<double>::digits10);
                    ToStream(pos,true);
                    pos.close();
                }
            }
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
