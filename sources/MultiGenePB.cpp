
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MultiGeneGTRSBDPMixture.h"
#include "MultiGenePoissonSBDPMixture.h"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneModel	{

	public:

	// MultiGeneGTRSBDPMixture* process;
	MultiGeneMixture* process;

	string type;
	string name;
	int every;
	int until;
	int size;
	int saveall;

	MultiGeneModel(string datafile, string treefile, int modeltype, int nratecat, int fixtopo, int kappaprior, double gibbsfactor, int dc, int inevery, int inuntil, int insaveall, string inname, int myid, int nprocs)	{

		every = inevery;
		until = inuntil;
		name = inname;
		saveall = insaveall;
		size = 0;

		if (modeltype == 1)	{
			if (myid == 0) {
				cerr << "cat poisson model\n";
			}
			type = "CATPOISSON";
			process = new MultiGenePoissonSBDPMixture(datafile,treefile,name,nratecat,fixtopo,kappaprior,gibbsfactor,dc,myid,nprocs); 
		}

		else if (modeltype == 2)	{
			if (myid == 0) {
				cerr << "catgtr model\n";
			}
			type = "CATGTRSBDP";
			process = new MultiGeneGTRSBDPMixture(datafile,treefile,name,nratecat,fixtopo,kappaprior,gibbsfactor,dc,myid,nprocs); 
		}
	}

	MultiGeneModel(string inname, int myid, int nprocs)	{

		name = inname;

		ifstream is((name + ".param").c_str());
		if (! is)	{
			cerr << "error: cannot open " << name << ".param\n";
			exit(1);
		}

		is >> type;
		is >> every >> until >> size;
		is >> saveall;
		
		if (type == "CATPOISSON")	{
			process = new MultiGenePoissonSBDPMixture(is,myid,nprocs); 
		}
		else if (type == "CATGTRSBDP")	{
			process = new MultiGeneGTRSBDPMixture(is,myid,nprocs); 
		}
		else	{
			cerr << "error, does not recognize model type : " << type << '\n';
			exit(1);
		}

		// cerr << "RESTORE SETSIZE\n";
		// process->SetSize(size);
		// cerr << "reset size to " << process->GetSize() << '\n';
	}

	void ToStream(ostream& os, bool header)	{
		if (header)	{
			os << type << '\n';
			os << every << '\t' << until << '\t' << size << '\n';
			os << saveall << '\n';
			process->ToStreamHeader(os);
		}
		process->ToStream(os);
	}

	~MultiGeneModel()	{
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
		return size;
	}

	void Run()	{

		ofstream ros((name + ".run").c_str());
		ros << 1 << '\n';
		ros.close();
	
		while (RunningStatus() && ((until == -1) || (size < until)))	{

			Move(1,every);
			
			size++;

			ofstream tos((name + ".trace").c_str(), ios_base::app);
			Trace(tos);
			tos.close();

			/*
			ofstream mos((name + ".monitor").c_str());
			process->Monitor(mos);
			mos.close();
			*/

			ofstream pos((name + ".param").c_str());
			pos.precision(10);
			ToStream(pos,true);
			pos.close();

			process->SaveTrees();

			if (saveall)	{
				ofstream cos((name + ".chain").c_str(),ios_base::app);
				ToStream(cos,false);
				cos.close();
			}

		}	
		cerr << name << ": stopping after " << GetSize() << " points.\n";
	}

	void TraceHeader(ostream& os)	{
		process->TraceHeader(os);
	}

	void Trace(ostream& os)	{
		process->Trace(os);
	}
};

int main(int argc, char* argv[])	{

	int myid,nprocs;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int blockcounts[2] = {1,3};
	MPI_Datatype types[2] = {MPI_DOUBLE,MPI_INT};
	MPI_Aint dtex,displacements[2];
	
	displacements[0] = (MPI_Aint) 0;
	MPI_Type_extent(MPI_DOUBLE,&dtex);
	displacements[1] = dtex;
	MPI_Type_struct(2,blockcounts,displacements,types,&Propagate_arg);
	MPI_Type_commit(&Propagate_arg); 

	string datafile = "";
	string treefile = "None";
	string name = "";

	int modeltype = 2;

	int every = 1;
	int until = -1;

	int dgam = 4;
	int dc = 0;
	int fixtopo = 0;
	int force = 0;

	int kappaprior = 0;
	double gibbsfactor = 10;

	int saveall = 0;
	int burnin = 0;
	int randfix = -1;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-v") || (s == "--version"))	{
				cerr << "\n";
				cerr << "pb_mpi version 1.5\n";
				cerr << "\n";
				exit(1);
			}
			else if ((s == "-h") || (s == "--help"))	{
				throw(0);
			}
			else if (s == "-f")	{
				force = 1;
			}
			else if (s == "-rnd")	{
				i++;
				randfix = atoi(argv[i]);
			}
			else if (s == "-d")	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-gibbs")	{
				i++;
				gibbsfactor = atof(argv[i]);
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
				if (s == "-T")	{
					fixtopo = 1;
				}
			}
			else if (s == "-dc")	{
				dc = 1;
			}
			else if (s == "-s")	{
				saveall = 1;
			}
			else if (s == "-poisson")	{
				modeltype = 1;
			}
			else if (s == "-gtr")	{
				modeltype = 2;
			}
			else if (s == "-dgam")	{
				i++;
				dgam = atoi(argv[i]);
			}
			else if (s == "-jeffkappa")	{
				kappaprior = 1;
			}
			else if (s == "-expkappa")	{
				kappaprior = 0;
			}
			else if (s == "-b")	{
				i++;
				burnin = atoi(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				until = atoi(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (nprocs <= 1)	{
			cerr << "error : pb_mpi requires at least 2 processes running in parallel (one master and at least one slave)\n";
			exit(1);
		}
	}
	catch(...)	{
		if (! myid)	{
			cerr << '\n';
			cerr << "error in command\n";
			exit(1);
		}
	}

	if (randfix != -1)	{
		rnd::init(1,randfix);
	}

	if (name == "")		{
		if (! myid)	{
			cerr << "error in command: no name was specified\n";
		}
		exit(1);
	}

	MultiGeneModel* model = 0;

	if (datafile != "")	{
		model = new MultiGeneModel(datafile,treefile,modeltype,dgam,fixtopo,kappaprior,gibbsfactor,dc,every,until,saveall,name,myid,nprocs);
		if (! myid)	{
			cerr << "create files\n";
			cerr << name << '\n';
			// MPI master only
			ofstream tos((name + ".trace").c_str());
			model->TraceHeader(tos);
			tos.close();
			ofstream pos((name + ".param").c_str());
			model->ToStream(pos,true);
			pos.close();
			if (saveall)	{
				ofstream cos((name + ".chain").c_str());
			}
			cerr << "create files ok\n";
		}
	}
	else	{
		model = new MultiGeneModel(name,myid,nprocs);
		if (until != -1)	{
			model->until = until;
		}
	}

	if (myid == 0) {
		cerr << "start\n";
		model->Trace(cerr);
		model->Run();
		MESSAGE signal = KILL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else {
		// MPI slave
		model->WaitLoop();
	}
	MPI_Finalize();
}
