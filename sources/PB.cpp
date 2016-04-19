
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Model.h"

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

	if (! myid)	{
		cerr << '\n';
	}

	string datafile = "";
	string treefile = "None";
	string name = "";
	GeneticCodeType type = Universal;

	int every = 1;
	int until = -1;

	int mixturetype = -1;
	// int mixturetype = 3;
	int modeltype = -1;
	// int modeltype = 2;
	int dgam = 4;
	int ncat = 100;
	int iscodon = 0;
	int omegaprior = 0;

	int dc = 0;
	int fixtopo = 0;
	int fixcodonprofile = 1;
	int fixomega = 1;
	int NSPR = 10;
	int NNNI = 0;
	int fixbl = 0;
	int fixncomp = 0;
	int force = 0;
	int empmix = 0;
	string mixtype = "None";
	string rrtype = "None";

	int kappaprior = 0;
	int dirweightprior=0;
	// int betaprior = 0;

	int suffstat = 1;

	int saveall = 1;
	int incinit = 0;

	int burnin = 0;

	int randfix = -1;

	double mintotweight = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-v") || (s == "--version"))	{
				if (! myid)	{
					cerr << "\n";
					cerr << "pb_mpi version 1.5\n";
					cerr << "\n";
				}
				MPI_Finalize();
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
			else if (s == "-iscodon")	{
				iscodon = 1;
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
				if (s == "-T")	{
					fixtopo = 1;
				}
			}
			else if (s == "-Tbl")	{
				i++;
				treefile = argv[i];
				fixtopo = 1;
				fixbl = 1;
			}
			else if (s == "-spr")	{
				i++;
				NSPR = atoi(argv[i]);
			}
			else if (s == "-nni")	{
				i++;
				NNNI = atoi(argv[i]);
			}
			else if (s == "-fixcodonprofile")	{
				fixcodonprofile = 1;
			}
			else if (s == "-freecodonprofile")	{
				fixcodonprofile = 0;
			}
			else if (s == "-fixomega")	{
				fixomega = 1;
			}
			else if (s == "-freeomega")	{
				fixomega = 0;
			}
			else if (s == "-dc")	{
				dc = 1;
			}
			else if (s == "-s")	{
				saveall = 1;
			}
			else if (s == "-S")	{
				saveall = 0;
			}
			else if (s == "-priorinit")	{
				incinit = 0;
			}
			else if (s == "-incinit")	{
				i++;
				incinit = atoi(argv[i]);
				// incinit = 0;
			}
			else if ((s == "-poisson") || (s == "-f81"))	{
				modeltype = 1;
			}
			else if (s == "-gtr")	{
				modeltype = 2;
			}
			else if (s == "-mutselc")	{
				modeltype = 4;
			}
			else if ((s == "-mutsel") || (s == "-mutselaa") || (s == "-mutselaac"))	{
				modeltype = 5;
			}
			else if (s == "-dgam")	{
				i++;
				dgam = atoi(argv[i]);
			}
			/*
			else if (s == "-genpath")	{
				suffstat = 0;
			}
			else if (s == "-olddp")	{
				mixturetype = 2;
			}
			*/
			else if ((s == "-finite") || (s == "-ncat"))	{
				mixturetype = 1;
				i++;
				ncat = atoi(argv[i]);
				fixncomp = 1;
			}
			else if (s == "-fixncomp")	{
				fixncomp = 1;
			}
			else if (s == "-freencomp")	{
				fixncomp = 0;
			}
			else if (s == "-catfix")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				i++;
				mixtype = argv[i];
			}
			else if (s == "-rr")	{
				modeltype = 2;
				i++;
				rrtype = argv[i];
			}
			else if (s == "-lg")	{
				modeltype = 2;
				rrtype = "lg";
			}
			else if (s == "-wag")	{
				modeltype = 2;
				rrtype = "wag";
			}
			else if (s == "-jtt")	{
				modeltype = 2;
				rrtype = "jtt";
			}
			else if (s == "-mtzoa")	{
				modeltype = 2;
				rrtype = "mtzoa";
			}
			else if (s == "-mtrev")	{
				modeltype = 2;
				rrtype = "mtrev";
			}
			else if (s == "-mtart")	{
				modeltype = 2;
				rrtype = "mtart";
			}
			else if (s == "-mtzoa")	{
				modeltype = 2;
				rrtype = "mtzoa";
			}
			else if (s == "-cg6")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG6";
			}
			else if (s == "-cg10")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG10";
			}
			else if (s == "-cg20")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG20";
			}
			else if (s == "-cg30")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG30";
			}
			else if (s == "-cg40")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG40";
			}
			else if (s == "-cg50")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG50";
			}
			else if (s == "-cg60")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG60";
			}
			else if (s == "-cgr6")	{
				modeltype = 2;
				rrtype = "CG6";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG6";
			}
			else if (s == "-cgr10")	{
				modeltype = 2;
				rrtype = "CG10";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG10";
			}
			else if (s == "-cgr20")	{
				modeltype = 2;
				rrtype = "CG20";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG20";
			}
			else if (s == "-cgr30")	{
				modeltype = 2;
				rrtype = "CG30";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG30";
			}
			else if (s == "-cgr40")	{
				modeltype = 2;
				rrtype = "CG40";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG40";
			}
			else if (s == "-cgr50")	{
				modeltype = 2;
				rrtype = "CG50";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG50";
			}
			else if (s == "-cgr60")	{
				modeltype = 2;
				rrtype = "CG60";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG60";
			}
			
			else if ((s == "-dp") || (s == "-sbdp")	|| (s == "-cat")){
				mixturetype = 3;
			}
			/*
			else if (s == "-tdp")	{
				mixturetype = 4;
				i++;
				ncat = atoi(argv[i]);
			}
			*/
			else if (s == "-ss")	{
				mixturetype = 5;
			}
			else if (s == "-uni")	{
				type = Universal;
			}
			else if ((s == "-mtmam") || (s == "-MtMam") || (s == "mtvert") || (s == "MtVert"))	{
				type = MtMam;
			}
			else if (s == "-jeffkappa")	{
				kappaprior = 1;
			}
			else if (s == "-expkappa")	{
				kappaprior = 0;
			}
			else if (s == "-rigidbaseprior")	{
				dirweightprior = 1;
			}
			else if (s == "-mintotweight")	{
				i++;
				mintotweight = atof(argv[i]);
			}
			/*
			else if (s == "-jeffbeta")	{
				betaprior = 1;
			}
			else if (s == "-expbeta")	{
				betaprior = 0;
			}
			*/
			else if (s == "-jeffomega")	{
				omegaprior = 1;
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
		/*
		if ((datafile == "") && (argc != 2))	{
			throw(0);
		}
		*/
		if (nprocs <= 1)	{
			if (! myid)	{
				cerr << "error : pb_mpi requires at least 2 processes running in parallel (one master and at least one slave)\n";
			}
			MPI_Finalize();
			exit(1);
		}
	}
	catch(...)	{
		if (! myid)	{
			cerr << '\n';
			cerr << "mpirun -np <np> pb_mpi -d <datafile> [options] <chainname>\n";
			cerr << "\tcreates a new chain, sampling from the posterior distribution, conditional on specified data\n";
			cerr << "\n";
			cerr << "mpirun -np <np> pb_mpi <chainname>\n";
			cerr << "\tstarts an already existing chain\n";
			cerr << "\n";
			cerr << "\tmpirun -np <np>     : number of parallel processes (should be at least 2)\n";
			cerr << "\n";
			cerr << "\t-cat -dp            : infinite mixture (Dirichlet process) of equilibirium frequency profiles\n";
			cerr << "\t-ncat <ncat>        : finite mixture of equilibirium frequency profiles\n";
			cerr << "\t-catfix <pr>        : specifying a fixed pre-defined mixture of profiles\n";
			cerr << '\n';
			cerr << "\t-lg                 : Le and Gascuel 2008\n";
			cerr << "\t-wag                : Whelan and Goldman 2001\n";	
			cerr << "\t-jtt                : Jones, Taylor, Thornton 1992\n";	
			cerr << "\t-gtr                : general time reversible\n";
			cerr << "\t-poisson            : Poisson matrix, all relative exchangeabilities equal to 1 (Felsenstein 1981)\n";
			cerr << '\n';
			cerr << "\t-dgam <ncat>        : discrete gamma. ncat = number of categories (4 by default, 1 = uniform rates model)\n";
			cerr << '\n';
			cerr << "\t-dc                 : excludes constant columns\n";
			cerr << "\t-t <treefile>       : starts from specified tree\n"; 
			cerr << "\t-T <treefile>       : chain run under fixed, specified tree\n"; 
			cerr << '\n';
			cerr << "\t-x <every> <until>  : saving frequency, and chain length (until = -1 : forever)\n";
			cerr << "\t-f                  : forcing checks\n";
			cerr << "\t-s/-S               : -s : save all / -S : save only the trees\n";
			cerr << '\n';
			
			cerr << '\n';
			cerr << "see manual for details\n";
			cerr << '\n';

		}
		MPI_Finalize();
		exit(1);
	}

	if ((modeltype == -1) && (mixturetype == -1))	{
		modeltype = 2;
		mixturetype = 3;
	}
	else	{
		if (modeltype == -1)	{
			if (!myid)	{
			cerr << '\n';
			cerr << "error: incompletely specified model\n";
			cerr << "exchangeability parameters should be explicitly given\n";
			cerr << "-gtr -poisson (-f81) -lg -wag -jtt -mtrev -mtart -mtzoa or custom (-rr <filename>)\n";
			cerr << '\n';
			}
			MPI_Finalize();
			exit(1);
		}
		if (mixturetype == -1)	{
			if (!myid)	{
			cerr << '\n';
			cerr << "error: incompletely specified model\n";
			cerr << "mixture of equilibrium frequency profiles should be explicitly chosen among:\n";
			cerr << "-cat (or -dp) : infinite mixture (Dirichlet process)\n";
			cerr << "-ncat 1 : one matrix model\n";
			cerr << "-catfix <empmix>: empirical mixture (see manual for details)\n";
			cerr << '\n';
			}
			MPI_Finalize();
			exit(1);
		}
	}
	if (randfix != -1)	{
		rnd::init(1,randfix);
	}

	Model* model = 0;
	if (name == "")		{
		if (! myid)	{
			cerr << "error in command: no name was specified\n";
		}
		MPI_Finalize();
		exit(1);
	}
	if (datafile != "")	{
		if (myid == 0) {
			cerr << "model:\n";
			if (mixturetype == 1)	{
				if (empmix)	{
					cerr << "empirical mixture: " << mixtype << '\n';
				}
				else if (ncat == 1)	{
					cerr << "one-matrix model\n";
				}
				else	{
					cerr << "finite mixture of " << ncat << " components\n";
				}
			}
			else if (mixturetype == 3)	{
				cerr << "stick-breaking Dirichlet process mixture (cat)\n";
			}
			if (modeltype == 3)	{
				cerr << "codon mutation selection model\n";
			}
			else if (modeltype == 1)	{
				cerr << "exchangeabilities : f81 (Poisson)\n";
			}
			else if (modeltype == 2)	{
				if (rrtype == "None")	{
					cerr << "exchangeabilities estimated from data (gtr)\n";
				}
				else	{
					cerr << "exchangeabilities : " << rrtype << '\n';
				}
			}
			if (dgam == 1)	{
				cerr << "uniform rates across sites\n";
			}
			else	{
				if (modeltype != 3 && modeltype != 5)	{
					cerr << "discrete gamma distribution of rates across sites (" << dgam << " categories)\n";
				}
			}
			cerr << '\n';
		}
		if (! force)	{
			if (ifstream((name + ".param").c_str()))	{
				if (!myid)	{
					cerr << "a chain named " << name << " already exists; use -f to override\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
		}
		model = new Model(datafile,treefile,modeltype,dgam,mixturetype,ncat,type,suffstat,fixncomp,empmix,mixtype,rrtype,iscodon,fixtopo,NSPR,NNNI,fixcodonprofile,fixomega,fixbl,omegaprior,kappaprior,dirweightprior,mintotweight,dc,every,until,saveall,incinit,name,myid,nprocs);
		if (! myid)	{
			// cerr << "create files\n";
			cerr << '\n';
			cerr << "chain name : " << name << '\n';
			// MPI master only
			ofstream os((name + ".treelist").c_str());
			ofstream tos((name + ".trace").c_str());
			model->TraceHeader(tos);
			tos.close();
			ofstream pos((name + ".param").c_str());
			model->ToStream(pos,true);
			pos.close();
			if (saveall)	{
				ofstream cos((name + ".chain").c_str());
			}
			// cerr << "create files ok\n";
		}
	}
	else	{
		model = new Model(name,myid,nprocs);
		if (until != -1)	{
			model->until = until;
		}
	}

	if (myid == 0) {
		cerr << "run started\n";
		cerr << '\n';
		// model->Trace(cerr);
		model->Run(burnin);
		MESSAGE signal = KILL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else {
		// MPI slave
		model->WaitLoop();
	}
	MPI_Finalize();
}
