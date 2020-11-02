
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixDPProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"


double MatrixDPProfileProcess::IncrementalDPMove(int nrep)	{

	// UpdateOccupancyNumbers();
	int NAccepted = 0;
	int Nrep = (GetNsite() * nrep )/ 10;

	for (int rep=0; rep<Nrep; rep++)	{

		int site = (int) (GetNsite() * rnd::GetRandom().Uniform());

	// for (int site=0; site<GetNsite(); site++)	{
		int bk = alloc[site];
		int k = occupancy[alloc[site]] > 1 ? Ncomponent : Ncomponent-1;
		int h = k + Nadd;

		double* mLogSamplingArray = new double[h];

		// draw a new matrix for Nmode <= i < h
		for (int i=Ncomponent; i<h ; i++)	{
			CreateComponent(i);
			// already called in CreateComponent
			// SampleStat(i);
		}

		RemoveSite(site,bk);

		// Gibbs

		double max = 0;
		for (int mode = 0; mode < h; mode++)	{
			mLogSamplingArray[mode] =  LogStatProb(site,mode);
			if ((!mode) || (max < mLogSamplingArray[mode]))	{
				max = mLogSamplingArray[mode];
			}
		}

		double* cumul = new double[h];
		double total = 0;

		for (int mode = 0; mode < h; mode++)	{

			double p = 0;			
			if (mode < Ncomponent)	{			
				if (occupancy[mode])	{
					p = (double) occupancy[mode];
				}
				else	{
					p = kappa / Nadd;
				}
				p *= exp(mLogSamplingArray[mode] - max);
			}
			else	{
				p = (kappa / Nadd) * exp(mLogSamplingArray[mode] - max);		
			}

			total += p;
			cumul[mode] = total;
		}

		double q = total * rnd::GetRandom().Uniform();
		int mode = 0;
		while ( (mode<h) && (q > cumul[mode])) mode++;
		if (mode == h)	{
			cerr << "error in switch mode: gibbs overflow\n";
			exit(1);
		}
		delete[] cumul;

		int Accepted = (mode != bk);
		if (Accepted)	{
			NAccepted ++;
		}
		AddSite(site,mode);

		if (mode >= Ncomponent)	{			// if it's a new one
			if (mode > Ncomponent)	{
				SwapComponents(mode, Ncomponent);
				mode = Ncomponent;
			}
			Ncomponent++;
		}
		if (! occupancy[bk])	{
			if (bk!= Ncomponent-1)	{
				SwapComponents(bk, Ncomponent-1);
			}
			Ncomponent--;
		}

		for (int k=Ncomponent; k<h; k++)	{
			DeleteComponent(k);
		}

		delete[] mLogSamplingArray;
	}
	
	// UpdateModeProfileSuffStat();
	return ((double) NAccepted) / GetNsite() / nrep * 10;
}

void MatrixDPProfileProcess::SlaveMixMove()	{

	int h = Ncomponent + Nadd*Ninc;
	double* logsamp = new double[h];
	int size = (GetSiteMax() - GetSiteMin()) * (h + 1 + Nadd*Ninc*GetDim());
	double* alloc = new double[size];

	int index = 0;
	for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{
	
		for (int i=Ncomponent; i<h; i++)	{
			CreateComponent(i);
			double totstat = 0;
			for (int k=0; k<GetDim(); k++)	{
				alloc[index] = profile[i][k];
				totstat += alloc[index];
				index++;
			}
			if (fabs(totstat - 1) > 1e-6)	{
				cerr << "slave : normalisation error : " << totstat << '\n';
			}
			UpdateComponent(i);
		}

		double max = 0;
		for (int k=0; k<h; k++)	{
			logsamp[k] = LogStatProb(site,k);
			if ((!k) || (max < logsamp[k]))	{
				max = logsamp[k];
			}
		}
		for (int k=0; k<h; k++)	{
			alloc[index] = exp(logsamp[k] - max);
			index++;
		}
		alloc[index] = max;
		index++;
		
		for (int i=Ncomponent; i<h; i++)	{
			DeleteComponent(i);
		}
	}

	if (index != size)	{
		cerr << "error in slave: non matching dim\n";
		exit(1);
	}

	MPI_Send(alloc,size,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	delete[] logsamp;
	delete[] alloc;
}

double MatrixDPProfileProcess::GlobalMixMove(int Nrep, int Nprofile)	{

	int Nextra = 1000;
	int Nmax = Ncomponent + Nextra;
	double** prob = new double*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		prob[i] = new double[Nmax];
		for (int j=0; j<Nmax; j++)	{
			prob[i][j] = -1;
		}
	}
	
	double** addprob = new double*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		addprob[i] = new double[Nadd*Ninc];
	}
	
	double* offset = new double[GetNsite()];

	double*** stat = new double**[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		stat[i] = new double*[Nadd * Ninc];
		for (int j=0; j<Nadd*Ninc; j++)	{
			stat[i][j] = new double[GetDim()];
		}
	}

	UpdateOccupancyNumbers();

	// split Nsite among GetNprocs()-1 slaves
	int width = GetNsite()/(GetNprocs()-1);
	int maxw = 0;
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
		int w = smax[i] - smin[i];
		if (maxw < w)	{
			maxw = w;
		}
	}
	for (int rep=0; rep<Nrep; rep++)	{


		assert(GetMyid() == 0);

		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		UpdateModeProfileSuffStat();
		// send PROFILE_MOVE Message with n and nrep and tuning
		
		MESSAGE signal = MIX_MOVE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		// mpi send message
		// mpi send Nmode
		// mpi send stats

		// mpi receive probs
		// move additional probs into addprobs
		// set all non valid probs to -1
		// mpi receive additional stats
		int h = Ncomponent + Nadd*Ninc;
		if (h >= Nmax)	{
			cerr << "overflow error\n";
			exit(1);
		}	
		int maxsize = maxw * (h + 1 + Nadd*Ninc*GetDim());
		double* tmpalloc = new double[maxsize];
		MPI_Status status;
		for(int i=1; i<GetNprocs(); ++i) {
			
			int size = (smax[i-1] - smin[i-1]) * (h + 1 + Nadd*Ninc*GetDim());
			MPI_Recv(tmpalloc,size,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&status);

			int index = 0;
			for (int site=smin[i-1]; site<smax[i-1]; site++)	{
				for (int j=0; j<Nadd*Ninc; j++)	{
					for (int k=0; k<GetDim(); k++)	{
						stat[site][j][k] = tmpalloc[index];
						index++;
					}
				}
				for (int k=0; k<Ncomponent; k++)	{
					prob[site][k] = tmpalloc[index];
					index++;
				}
				for (int j=0; j<Nadd*Ninc; j++)	{
					addprob[site][j] = tmpalloc[index];
					index++;
				}
				offset[site] = tmpalloc[index];
				index++;
			}
			if (index != size)	{
				cerr << "error : non matching alloc size\n";
				exit(1);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		delete[] tmpalloc;


		// checks
		for (int i=0; i<GetNsite(); i++)	{
			for (int j=0; j<Ncomponent; j++)	{
				if (prob[i][j] == -1)	{
					cerr << "prob is -1: " << i << '\t' << j << '\n';
					cerr << Ncomponent << '\n';
					exit(1);
				}
			}
		}

		for (int i=0; i<GetNsite(); i++)	{
			for (int j=0; j<Nadd*Ninc; j++)	{
				if (addprob[i][j] == -1)	{
					cerr << "add prob is -1\n";
					cerr << i << '\t' << j << '\n';
					exit(1);
				}
			}
			if (offset[i] == -1)	{
				cerr << "offset is -1\n";
				exit(1);
			}
		}

		// reset all unoccupied components
		for (int i=0; i<GetNsite(); i++)	{
			for (int j=Ncomponent; j<Nmax; j++)	{
				prob[i][j] = -1;
			}
		}

		int NAccepted = 0;

		for (int inc=0; inc<Ninc; inc++)	{

			for (int siterep=0; siterep<GetNsite(); siterep++)	{
				// int site = (int) (GetNsite() * rnd::GetRandom().Uniform());
				int site = siterep;

				int bk = alloc[site];
				int k = occupancy[alloc[site]] > 1 ? Ncomponent : Ncomponent-1;
				int h = k + Nadd;
				if (h == Nmax)	{
					cerr << "overflow in global mix move\n";
					exit(1);
				}

				// draw a new component for Nmode <= i < h
				for (int i=Ncomponent; i<h ; i++)	{
					CreateComponent(i);
					double totstat = 0;
					for (int k=0; k<GetDim(); k++)	{
						profile[i][k] = stat[site][Nadd*inc + i-Ncomponent][k];
						totstat += profile[i][k];
					}
					if (std::isnan(totstat))	{
						cerr << "totstat is nan\n";
						exit(1);
					}
					if (fabs(totstat - 1) > 1e-6)	{
						cerr << "normalisation error\n";
						cerr << totstat << '\n';
						exit(1);
					}
					prob[site][i] = addprob[site][Nadd*inc + i-Ncomponent];
					UpdateComponent(i);
				}

				RemoveSite(site,bk);

				// Gibbs
				double max = offset[site];
				double* cumul = new double[h];
				double total = 0;
				for (int mode = 0; mode < h; mode++)	{
					if (prob[site][mode] == -1)	{
						prob[site][mode] = exp(LogStatProb(site,mode) - max);
						if (std::isnan(prob[site][mode]))	{
							cerr << site << '\t' << mode << '\t' << Ncomponent << '\t' << prob[site][mode] << '\t' << LogStatProb(site,mode) << '\n';
							for (int k=0; k<GetDim(); k++)	{
								cerr << profile[mode][k] << '\n';
							}
							exit(1);
						}
					}
					// check
					/*
					else	{
						if (prob[site][mode] != exp(LogStatProb(site,mode) - max))	{
							cerr << "error: non matching prob\n";
							cerr << prob[site][mode] << '\t' << exp(LogStatProb(site,mode) - max) << '\n';
							cerr << site << '\t' << mode << '\t' << max << '\n';
							exit(1);
						}
					}
					*/

					double p = 0;			
					if (occupancy[mode])	{
						p = occupancy[mode] * prob[site][mode];
					}
					else	{
						p = kappa / Nadd * prob[site][mode];
					}
					total += p;
					cumul[mode] = total;
				}

				double q = total * rnd::GetRandom().Uniform();
				int mode = 0;
				while ( (mode<h) && (q > cumul[mode])) mode++;
				if (mode == h)	{
					cerr << "error in switch mode: gibbs overflow\n";
					exit(1);
				}
				delete[] cumul;

				int Accepted = (mode != bk);
				if (Accepted)	{
					NAccepted ++;
				}
				AddSite(site,mode);

				if (mode >= Ncomponent)	{			// if it's a new one
					if (mode > Ncomponent)	{
						SwapComponents(mode, Ncomponent);
						for (int l=0; l<GetNsite(); l++)	{
							double tmp = prob[l][mode];
							prob[l][mode] = prob[l][Ncomponent];
							prob[l][Ncomponent] = tmp;
						}
						mode = Ncomponent;
					}
					Ncomponent++;
				}
				if (! occupancy[bk])	{
					if (bk!= Ncomponent-1)	{
						SwapComponents(bk, Ncomponent-1);
						for (int l=0; l<GetNsite(); l++)	{
							// double tmp = prob[l][bk];
							prob[l][bk] = prob[l][Ncomponent-1];
							prob[l][Ncomponent-1] = -1;
						}
					}
					Ncomponent--;
				}

				for (int k=Ncomponent; k<h; k++)	{
					DeleteComponent(k);
					prob[site][k] = -1;
				}
			}
		}

		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		UpdateModeProfileSuffStat();

		GlobalMoveProfile(1,1,Nprofile);
		GlobalMoveProfile(1,3,Nprofile);
		GlobalMoveProfile(0.1,3,Nprofile);

	}

	for (int i=0; i<GetNsite(); i++)	{
		delete[] prob[i];
	}
	delete[] prob;
	for (int i=0; i<GetNsite(); i++)	{
		delete[] addprob[i];
	}
	delete[] addprob;
	
	delete[] offset;

	for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<Nadd*Ninc; j++)	{
			delete[] stat[i][j];
		}
		delete[] stat[i];
	}
	delete[] stat;

	return 1;
}

