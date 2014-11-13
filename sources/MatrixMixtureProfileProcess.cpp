
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixMixtureProfileProcess.h"
#include "Random.h"
#include <cassert>
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MatrixMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void MatrixMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);

	UpdateMatrix(cat1);
	UpdateMatrix(cat2);

	// useful?
	// in expo gtr: null pointers anyway
	// with GeneralSuffStat: might be incorrect
	// CHECK

	/*
	SubMatrix* tmp = matrixarray[cat1];
	matrixarray[cat1] = matrixarray[cat2];
	matrixarray[cat2] = tmp;
	*/

}

void MatrixMixtureProfileProcess::Create(int innsite, int indim)	{
	if (! matrixarray)	{
		MixtureProfileProcess::Create(innsite,indim);
		matrixarray = new SubMatrix*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			matrixarray[i] = 0;
		}
		// SampleProfile();
	}
}

void MatrixMixtureProfileProcess::Delete() {
	if (matrixarray)	{
		for (int i=0; i<GetNmodeMax(); i++)	{
			delete matrixarray[i];
		}
		delete[] matrixarray;
		matrixarray = 0;
		MixtureProfileProcess::Delete();
	}
}

double MatrixMixtureProfileProcess::GlobalMoveProfile(double tuning, int n, int nrep)	{

	UpdateOccupancyNumbers();

	assert(GetMyid() == 0);

	// send PROFILE_MOVE Message with n and nrep and tuning
	MESSAGE signal = PROFILE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int* itmp = new int[3+GetNsite()];
	itmp[0] = n;
	itmp[1] = nrep;
	itmp[2] = Ncomponent;
	for (int i=0; i<GetNsite(); i++)	{
		itmp[3+i] = alloc[i];
	}
	MPI_Bcast(itmp,3+GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
	delete[] itmp;

	int Nocc = GetNOccupiedComponent();

	double* dtmp = new double[1+Nocc*GetDim()];
	dtmp[0] = tuning;
	int k = 1;
	for (int i=0; i<Ncomponent; i++)	{
		if (occupancy[i])	{
			for (int j=0; j<GetDim(); j++)	{
				dtmp[k] = profile[i][j];
				k++;
			}
		}
	}
	MPI_Bcast(dtmp,1+Nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	delete[] dtmp;

	// split Ncomponent items among GetNprocs() - 1 slaves
	int width = Nocc/(GetNprocs()-1);
	int maxwidth = 0;
	int cmin[GetNprocs()-1];
	int cmax[GetNprocs()-1];
	int dmin[GetNprocs()-1];
	int dmax[GetNprocs()-1];

	for(int i=0; i<GetNprocs()-1; ++i) {

		int ddmin = width * i;
		int ddmax = (i == GetNprocs() - 2) ? Nocc : width * (i+1);
		if (maxwidth < (ddmax - ddmin))	{
			maxwidth = ddmax - ddmin;
		}

		int k = -1;
		int ccmin = -1;
		while ((ccmin<Ncomponent) && (k<ddmin))	{
			ccmin++;
			if (ccmin == Ncomponent)	{
				cerr << "error in matmixslavemoveprofile: overflow\n";
				exit(1);
			}
			if (occupancy[ccmin])	{
				k++;
			}
		}
		int ccmax = ccmin;
		if (ddmax == Nocc)	{
			ccmax = Ncomponent;
		}
		else	{
			while ((ccmax<Ncomponent) && (k<ddmax))	{
				ccmax++;
				if (occupancy[ccmax])	{
					k++;
				}
			}
		}

		cmin[i] = ccmin;
		cmax[i] = ccmax;
		dmin[i] = ddmin;
		dmax[i] = ddmax;
		int nocc = 0;
		for (int j=ccmin; j<ccmax; j++)	{
			if (occupancy[j])	{
				nocc++;
			}
		}
		if (nocc != (dmax[i] - dmin[i]))	{
			cerr << "error: non matching numbers: " << i << '\t' << nocc << '\t' << ccmin << '\t' << ccmax << '\t' << ddmin << '\t' << ddmax << '\n';
			for (int j=ccmin; j<ccmax; j++)	{
				cerr << occupancy[j] << '\t';
			}
			cerr << '\n';
			exit(1);
		}
	}

	/*
	// split Ncomponent items among GetNprocs() - 1 slaves
	int width = Nocc/(GetNprocs()-1);
	// int width = GetNcomponent()/(GetNprocs()-1);
	int maxwidth = 0;
	int cmin[GetNprocs()-1];
	int cmax[GetNprocs()-1];
	int dmin[GetNprocs()-1];
	int dmax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		// determine the range of components to move
		int ddmin = width * i;
		int ddmax = (i == GetNprocs() - 2) ? Nocc : width * (i+1);
		if (maxwidth < (ddmax - ddmin))	{
			maxwidth = ddmax - ddmin;
		}

		k = 0;
		int ccmin = 0;
		while ((ccmin<Ncomponent) && (k<ddmin))	{
			ccmin++;
			if (ccmin == Ncomponent)	{
				cerr << "error in matmixslavemoveprofile: overflow\n";
				exit(1);
			}
			if (occupancy[ccmin])	{
				k++;
			}
		}
		int ccmax = ccmin;
		if (ddmax == Nocc)	{
			ccmax = Ncomponent;
		}
		else	{
			while ((ccmax<Ncomponent) && (k<ddmax))	{
				ccmax++;
				if (ccmax == Ncomponent)	{
					cerr << "error in matmixslavemoveprofile: overflow\n";
					exit(1);
				}
				if (occupancy[ccmax])	{
					k++;
				}
			}
		}

		// cerr << ":::" << '\t' << i+1 << '\t' << ccmin << '\t' << ccmax << '\t' << dmin << '\t' << dmax << '\n';
		cmin[i] = ccmin;
		cmax[i] = ccmax;
		dmin[i] = ddmin;
		dmax[i] = ddmax;
	}
	*/

	// collect final values of profiles (+ total acceptance rate) from slaves
	MPI_Status stat;
	int bigdim = maxwidth * GetDim();
	double* tmp = new double[bigdim+1]; // (+1 for the acceptance rate)
	double total = 0;
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(tmp,(dmax[i-1]-dmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		int l = 0;
		for(int j=cmin[i-1]; j<cmax[i-1]; ++j) {
			if (occupancy[j])	{
				for (int k=0; k<GetDim(); k++)	{
					profile[j][k] = tmp[l];
					l++;
				}
			}
		}
		total += tmp[l]; // (sum all acceptance rates)
	}
	delete[] tmp;
	// return average acceptance rate
	return total / GetNOccupiedComponent();
	// return total / GetNOccupiedComponent();
}

void MatrixMixtureProfileProcess::SlaveMoveProfile()	{

	assert(GetMyid() > 0);

	// parse arguments sent by master

	int* itmp = new int[3+GetNsite()];
	MPI_Bcast(itmp,3+GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
	int n = itmp[0];
	int nrep = itmp[1];
	Ncomponent = itmp[2];
	for (int i=0; i<GetNsite(); i++)	{
		alloc[i] = itmp[3+i];
	}
	delete[] itmp;

	UpdateOccupancyNumbers();
	int Nocc = GetNOccupiedComponent();

	double* dtmp = new double[1 + Nocc*GetDim()];
	MPI_Bcast(dtmp,1+Nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	double tuning = dtmp[0];
	int k = 1;
	for (int i=0; i<Ncomponent; i++)	{
		if (occupancy[i])	{
			for (int j=0; j<GetDim(); j++)	{
				profile[i][j] = dtmp[k];
				k++;
			}
		}
	}
	delete[] dtmp;

	// determine the range of components to move
	int width = Nocc/(GetNprocs()-1);
	int dmin = width * (GetMyid() - 1);
	int dmax = (GetMyid() == GetNprocs() - 1) ? Nocc : width * GetMyid();

	k = -1;
	int cmin = -1;
	while ((cmin<Ncomponent) && (k<dmin))	{
		cmin++;
		if (cmin == Ncomponent)	{
			cerr << "error in matmixslavemoveprofile: overflow\n";
			exit(1);
		}
		if (occupancy[cmin])	{
			k++;
		}
	}
	int cmax = cmin;
	if (dmax == Nocc)	{
		cmax = Ncomponent;
	}
	else	{
		while ((cmax<Ncomponent) && (k<dmax))	{
			cmax++;
			if (occupancy[cmax])	{
				k++;
			}
		}
	}
	int nocc = 0;
	for (int j=cmin; j<cmax; j++)	{
		if (occupancy[j])	{
			nocc++;
		}
	}
	if (nocc != (dmax - dmin))	{
		cerr << "error : mismatch in nocc\n";
		exit(1);
	}

	// determine the range of components to move
	/*
	k = 0;
	int cmin = 0;
	while ((cmin<Ncomponent) && (k<dmin))	{
		cmin++;
		if (cmin == Ncomponent)	{
			cerr << "error in matmixslavemoveprofile: overflow\n";
			exit(1);
		}
		if (occupancy[cmin])	{
			k++;
		}
	}
	int cmax = cmin;
	if (dmax == Nocc)	{
		cmax = Ncomponent;
	}
	else	{
		while ((cmax<Ncomponent) && (k<dmax))	{
			cmax++;
			if (cmax == Ncomponent)	{
				cerr << "error in matmixslavemoveprofile: overflow\n";
				exit(1);
			}
			if (occupancy[cmax])	{
				k++;
			}
		}
	}

	if (GetMyid() == 15)	{
		cerr << cmin << '\t' << cmax << '\t' << dmin << '\t' << dmax << '\t' << Ncomponent << '\t' << GetNOccupiedComponent() << '\n';
		for (int i=cmin; i<cmax; i++)	{
			cerr << occupancy[i] << '\t';
		}
		cerr << '\n';
		for (int i=0; i<Ncomponent; i++)	{
			cerr << occupancy[i] << '\t';
		}
		cerr << '\n';
	}
	// cerr << GetMyid() << '\t' << cmin << '\t' << cmax << '\n';
	*/

	// update sufficient statistics
	UpdateModeProfileSuffStat();
	UpdateOccupancyNumbers();
	// move components in the range just computed
	double total = 0;
	for (int i=cmin; i<cmax; i++)	{
		if (occupancy[i])	{
			total += MoveProfile(i,tuning,n,nrep);
		}
		/*
		else	{
			SampleStat(i);
			total++;
		}
		*/
	}

	// send the new values of the profiles, plus the total success rate (total)
	double* tmp = new double[(dmax - dmin) * GetDim() + 1];
	int l = 0;
	for (int i=cmin; i<cmax; i++)	{
		if (occupancy[i])	{
			for (int k=0; k<GetDim(); k++)	{
				tmp[l] = profile[i][k];
				l++;
			}
		}
	}
	tmp[l] = total;
	MPI_Send(tmp,(dmax-dmin)*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	delete[] tmp;
}

double MatrixMixtureProfileProcess::MoveProfile(double tuning, int n, int nrep)	{
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += MoveProfile(i,tuning,n,nrep);
	}
	return total / GetNcomponent();
}

double MatrixMixtureProfileProcess::MoveProfile(int cat, double tuning, int n, int nrep)	{
// double MatrixMixtureProfileProcess::MoveProfile(int cat, double tuning, int n, int nrep, RngStream& g)	{
	int naccepted = 0;
	double* bk = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bk[k] = profile[cat][k];
	} 
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogStatPrior(cat) - ProfileSuffStatLogProb(cat);
		double loghastings = ProfileProposeMove(profile[cat],tuning,n,0);
		// double loghastings = ProfileProposeMove(profile[cat],tuning,n,0,g);
		UpdateComponent(cat);
		deltalogprob += LogStatPrior(cat) + ProfileSuffStatLogProb(cat);
		deltalogprob += loghastings;
		// int accepted = (g.RandU01() < exp(deltalogprob));
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted ++;
			for (int k=0; k<GetDim(); k++)	{
				bk[k] = profile[cat][k];
			} 
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[cat][k] = bk[k];
			} 
			UpdateComponent(cat);
		}
	}
	delete[] bk;
	return naccepted / nrep;
}

