
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
#include "AAMutSelSBDPPhyloProcess.h"
#include "Parallel.h"
#include <string.h>


// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
// are model dependent
// should be implemented in .cpp file
void AAMutSelSBDPPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nnucrr = GetNnucrr(),nnucstat = 4;
	L1 = GetNmodeMax();
	L2 = GetDim();
	//nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1); // check if these last terms are correct in this context...
	nd = 2 + nbranch + nnucrr + nnucstat + L1*L2 + GetDim() + 1;
	ni = 1 + ProfileProcess::GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int index = 0;
	branchalpha = dvector[index];
	index++;
	branchbeta = dvector[index];
	index++;
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
	for(i=0; i<nnucrr; ++i) {
		nucrr[i] = dvector[index];
		index++;
	}
	for(i=0; i<nnucstat; ++i) {
		nucstat[i] = dvector[index];
		index++;
	}
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}
	kappa = dvector[index];
	index++;
	
	Ncomponent = ivector[0];
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		SBDPProfileProcess::alloc[i] = ivector[1+i];
	}
	//GetBranchLengthsFromArray();
	delete[] dvector;
	delete[] ivector;
	// this one is really important
	// in those cases where new components have appeared, or some old ones have disappeared
	// during allocation move on the master node.
	// 
	// note that CreateMatrices() in fact creates only those that are not yet allocated
	// and also deletes those that are now obsolete
	CreateMatrices();
	UpdateMatrices();

	/*
	SlaveBroadcastTree();

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nnucrr = GetNnucrr(),nnucstat = 4,k = 0;
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1); // check if these last terms are correct in this context...
	ni = 1 + ProfileProcess::GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[i];
	}
	for(i=0; i<nnucrr; ++i) {
		nucrr[i] = dvector[nbranch+i];
	}
	for(i=0; i<nnucstat; ++i) {
		nucstat[i] = dvector[nbranch+nnucrr+i];
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[1+nbranch+nnucrr+nnucstat+i];
	}
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[nbranch+nnucrr+nnucstat+L2+k];
			k++;
		}
		weight[i] = dvector[nbranch+nnucrr+nnucstat+L2+k];
		k++;
	}
	Ncomponent = ivector[0];
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		SBDPProfileProcess::alloc[i] = ivector[1+i];
	}
	//GetBranchLengthsFromArray();
	delete[] dvector;
	delete[] ivector;
	// this one is really important
	// in those cases where new components have appeared, or some old ones have disappeared
	// during allocation move on the master node.
	// 
	// note that CreateMatrices() in fact creates only those that are not yet allocated
	// and also deletes those that are now obsolete
	CreateMatrices();
	UpdateMatrices();
	*/
}


void AAMutSelSBDPPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

	/*
	case PRINT_TREE:
		SlavePrintTree();
		break;
	*/
	case REALLOC_MOVE:
		SlaveIncrementalDPMove();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	case MIX_MOVE:
		SlaveMixMove();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void AAMutSelSBDPPhyloProcess::GlobalUpdateParameters() {
	// MPI2
	// should send the slaves the relevant information
	// about model parameters
	// for this model, should broadcast
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	assert(myid == 0);
	int i,j,nnucrr,nnucstat,nbranch = GetNbranch(),ni,nd,L1,L2;
	nnucrr = GetNnucrr();
	nnucstat = 4;	
	L1 = GetNmodeMax();
	L2 = GetDim();
	//nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1);  // check if these last terms are correct in this context...
	nd = 2 + nbranch + nnucrr + + nnucstat + L1*L2 + GetDim() + 1;
	ni = 1 + ProfileProcess::GetNsite(); // 1 for the number of componenets, and the rest for allocations
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	// GlobalBroadcastTree();
	// First we assemble the vector of doubles for distribution
	int index = 0;
	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	index++;
	
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}

	for(i=0; i<nnucrr; ++i) {
		dvector[index] = nucrr[i];
		index++;
	}
	for(i=0; i<nnucstat; ++i) {
		dvector[index] = nucstat[i];
		index++;
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}
	dvector[index] = kappa;
	index++;

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		ivector[1+i] = SBDPProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	/*
	// MPI2
	// should send the slaves the relevant information
	// about model parameters
	// for this model, should broadcast
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	assert(myid == 0);
	int i,j,nnucrr,nnucstat,nbranch = GetNbranch(),ni,nd,L1,L2,k = 0;
	nnucrr = GetNnucrr();
	nnucstat = 4;	
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1);  // check if these last terms are correct in this context...
	ni = 1 + ProfileProcess::GetNsite(); // 1 for the number of componenets, and the rest for allocations
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	GlobalBroadcastTree();
	// First we assemble the vector of doubles for distribution
	for(i=0; i<nbranch; ++i) {
		dvector[i] = blarray[i];
	}

	for(i=0; i<nnucrr; ++i) {
		dvector[nbranch+i] = nucrr[i];
	}
	for(i=0; i<nnucstat; ++i) {
		dvector[nbranch+nnucrr+i] = nucstat[i];
	}

	for (int i=0; i<GetDim(); i++)	{
		dvector[nbranch+nnucrr+nnucstat+i] = dirweight[i];
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[nbranch+nnucrr+nnucstat+L2+k] = profile[i][j];
			k++;
		}
		dvector[nbranch+nnucrr+nnucstat+L2+k] = weight[i];
		k++;
	}

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		ivector[1+i] = SBDPProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/
}

void AAMutSelSBDPPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = 0;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	int sel = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-sel")	{
				sel = 1;
			}
			//else if (s == "-ppred")	{
			//	ppred = 1;
			//}
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

	//if (sel)	{
	//	ReadSDistributions(name,burnin,every,until);
	//}
	//else if (ppred)	{
	//	PostPred(ppred,name,burnin,every,until);
	//}
	//else	{
		Read(name,burnin,every,until);
	//}
}

void AAMutSelSBDPPhyloProcess::Read(string name, int burnin, int every, int until)	{

	cerr << "in read\n";
	cerr.flush();

	ifstream is((name + ".chain").c_str());
	if (! is)	{
		cerr << "error: did not find " << name << ".chain\n";
		exit(1);
	}
	//cerr << "In AAMutSelDPPhyloProcess. GetDim() is : " << GetDim() << "\n";
	int Nstate = AAMutSelSBDPSubstitutionProcess::GetNstate();
	//cerr << "Nstate is: " << Nstate << "\n";
	//cerr.flush();
	double TOOSMALL = 1e-10;

	double meanlength = 0;
	double varlength = 0;
	int Ncat = 241;
	double min = -30;
	double max = 30;
	double step = 0.25;
	double** PosteriorMeanSiteAAP = new double*[ProfileProcess::GetNsite()];
	double** sshistoMut = new double*[ProfileProcess::GetNsite()];
	double** sshistoSub = new double*[ProfileProcess::GetNsite()];
	double** sshistoNonsynMut = new double*[ProfileProcess::GetNsite()];
	double** sshistoNonsynSub = new double*[ProfileProcess::GetNsite()];
	for (int site = 0; site < ProfileProcess::GetNsite(); site++)        {
		PosteriorMeanSiteAAP[site] = new double[GetDim()];
		sshistoMut[site] = new double[Ncat];
		sshistoSub[site] = new double[Ncat];
		sshistoNonsynMut[site] = new double[Ncat];
		sshistoNonsynSub[site] = new double[Ncat];
		for (int a = 0; a < GetDim(); a++)   {
			PosteriorMeanSiteAAP[site][a] = 0;
		}
		for (int a = 0; a < Ncat; a++)	{
			sshistoMut[site][a] = 0;
			sshistoSub[site][a] = 0;
			sshistoNonsynMut[site][a] = 0;
			sshistoNonsynSub[site][a] = 0;
		}
	}
	double* meanNucStat = new double[Nnuc];
	for (int i=0; i<Nnuc; i++)	{
		meanNucStat[i]=0.0;
	}
	double* meanNucRR = new double[Nnucrr];
	for (int i=0; i<Nnucrr; i++)	{
		meanNucRR[i] = 0.0;
	}
	
	double* ghistoMut = new double[Ncat];
	double* ghistoSub = new double[Ncat];
	double* ghistoNonsynMut = new double[Ncat];
	double* ghistoNonsynSub = new double[Ncat];

	double* shistoMut = new double[Ncat];
	double* shistoSub = new double[Ncat];
	double* shistoNonsynMut = new double[Ncat];
	double* shistoNonsynSub = new double[Ncat];
	
	double* tsshistoMut = new double[Ncat];
	double* tsshistoSub = new double[Ncat];
	double* tsshistoNonsynMut = new double[Ncat];
	double* tsshistoNonsynSub = new double[Ncat];
	for (int c = 0; c < Ncat; c++)	{
		ghistoMut[c] = 0;
		ghistoSub[c] = 0;
		ghistoNonsynMut[c] = 0;
		ghistoNonsynSub[c] = 0;
	}
	//const double* stat;
	double* stat = new double[Nstate];
	int pos, nucFrom, nucTo, nucRRIndex, aaFrom, aaTo, c;
	double statMutRate, S, statSubRate, totalMut, totalSub, totalNonsynMut, totalNonsynSub, siteTotalMut, siteTotalSub, siteTotalNonsynMut, siteTotalNonsynSub, Z;
	cerr << "Ncat is " << Ncat << "\n";
	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		cerr << ".";
		cerr.flush();
		FromStream(is);
		i++;
	}
	cerr << "\nburnin complete\n";
	cerr.flush();
	int samplesize = 0;
	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		//UpdateMatrices();
		//Trace(cerr);
		double length = GetTotalLength();
		// int nocc = process->GetNOccupiedComponent();
		meanlength += length;
		varlength += length * length;
		for (int a=0; a<Nnuc; a++)	{
			meanNucStat[a] += GetNucStat(a);
		}
		for (int a=0; a<Nnucrr; a++)	{
			meanNucRR[a] += GetNucRR(a);
		}
		totalMut = 0;
		totalSub = 0;
		totalNonsynMut = 0;
		totalNonsynSub = 0;
		for (c = 0; c < Ncat; c++)	{
			shistoMut[c] = 0;
			shistoSub[c] = 0;
			shistoNonsynMut[c] = 0;
			shistoNonsynSub[c] = 0;
		}
		for (int site=0; site<ProfileProcess::GetNsite(); site++)	{
			siteTotalMut = 0;
			siteTotalSub = 0;
			siteTotalNonsynMut = 0;
			siteTotalNonsynSub = 0;
			
			for (c = 0; c < Ncat; c++)	{
				tsshistoMut[c] = 0;
				tsshistoSub[c] = 0;
				tsshistoNonsynMut[c] = 0;
				tsshistoNonsynSub[c] = 0;
			}
			for (int a=0; a<GetDim(); a++)	{
				PosteriorMeanSiteAAP[site][a] += profile[alloc[site]][a]; 
			}
			Z = 0;
			for (int s=0; s<Nstate; s++)	{
				stat[s] = 	GetNucStat(AAMutSelProfileProcess::statespace->GetCodonPosition(0, s)) *
						GetNucStat(AAMutSelProfileProcess::statespace->GetCodonPosition(1, s)) *
						GetNucStat(AAMutSelProfileProcess::statespace->GetCodonPosition(2, s)) *
						profile[alloc[site]][AAMutSelProfileProcess::statespace->Translation(s)];
				Z += stat[s];
			}
			for (int s=0; s<Nstate; s++)	{
				stat[s] /= Z;
			}
			//stat = GetStationary(site);
			//stat = matrixarray[alloc[site]]->GetStationary();
			for (int codonFrom = 0; codonFrom<Nstate; codonFrom++)	{
				for (int codonTo = 0; codonTo<Nstate; codonTo++)	{
					pos = AAMutSelProfileProcess::statespace->GetDifferingPosition(codonFrom, codonTo);				
					if ((pos != -1) && (pos != 3))  {
						nucFrom = AAMutSelProfileProcess::statespace->GetCodonPosition(pos, codonFrom);
						nucTo = AAMutSelProfileProcess::statespace->GetCodonPosition(pos, codonTo);
						if (nucFrom<nucTo)	{
							nucRRIndex = (2 * Nnuc - nucFrom - 1) * nucFrom / 2 + nucTo - nucFrom - 1;
						}
						else {
							nucRRIndex = (2 * Nnuc - nucTo - 1) * nucTo / 2 + nucFrom - nucTo - 1;
						}
						statMutRate = GetNucRR(nucRRIndex) * GetNucStat(nucTo) * stat[codonFrom];
						if (!AAMutSelProfileProcess::statespace->Synonymous(codonFrom, codonTo))	{
							aaFrom = AAMutSelProfileProcess::statespace->Translation(codonFrom);
							aaTo = AAMutSelProfileProcess::statespace->Translation(codonTo);
							S = log(profile[alloc[site]][aaTo]/profile[alloc[site]][aaFrom]);

							if (fabs(S) < TOOSMALL)	{
								statSubRate = statMutRate * 1.0/( 1.0 - (S / 2) );
							}
							else {
								statSubRate = statMutRate * (S/(1-(profile[alloc[site]][aaFrom]/profile[alloc[site]][aaTo])));
							}
							if (S < min)	{
								c = 0;
							}
							else if (S > max)	{
								c = Ncat-1;
							}
							else {
								c = 0;
								double tmp = min + ((double)c * step) - step/2 + step;
								do	{
									c++;
									tmp = min + ((double)(c) * step) - step/2 + step;
								} while (tmp < S );
							}
							if (c == Ncat)	{
								cout << "error, c==Ncat.\n";
								cout.flush();
							}
	
							//while ((S > (min + c * step)) && (c<Ncat)) c++;
							//while (S > (min + c * step)) c++;
							shistoMut[c] += statMutRate;
							shistoSub[c] += statSubRate;
							shistoNonsynMut[c] += statMutRate;
							shistoNonsynSub[c] += statSubRate;

							tsshistoMut[c] += statMutRate;
							tsshistoSub[c] += statSubRate;
							tsshistoNonsynMut[c] += statMutRate;
							tsshistoNonsynSub[c] += statSubRate;

							totalNonsynMut += statMutRate;
							totalNonsynSub += statSubRate;
							siteTotalNonsynMut += statMutRate;
							siteTotalNonsynSub += statSubRate;
						}
						else {
							S = 0.0;
							statSubRate = statMutRate;
							shistoMut[(int)(Ncat/2)] += statMutRate; // change this to make more flexible
							shistoSub[(int)(Ncat/2)] += statMutRate; // change this to make more flexible
							tsshistoMut[(int)(Ncat/2)] += statMutRate;
							tsshistoSub[(int)(Ncat/2)] += statMutRate;
						}
						totalMut += statMutRate;
						totalSub += statSubRate;
						siteTotalMut += statMutRate;
						siteTotalSub += statSubRate;
						//cerr << "nucFrom is: " << nucFrom << "\n";
						//cerr << "mutRate is: " << statMutRate << "\n";
						//cerr << "subRate is: " << statSubRate << "\n";
					}

				}
				//cerr << "codonStat[" << codonFrom << "]: " << stat[codonFrom];
				//cerr << ", amino acid " << AAMutSelProfileProcess::statespace->Translation(codonFrom);
				//cerr << ", nuc pos 1: " << AAMutSelProfileProcess::statespace->GetCodonPosition(0, codonFrom);
				//cerr << ", nuc pos 2: " << AAMutSelProfileProcess::statespace->GetCodonPosition(1, codonFrom);
				//cerr << ", nuc pos 3: " << AAMutSelProfileProcess::statespace->GetCodonPosition(2, codonFrom);
				//cerr << "\n";
			}

			for (c=0; c<Ncat; c++)	{
				sshistoMut[site][c] += tsshistoMut[c]/siteTotalMut;
				sshistoSub[site][c] += tsshistoSub[c]/siteTotalSub;
				sshistoNonsynMut[site][c] += tsshistoNonsynMut[c]/siteTotalNonsynMut;
				sshistoNonsynSub[site][c] += tsshistoNonsynSub[c]/siteTotalNonsynSub;
			}

			//exit(1);
		}	
		// cerr << process->GetLogLikelihood() << '\t' << logl << '\t' << length << '\n';

		for (c=0; c<Ncat; c++)	{
			ghistoMut[c] += shistoMut[c]/totalMut;
			ghistoSub[c] += shistoSub[c]/totalSub;
			ghistoNonsynMut[c] += shistoNonsynMut[c]/totalNonsynMut;
			ghistoNonsynSub[c] += shistoNonsynSub[c]/totalNonsynSub;
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	meanlength /= samplesize;
	varlength /= samplesize;
	varlength -= meanlength * meanlength;

	cerr << "mean length : " << meanlength << " +/- " << sqrt(varlength) << '\n';

	double* ppvhistoMut = new double[Ncat];
	double* ppvhistoSub = new double[Ncat];
	double* ppvhistoNonsynMut = new double[Ncat];
	double* ppvhistoNonsynSub = new double[Ncat];
	for (c=0; c<Ncat; c++)	{
		ppvhistoMut[c] = 0;
		ppvhistoSub[c] = 0;
		ppvhistoNonsynMut[c] = 0;
		ppvhistoNonsynSub[c] = 0;
	}
	
	totalMut = 0;
	totalSub = 0;
	totalNonsynMut = 0;
	totalNonsynSub = 0;

	cerr << "final treatment\n";

	ofstream aap_os( (name + ".aapforlogo").c_str(), std::ios::out);
	for (int site=0; site<ProfileProcess::GetNsite(); site++)	{
		//cerr << ".";
		cerr.flush();
		for (int a=0; a<GetDim(); a++)	{
			aap_os << PosteriorMeanSiteAAP[site][a] / samplesize << "\t";
		}
		aap_os << "\n";

		Z = 0;
		for (int s=0; s<Nstate; s++)	{
			stat[s] = 	(meanNucStat[AAMutSelProfileProcess::statespace->GetCodonPosition(0, s)] *
					meanNucStat[AAMutSelProfileProcess::statespace->GetCodonPosition(1, s)] *
					meanNucStat[AAMutSelProfileProcess::statespace->GetCodonPosition(2, s)] *
					PosteriorMeanSiteAAP[site][AAMutSelProfileProcess::statespace->Translation(s)]) / samplesize;
			Z += stat[s];
		}
		for (int s=0; s<Nstate; s++)	{
			stat[s] /= Z;
		}
		for (int codonFrom = 0; codonFrom<Nstate; codonFrom++)	{
			for (int codonTo = 0; codonTo<Nstate; codonTo++)	{
				pos = AAMutSelProfileProcess::statespace->GetDifferingPosition(codonFrom, codonTo);				
				if ((pos != -1) && (pos != 3))  {
					nucFrom = AAMutSelProfileProcess::statespace->GetCodonPosition(pos, codonFrom);
					nucTo = AAMutSelProfileProcess::statespace->GetCodonPosition(pos, codonTo);
					if (nucFrom<nucTo)	{
						nucRRIndex = (2 * Nnuc - nucFrom - 1) * nucFrom / 2 + nucTo - nucFrom - 1;
					}
					else {
						nucRRIndex = (2 * Nnuc - nucTo - 1) * nucTo / 2 + nucFrom - nucTo - 1;
					}
					statMutRate = meanNucRR[nucRRIndex] * meanNucStat[nucTo] * stat[codonFrom];
					if (!AAMutSelProfileProcess::statespace->Synonymous(codonFrom, codonTo))	{
						aaFrom = AAMutSelProfileProcess::statespace->Translation(codonFrom);
						aaTo = AAMutSelProfileProcess::statespace->Translation(codonTo);
						double expS = (PosteriorMeanSiteAAP[site][aaTo]/samplesize)/(PosteriorMeanSiteAAP[site][aaFrom]/samplesize);
						//S = log((PosteriorMeanAAP[site][aaTo]/samplesize)/(PosteriorMeanAAP[site][aaFrom]/samplesize));
						S = log(expS);
						statSubRate = statMutRate * (S/(1-((PosteriorMeanSiteAAP[site][aaFrom]/samplesize)/(PosteriorMeanSiteAAP[site][aaTo]/samplesize))));
						if (S < min)	{
							c = 0;
						}
						else if (S > max)	{
							c = Ncat-1;
						}
						else {
							c = 0;
							double tmp = min + ((double)c * step) - step/2 + step;
							do	{
								c++;
								tmp = min + ((double)(c) * step) - step/2 + step;
							} while (tmp < S );
						}
						if (c == Ncat)	{
							cout << "error, c==Ncat.\n";
							cout.flush();
						}

						ppvhistoMut[c] += statMutRate;
						ppvhistoSub[c] += statSubRate;
						ppvhistoNonsynMut[c] += statMutRate;
						ppvhistoNonsynSub[c] += statSubRate;
						totalNonsynMut += statMutRate;
						totalNonsynSub += statSubRate;
					}
					else {
						S = 0.0;
						statSubRate = statMutRate;
						ppvhistoMut[(int)(Ncat/2)] += statMutRate;
						ppvhistoSub[(int)(Ncat/2)] += statMutRate;
					}
					totalMut += statMutRate;
					totalSub += statSubRate;
				}
			}
		}
	}	

	
	ofstream mutmutsel_os( (name + ".mutsel").c_str(), std::ios::out);
	ofstream mutsubsel_os( (name + ".subsel").c_str(), std::ios::out);
	ofstream nonsynmutmutsel_os( (name + ".nonsynmutsel").c_str(), std::ios::out);
	ofstream nonsynmutsubsel_os( (name + ".nonsynsubsel").c_str(), std::ios::out);
	ofstream sitemutmutsel_os( (name + ".sitemutsel").c_str(), std::ios::out);
	ofstream sitemutsubsel_os( (name + ".sitesubsel").c_str(), std::ios::out);
	ofstream sitenonsynmutmutsel_os( (name + ".sitenonsynmutsel").c_str(), std::ios::out);
	ofstream sitenonsynmutsubsel_os( (name + ".sitenonsynsubsel").c_str(), std::ios::out);
	ofstream ppvmutmutsel_os( (name + ".ppvmutsel").c_str(), std::ios::out);
	ofstream ppvmutsubsel_os( (name + ".ppvsubsel").c_str(), std::ios::out);
	ofstream ppvnonsynmutmutsel_os( (name + ".ppvnonsynmutsel").c_str(), std::ios::out);
	ofstream ppvnonsynmutsubsel_os( (name + ".ppvnonsynsubsel").c_str(), std::ios::out);

	for (c = 0; c < Ncat; c++)	{
		mutmutsel_os << (min + (c * step)) << "\t" << (ghistoMut[c]/samplesize) << '\n';
		mutsubsel_os << (min + (c * step)) << "\t" << (ghistoSub[c]/samplesize) << '\n';
		nonsynmutmutsel_os << (min + (c * step)) << "\t" << (ghistoNonsynMut[c]/samplesize) << '\n';
		nonsynmutsubsel_os << (min + (c * step)) << "\t" << (ghistoNonsynSub[c]/samplesize) << '\n';

		ppvmutmutsel_os << (min + (c * step)) << "\t" << (ppvhistoMut[c]/totalMut) << '\n';
		ppvmutsubsel_os << (min + (c * step)) << "\t" << (ppvhistoSub[c]/totalSub) << '\n';
		ppvnonsynmutmutsel_os << (min + (c * step)) << "\t" << (ppvhistoNonsynMut[c]/totalNonsynMut) << '\n';
		ppvnonsynmutsubsel_os << (min + (c * step)) << "\t" << (ppvhistoNonsynSub[c]/totalNonsynSub) << '\n';
		for (int site = 0; site < ProfileProcess::GetNsite(); site++)	{
			if (site == 0) 	{
				sitemutmutsel_os << (min + (c * step)) << '\t' << (sshistoMut[site][c]/samplesize) << '\t';
				sitemutsubsel_os << (min + (c * step)) << '\t' << (sshistoSub[site][c]/samplesize) << '\t';
				sitenonsynmutmutsel_os << (min + (c * step)) << '\t' << (sshistoNonsynMut[site][c]/samplesize) << '\t';
				sitenonsynmutsubsel_os << (min + (c * step)) << '\t' << (sshistoNonsynSub[site][c]/samplesize) << '\t';
			}
			else if (site == ProfileProcess::GetNsite() - 1)	{
				sitemutmutsel_os << (sshistoMut[site][c]/samplesize) << '\n';
				sitemutsubsel_os << (sshistoSub[site][c]/samplesize) << '\n';
				sitenonsynmutmutsel_os << (sshistoNonsynMut[site][c]/samplesize) << '\n';
				sitenonsynmutsubsel_os << (sshistoNonsynSub[site][c]/samplesize) << '\n';
			}
			else {
				sitemutmutsel_os << (sshistoMut[site][c]/samplesize) << '\t';
				sitemutsubsel_os << (sshistoSub[site][c]/samplesize) << '\t';
				sitenonsynmutmutsel_os << (sshistoNonsynMut[site][c]/samplesize) << '\t';
				sitenonsynmutsubsel_os << (sshistoNonsynSub[site][c]/samplesize) << '\t';

			}
		}
	}

	cerr << samplesize << '\n';
	for (int site = 0; site < GetNmodeMax(); site++)        {
		delete[] PosteriorMeanSiteAAP[site];
		delete[] sshistoMut[site];
		delete[] sshistoSub[site];
		delete[] sshistoNonsynMut[site];
		delete[] sshistoNonsynSub[site];
	}
	delete[] PosteriorMeanSiteAAP;
	delete[] sshistoMut;
	delete[] sshistoSub;
	delete[] sshistoNonsynMut;
	delete[] sshistoNonsynSub;
	delete[] tsshistoMut;
	delete[] tsshistoSub;
	delete[] tsshistoNonsynMut;
	delete[] tsshistoNonsynSub;
	delete[] ghistoMut;
	delete[] ghistoMut;
	delete[] ghistoSub;
	delete[] ghistoNonsynMut;
	delete[] ghistoNonsynSub;
	delete[] shistoMut;
	delete[] shistoSub;
	delete[] shistoNonsynMut;
	delete[] shistoNonsynSub;
	delete[] ppvhistoMut;
	delete[] ppvhistoSub;
	delete[] ppvhistoNonsynMut;
	delete[] ppvhistoNonsynSub;
	delete[] stat;
	delete[] meanNucStat;
	delete[] meanNucRR;
}
