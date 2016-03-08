
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef SEQUENCEALIGNMENT_H
#define SEQUENCEALIGNMENT_H

#include "StateSpace.h"
#include "TaxonSet.h"
#include "Random.h"

// this class works like an interface
// it does not do any job
class SequenceAlignment	{

	public:

	SequenceAlignment() : Data(0), BKData(0) {}

	SequenceAlignment(SequenceAlignment* from)	{

		Ntaxa = from->Ntaxa;
		Nsite = from->Nsite;
		taxset = from->taxset;
		statespace = from->statespace;

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = from->Data[i][j];
			}
		}
		BKData = 0;
	}
	
	SequenceAlignment(SequenceAlignment* from, int start, int length)	{

		Ntaxa = from->Ntaxa;
		if ((start < 0) || (length < 0) || (start + length > from->Nsite))	{
			cerr << "error in sequence alignment: overflow\n";
			exit(1);
		}

		Nsite = length;
		taxset = from->taxset;
		statespace = from->statespace;
		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = from->Data[i][j+start];
			}
		}
		BKData = 0;
	}
	
	SequenceAlignment(SequenceAlignment* from, int N, int Ngene, int* genesize, int* exclude)	{

		Ntaxa = from->Ntaxa;

		int* include = new int[Ngene];
		int nex = 0;
		for (int i=0; i<Ngene; i++)	{
			nex += exclude[i];
			include[i] = exclude[i];
		}

		int Nleft = Ngene - nex;
		for (int i=0; i<N; i++)	{
			int gene = (int) (Nleft * rnd::GetRandom().Uniform());
			int k = 0;
			while (gene)	{
				if (!include[k])	{
					gene--;
				}
				k++;
			}
			if (include[k])	{
				k++;
			}

			if (k >= Ngene)	{
				cerr << "error : Ngene overflow\n";
				exit(1);
			}

			include[k] = 1;
			Nleft--;
		}

		int Ncheck = 0;
		int totsize = 0;
		for (int i=0; i<Ngene; i++)	{
			include[i] -= exclude[i];
			if (include[i])	{
				Ncheck++;
				totsize += genesize[i];
			}
		}

		Nsite = totsize;

		taxset = from->taxset;
		statespace = from->statespace;

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}

		int kfrom = 0;
		int kto = 0;
		for (int i=0; i<Ngene; i++)	{
			if (include[i])	{
				for (int k=0; k<genesize[i]; k++)	{
					for (int j=0; j<Ntaxa; j++)	{
						Data[j][kto+k] = from->Data[j][kfrom+k];
					}
				}
				kto += genesize[i];
			}
			kfrom += genesize[i];
		}

		delete[] include;

		BKData = 0;
	}

	SequenceAlignment(SequenceAlignment* from, int Ngene, int* genesize, int* exclude, double* frac, double minfrac)	{

		cerr << "sub alignment, missing frac per gene less than " << minfrac << '\n';

		Ntaxa = from->Ntaxa;

		int Ninclude = 0;
		int totsize = 0;
		for (int i=0; i<Ngene; i++)	{
			if ((! exclude[i]) && (frac[i] <= minfrac))	{
				Ninclude++;
				totsize += genesize[i];
			}
		}

		cerr << "alignment size : " << totsize << '\n';

		Nsite = totsize;
		taxset = from->taxset;
		statespace = from->statespace;

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}

		int kfrom = 0;
		int kto = 0;
		for (int i=0; i<Ngene; i++)	{
			if ((! exclude[i]) && (frac[i] <= minfrac))	{
				for (int k=0; k<genesize[i]; k++)	{
					for (int j=0; j<Ntaxa; j++)	{
						Data[j][kto+k] = from->Data[j][kfrom+k];
					}
				}
				kto += genesize[i];
			}
			kfrom += genesize[i];
		}

		BKData = 0;
	}

	SequenceAlignment(SequenceAlignment* from, const TaxonSet* subset)	{

		Ntaxa = subset->GetNtaxa();
		Nsite = from->GetNsite();
		statespace = from->statespace;
		taxset = subset;
		Data = new int*[Ntaxa];
		for (int k=0; k<Ntaxa; k++)	{
			Data[k] = new int[Nsite];
		}
		BKData = 0;
		for (int i=0; i<from->GetNtaxa(); i++)	{
			int k  = subset->GetTaxonIndex(from->GetTaxonSet()->GetTaxon(i));
			if (k != -1)	{
				for (int j=0; j<Nsite; j++)	{
					Data[k][j] = from->Data[i][j];
				}
				if (k >= Ntaxa)	{
					cerr << "error in sequence alignment subset\n";
					cerr << Ntaxa << '\n';
					exit(1);
				}
			}
		}
	}

	SequenceAlignment(int** inData, string* names, int inNsite, StateSpace* instatespace, const TaxonSet* intaxset)	{

		Nsite = inNsite;
		taxset = intaxset;
		Ntaxa = taxset->GetNtaxa();
		statespace = instatespace;

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}
		for (int i=0; i<Ntaxa; i++)	{
			int mapi = taxset->GetTaxonIndex(names[i]);
			for (int j=0; j<Nsite; j++)	{
				Data[mapi][j] = inData[i][j];
			}
		}
	}

	void MissingFractionPerGene(int Ngene, int* genesize, int* exclude, double* frac)	{

		int offset = 0;

		for (int gene=0; gene<Ngene; gene++)	{
			

			int tot = 0;
			int mis = 0;
			for (int j=0; j<taxset->GetNtaxa(); j++)	{
				for (int i=0; i<genesize[gene]; i++)	{
					if (Data[j][offset+i] == unknown)	{
						mis++;
					}
					tot++;
				}
			}
			frac[gene] = ((double) mis) / tot;
			offset += genesize[gene];
		}
	}

	virtual ~SequenceAlignment() {

		if (Data)	{
			for (int k=0; k<Ntaxa; k++)	{
				delete[] Data[k];
			}
			delete[] Data;
			Data = 0;
		}
		if (BKData)	{
			for (int k=0; k<Ntaxa; k++)	{
				delete[] BKData[k];
			}
			delete[] BKData;
			BKData = 0;
		}

	}


	SequenceAlignment& operator=(SequenceAlignment& from)	{
		if (from.GetNsite() != GetNsite())	{
			cerr << "error in SequenceAlignment::operator=\n";
			exit(1);
		}
		if (from.GetNtaxa() != GetNtaxa())	{
			cerr << "error in SequenceAlignment::operator=\n";
			exit(1);
		}
		for (int i=0; i<Ntaxa; i++)	{
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = from.Data[i][j];
			}
		}

	}

	// the set of characters (A,C,G,T for nucleotides, etc..)
	StateSpace*  GetStateSpace()	{
		return statespace;
	}

	void Mask(SequenceAlignment* from)	{
		for (int i=0; i<from->GetNsite(); i++)	{
			for (int j=0; j<Ntaxa; j++)	{
				if (from->Data[j][i] == unknown)	{
					Data[j][i] = unknown;
				}
			}
		}
	}

	void Unclamp()	{
		if (! BKData)	{
			BKData = new int*[Ntaxa];
			for (int i=0; i<Ntaxa; i++)	{
				BKData[i] = new int[Nsite];
			}
		}
		for (int i=0; i<Ntaxa; i++)	{
			for (int j=0; j<Nsite; j++)	{
				BKData[i][j] = Data[i][j];
				Data[i][j] = unknown;
			}
		}
	}

	void Restore()	{
		if (! BKData)	{
			cerr << "error : cant restore data without backup\n";
			exit(1);
		}
		for (int i=0; i<Ntaxa; i++)	{
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = BKData[i][j];
			}
		}
	}

	int  GetNstate()	{
		return statespace->GetNstate();
	}
	
	// the list of taxa
	const TaxonSet* GetTaxonSet() const {
		return taxset;
	}

	int GetNsite()	{
		return Nsite;
	}

	int GetNtaxa()	{
		return taxset->GetNtaxa();
	}
	
	bool isMissing(int taxon, int site)	{
		return Data[taxon][site] == -1;
	}

	bool NoMissingColumn(int site)	{
		bool ret = true;
		int tax = 0;
		while ((tax < GetNtaxa()) && ret)	{
			ret &= (Data[tax][site] != unknown);
			tax++;
		}
		return ret;
	}
	
	bool ConstantColumn(int site)	{
		bool ret = true;
		int tax = 0;
		while ((tax < GetNtaxa()) && (Data[tax][site] == unknown))	{
			tax++;
		}

		if (tax < GetNtaxa())	{
			int refstate = Data[tax][site];

			while ((tax < GetNtaxa()) && ret)	{
				if (Data[tax][site] != -1)	{
					ret &= (Data[tax][site] == refstate);
				}
				tax++;
			}
		}
		return ret;
	}

	void SetState(int taxon, int site, int state)	{
		Data[taxon][site] = state;
	}

	int GetState(int taxon, int site)	{
		return Data[taxon][site];
	}

	int GetBKState(int taxon, int site)	{
		if (! BKData)	{
			cerr << "error in get bk state\n";
			exit(1);
		}
		return BKData[taxon][site];
	}

	const int* GetState(int taxon)	{
		return Data[taxon];
	}

	/*
	const int** GetDataMatrix()	{
		return data;
	}
	*/

	void GetEmpiricalFreq(double* in);

	void GetSiteEmpiricalFreq(double** in);

	void ToStream(ostream& os);
	void ToFasta(ostream& os);

	void SetTestData(int testnsite, int offset, int sitemin, int sitemax, int* tmp)	{

		int index = 0;
		int tmpdata[Ntaxa][testnsite];
		for (int k=0; k<Ntaxa; k++)	{
			for (int i=0; i<testnsite; i++)	{
				tmpdata[k][i] = tmp[index];
				index++;
			}
		}
		
		for (int k=0; k<Ntaxa; k++)	{
			for (int i=sitemin; i<sitemax; i++)	{
				Data[k][offset + (i-sitemin)] = tmpdata[k][i];
			}
		}
	}

	void GetDataVector(int* tmp)	{
		int index = 0;
		for (int k=0; k<Ntaxa; k++)	{
			for (int i=0; i<Nsite; i++)	{
				tmp[index] = Data[k][i];
				index++;
			}
		}
	}

	void DeleteConstantSites()	{
		int i=0;
		int j=0;
		int Eliminated = 0;
		while (i<Nsite)	{
			int k = 0;
			while ((k<Ntaxa) && (Data[k][i] == unknown)) k++;
			if (k<Ntaxa)	{
				int a = Data[k][i];
				k++;
				while ((k<Ntaxa) && ((Data[k][i] == unknown) || (Data[k][i] == a))) k++;
				if (k==Ntaxa)	{
					Eliminated ++;
				}
				else	{
					for (int k=0; k<Ntaxa; k++)	{
						Data[k][j] = Data[k][i];
					}
					j++;
				}
			}
			i++;
		}
		Nsite -= Eliminated;
		cout << "number of positions eliminated : " << Eliminated << '\n';
	}

	virtual double GetTotalDiversity(int sitemin, int sitemax)	{
		double total = 0;
		int obs[GetNstate()];
		for (int i=sitemin; i<sitemax; i++)	{
			for (int k=0; k<GetNstate(); k++)	{
				obs[k] = 0;
			}
			for (int j=0; j<Ntaxa; j++)	{
				int state = GetState(j,i);
				if (state != unknown)	{
					obs[state]++;
				}
			}
			int div = 0;
			for (int k=0; k<GetNstate(); k++)	{
				if (obs[k])	{
					div++;
				}
			}
			total += div;
		}
		// total /= (sitemax - sitemin);
		return total;
	}

	double GetMeanDiversity()	{
		return GetTotalDiversity(0,GetNsite()) / GetNsite();
	}

	double CompositionalHeterogeneity(ostream* os)	{

		int Nstate = GetNstate();
		double** taxfreq = new double*[Ntaxa];
		for (int j=0; j<Ntaxa; j++)	{
			taxfreq[j] = new double[Nstate];
			for (int k=0; k<Nstate; k++)	{
				taxfreq[j][k] = 0;
			} 
		}

		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Ntaxa; j++)	{
				int state = GetState(j,i);
				if (state != unknown)	{
					taxfreq[j][state]++;
				}
			}
		}
				
		for (int j=0; j<Ntaxa; j++)	{
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += taxfreq[j][k];
			}
			for (int k=0; k<Nstate; k++)	{
				taxfreq[j][k] /= total;
			}
		}

		// make global freqs out of tax-specific freqs
		double* globalfreq = new double[Nstate];
		for (int k=0; k<Nstate; k++)	{
			globalfreq[k] = 0;
			for (int j=0; j<Ntaxa; j++)	{
				globalfreq[k] += taxfreq[j][k];
			}
			globalfreq[k] /= Ntaxa;
		} 

		for (int j=0; j<Ntaxa; j++)	{
			for (int k=0; k<Nstate; k++)	{
				if (os)	{
					(*os) << taxfreq[j][k] << '\t';
				}
			}
			if (os)	{
				(*os) << '\n';
			}
		}
		if (os)	{
			(*os) << '\n';
		}

		// compute max distance
		double maxdist = 0;
		for (int j=0; j<Ntaxa; j++)	{
			double dist = 0;
			for (int k=0; k<Nstate; k++)	{
				double tmp = (taxfreq[j][k] - globalfreq[k]);
				dist += tmp * tmp;
			}
			if (maxdist < dist)	{
				maxdist = dist;
			}
		}

		delete[] globalfreq;
		for (int j=0; j<Ntaxa; j++)	{
			delete[] taxfreq[j];
		}
		delete[] taxfreq;

		return maxdist;
	}

	// data fields

	int Ntaxa;
	int Nsite;
	const TaxonSet* taxset;
	StateSpace* statespace;
	int** Data;
	int** BKData;
	
};

class FileSequenceAlignment : public SequenceAlignment	{


	public:
		FileSequenceAlignment(istream& is);
		FileSequenceAlignment(string filename,int fullline,int myid);

	private:

	int 			ReadDataFromFile(string filename, int forceinterleaved = 0);
	int 			ReadNexus(string filename);
	int 			TestPhylipSequential(string filename);
	void 			ReadPhylipSequential(string filename);
	int 			TestPhylip(string filename, int repeattaxa);
	void 			ReadPhylip(string filename, int repeattaxa);
	int			ReadSpecial(string filename);

	string* SpeciesNames;
};



#endif // SEQUENCEALIGNMENT_H
