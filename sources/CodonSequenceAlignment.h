
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef CODONSEQUENCEALIGNMENT_H
#define CODONSEQUENCEALIGNMENT_H


//#include "ContinuousData.h"
#include "SequenceAlignment.h"
#include "CodonStateSpace.h"
#include <cmath>

class CodonSequenceAlignment : public SequenceAlignment	{

	public:

	CodonSequenceAlignment(CodonSequenceAlignment* from) : SequenceAlignment((SequenceAlignment*) from) {}

	CodonSequenceAlignment(SequenceAlignment* from, bool force_stops = false,GeneticCodeType type = Universal);

	~CodonSequenceAlignment() {}

	void DeleteAAConstantSites()	{
		int i=0;
		int j=0;
		int Eliminated = 0;
		while (i<Nsite)	{
			int k = 0;
			while ((k<Ntaxa) && (Data[k][i] == unknown)) k++;
			if (k<Ntaxa)	{
				int a = GetCodonStateSpace()->Translation(Data[k][i]);
				k++;
				while ((k<Ntaxa) && ((Data[k][i] == unknown) || (GetCodonStateSpace()->Translation(Data[k][i]) == a))) k++;
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
		int obs[Naa];
		for (int i=sitemin; i<sitemax; i++)	{
			for (int k=0; k<Naa; k++)	{
				obs[k] = 0;
			}
			for (int j=0; j<Ntaxa; j++)	{
				int state = GetState(j,i);
				if (state != unknown)	{
					int aa = GetCodonStateSpace()->Translation(state);
					obs[aa]++;
				}
			}
			int div = 0;
			for (int k=0; k<Naa; k++)	{
				if (obs[k])	{
					div++;
				}
			}
			total += div;
		}
		return total;
	}

	double Nucleotide123CompositionalHeterogeneity()	{

		double** taxfreq12 = new double*[Ntaxa];
		double** taxfreq3 = new double*[Ntaxa];
		for (int j=0; j<Ntaxa; j++)	{
			taxfreq12[j] = new double[Nnuc];
			taxfreq3[j] = new double[Nnuc];
			for (int k=0; k<Nnuc; k++)	{
				taxfreq12[j][k] = 0;
				taxfreq3[j][k] = 0;
			} 
		}

		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Ntaxa; j++)	{
				int state = GetState(j,i);
				if (state != unknown)	{
					taxfreq12[j][GetCodonStateSpace()->GetCodonPosition(0,state)]++;
					taxfreq12[j][GetCodonStateSpace()->GetCodonPosition(1,state)]++;
					taxfreq3[j][GetCodonStateSpace()->GetCodonPosition(2,state)]++;
				}
			}
		}
				
		for (int j=0; j<Ntaxa; j++)	{
			double total = 0;
			for (int k=0; k<Nnuc; k++)	{
				total += taxfreq12[j][k];
			}
			for (int k=0; k<Nnuc; k++)	{
				taxfreq12[j][k] /= total;
			}
		}

		for (int j=0; j<Ntaxa; j++)	{
			double total = 0;
			for (int k=0; k<Nnuc; k++)	{
				total += taxfreq3[j][k];
			}
			for (int k=0; k<Nnuc; k++)	{
				taxfreq3[j][k] /= total;
			}
		}

		// compute max distance
		double maxdist = 0;
		for (int j=0; j<Ntaxa; j++)	{
			double dist = 0;
			for (int k=0; k<Nnuc; k++)	{
				double tmp = (taxfreq12[j][k] - taxfreq3[j][k]);
				dist += tmp * tmp;
			}
			if (maxdist < dist)	{
				maxdist = dist;
			}
		}

		for (int j=0; j<Ntaxa; j++)	{
			delete[] taxfreq12[j];
			delete[] taxfreq3[j];
		}
		delete[] taxfreq12;
		delete[] taxfreq3;

		return maxdist;
	}

	double NucleotideCompositionalHeterogeneity(ostream* os, int pos = -1, double** comp = 0)	{

		double** taxfreq = 0;
		if (comp)	{
			taxfreq = comp;
		}
		else	{
			taxfreq = new double*[Ntaxa];
			for (int j=0; j<Ntaxa; j++)	{
				taxfreq[j] = new double[Nnuc];
				for (int k=0; k<Nnuc; k++)	{
					taxfreq[j][k] = 0;
				} 
			}
		}

		for (int j=0; j<Ntaxa; j++)	{
			for (int k=0; k<Nnuc; k++)	{
				taxfreq[j][k] = 0;
			} 
		}

		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Ntaxa; j++)	{
				int state = GetState(j,i);
				if (state != unknown)	{
					if (pos == -1)	{
						taxfreq[j][GetCodonStateSpace()->GetCodonPosition(0,state)]++;
						taxfreq[j][GetCodonStateSpace()->GetCodonPosition(1,state)]++;
						taxfreq[j][GetCodonStateSpace()->GetCodonPosition(2,state)]++;
					}
					else if (pos > 2)	{
						cerr << "error in CodonSequenceAlignment::NucleotideCompositionHeterogeneity : " << pos << '\n';
						exit(1);
					}
					else	{
						taxfreq[j][GetCodonStateSpace()->GetCodonPosition(pos,state)]++;
					}
				}
			}
		}
				
		// make global freqs out of tax-specific freqs
		double* globalfreq = new double[Nnuc];
		for (int k=0; k<Nnuc; k++)	{
			globalfreq[k] = 0;
			for (int j=0; j<Ntaxa; j++)	{
				globalfreq[k] += taxfreq[j][k];
			}
		} 

		// normalise
		double total = 0;
		for (int k=0; k<Nnuc; k++)	{
			total += globalfreq[k];
		}
		for (int k=0; k<Nnuc; k++)	{
			globalfreq[k] /= total;
		}
		for (int j=0; j<Ntaxa; j++)	{
			double total = 0;
			for (int k=0; k<Nnuc; k++)	{
				total += taxfreq[j][k];
			}
			for (int k=0; k<Nnuc; k++)	{
				taxfreq[j][k] /= total;
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
			for (int k=0; k<Nnuc; k++)	{
				double tmp = (taxfreq[j][k] - globalfreq[k]);
				dist += tmp * tmp;
			}
			if (maxdist < dist)	{
				maxdist = dist;
			}
		}

		delete[] globalfreq;
		if (!comp)	{
			for (int j=0; j<Ntaxa; j++)	{
				delete[] taxfreq[j];
			}
			delete[] taxfreq;
		}

		return maxdist;
	}

	double AminoAcidCompositionalHeterogeneity(ostream* os)	{

		double** taxfreq = new double*[Ntaxa];
		for (int j=0; j<Ntaxa; j++)	{
			taxfreq[j] = new double[Naa];
			for (int k=0; k<Naa; k++)	{
				taxfreq[j][k] = 0;
			} 
		}

		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Ntaxa; j++)	{
				int state = GetState(j,i);
				if (state != unknown)	{
					taxfreq[j][GetCodonStateSpace()->Translation(state)]++;
				}
			}
		}
				
		// make global freqs out of tax-specific freqs
		double* globalfreq = new double[Naa];
		for (int k=0; k<Naa; k++)	{
			globalfreq[k] = 0;
			for (int j=0; j<Ntaxa; j++)	{
				globalfreq[k] += taxfreq[j][k];
			}
		} 

		// normalise
		double total = 0;
		for (int k=0; k<Naa; k++)	{
			total += globalfreq[k];
		}
		for (int k=0; k<Naa; k++)	{
			globalfreq[k] /= total;
		}
		for (int j=0; j<Ntaxa; j++)	{
			double total = 0;
			for (int k=0; k<Naa; k++)	{
				total += taxfreq[j][k];
			}
			for (int k=0; k<Naa; k++)	{
				taxfreq[j][k] /= total;
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
			for (int k=0; k<Naa; k++)	{
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

	CodonStateSpace* GetCodonStateSpace()	{
		// return static_cast<CodonStateSpace*>(statespace);
		return (CodonStateSpace*) (statespace);
	}

	void ToStream(ostream& os);

	private:

	SequenceAlignment* DNAsource;

};

/*
class GCContinuousData : public ContinuousData {

	public:

	GCContinuousData(CodonSequenceAlignment* from, int pos)	{
		taxset = from->GetTaxonSet();
		double** freq = new double*[taxset->GetNtaxa()];
		for (int i=0; i<taxset->GetNtaxa(); i++)	{
			freq[i] = new double[Nnuc];
		}
		from->NucleotideCompositionalHeterogeneity(0,pos,freq);
		Data = new double*[taxset->GetNtaxa()];
		Nsite = 1;
		for (int i=0; i<taxset->GetNtaxa(); i++)	{
			Data[i] = new double[1];
			double tmp = freq[i][1] + freq[i][2];
			Data[i][0] = tmp;
			// Data[i][0] = log(tmp / (1-tmp));
			cerr << taxset->GetTaxon(i) << '\t' << tmp << '\t' << Data[i][0] << '\n';
		}
	}
};
*/

#endif

