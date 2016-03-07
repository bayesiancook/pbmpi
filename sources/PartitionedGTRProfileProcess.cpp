
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Parallel.h"
#include "PartitionedGTRProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionedGTRProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PartitionedGTRProfileProcess::Create(int indim, PartitionScheme inscheme)	{
	if (! rr)	{
		ProfileProcess::Create(inscheme.GetNsite(),indim);
		PartitionProcess::Create(inscheme);
		Nrr = GetDim() * (GetDim()-1) / 2;

		rr = new double* [GetNpart()];
		fixrr = new bool[GetNpart()];

		for(int i = 0; i < GetNpart(); i++)
		{
			fixrr[i] = false;
			rr[i] = new double[Nrr];
		}
		SampleRR();

		for(int p = 0; p < inscheme.Npart; p++)
			SetRR(p, inscheme.partType[p]);
	}
}

void PartitionedGTRProfileProcess::Delete()	{
	if (rr)	{
		DeleteMatrices();
		for(int i = 0; i < GetNpart(); i++)
		{
			delete[] rr[i];
		}
		delete[] rr;
		delete[] fixrr;

		rr = 0;
		fixrr = 0;

		PartitionProcess::Delete();
		ProfileProcess::Delete();
	}
}

double PartitionedGTRProfileProcess::LogRRPrior()	{
	double total = 0;
	for(int p = 0; p < GetNpart(); p++)
	{
		if(!fixrr[p])
			for (int i=0; i<GetNrr(); i++)	{
				total -= rr[p][i];
			}
	}
	return total;
}

void PartitionedGTRProfileProcess::SampleRR()	{
	for(int p = 0; p < GetNpart(); p++)
	{
		if(!fixrr[p])
			for (int i=0; i<GetNrr(); i++)	{
				// rr[i] = LG_RR[i];
				rr[p][i] = rnd::GetRandom().sExpo();
			}
	}
}

void PartitionedGTRProfileProcess::SetRR(int inpart, string type)	{

	int p = inpart;
	if (type != "None")	{
		fixrr[inpart] = true;

		const double * RR = 0;

		if (type == "dayhoff")	{
			RR = DAYHOFF_RR;
		}

		else if (type == "dcmut")	{
			RR = DCMUT_RR;
		}

		else if (type == "jtt")	{
			RR = JTT_RR;
		}

		else if (type == "mtrev")	{
			RR = MTREV_RR;
		}

		else if (type == "wag")	{
			RR = WAG_RR;
		}

		else if (type == "rtrev")	{
			RR = RTREV_RR;
		}

		else if (type == "cprev")	{
			RR = CPREV_RR;
		}

		else if (type == "vt")	{
			RR = VT_RR;
		}

		else if (type == "blosum62")	{
			RR = BLOSUM62_RR;
		}

		else if (type == "mtmam")	{
			RR = MTMAM_RR;
		}

		else if (type == "lg")	{
			RR = LG_RR;
		}

		else if (type == "mtart")	{
			RR = MTART_RR;
		}

		else if (type == "mtzoa")	{
			RR = MTZOA_RR;
		}

		else if (type == "pmb")	{
			RR = PMB_RR;
		}

		else if (type == "hivb")	{
			RR = HIVB_RR;
		}

		else if (type == "hivw")	{
			RR = HIVW_RR;
		}

		else if (type == "jttdcmut")	{
			RR = JTTDCMUT_RR;
		}

		else if (type == "flu")	{
			RR = FLU_RR;
		}

		else if (type == "cg6")	{
			RR = CG6RR;
		}

		else if (type == "cg10")	{
			RR = CG10RR;
		}

		else if (type == "cg20")	{
			RR = CG20RR;
		}

		else if (type == "cg30")	{
			RR = CG30RR;
		}

		else if (type == "cg40")	{
			RR = CG40RR;
		}

		else if (type == "cg50")	{
			RR = CG50RR;
		}

		else if (type == "cg60")	{
			RR = CG60RR;
		}



		if(RR != 0){
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[p][i]= RR[i];
				total += rr[p][i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[p][i] /= total / Nrr;
			}
		}
		else	{
			ifstream is(type.c_str());
			if (!is)	{
				cerr << "unrecognized file for relative rates : " << type << '\n';
				exit(1);
			}

			int permut[GetDim()];
			for (int k=0; k<GetDim(); k++)	{
				string c;
				is >> c;
				permut[k] = GetStateSpace()->GetState(c);
			}
			for (int k=0; k<GetDim()-1; k++)	{
				for (int l=k+1; l<GetDim(); l++)	{
					double tmp;
					is >> tmp;
					if (tmp < 0)	{
						if (! GetMyid())	{
							cerr << "error when reading exchangeabilities from " << type << '\n';
							cerr << tmp << '\n';
							cerr << '\n';
						}
						MPI_Finalize();
						exit(1);
					}
					rr[p][rrindex(permut[k],permut[l],GetDim())] = tmp;
					
				}
			}
		}	
	}
	else
	{
		fixrr[inpart] = false;
	}

	nfreerr = GetNpart();
	for(int p = 0; p < GetNpart(); p++)
	{
		nfreerr -= fixrr[p];
	}
}


