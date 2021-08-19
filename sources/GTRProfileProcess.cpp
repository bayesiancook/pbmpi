
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
#include "GTRProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GTRProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void GTRProfileProcess::Create(int innsite, int indim)	{
	if (! rr)	{
		ProfileProcess::Create(innsite,indim);
		Nrr = GetDim() * (GetDim()-1) / 2;
		rr = new double[Nrr];
        emprralpha = new double[Nrr];
        emprrbeta = new double[Nrr];
		SampleRR();
	}
}

void GTRProfileProcess::Delete()	{
	if (rr)	{
		DeleteMatrices();
        delete[] emprrbeta;
        delete[] emprralpha;
		delete[] rr;
		rr = 0;
		ProfileProcess::Delete();
	}
}

double GTRProfileProcess::LogRRPrior()	{
    if (fixrr)  {
        cerr << "error: in GTRProfileProcess::LogRRPrior\n";
        exit(1);
    }
	double total = 0;
    if (profilefrac == 1.0) {
        for (int i=0; i<GetNrr(); i++)	{
            total -= rr[i];
        }
    }
    else    {
        for (int i=0; i<GetNrr(); i++)	{
            double a = profilefrac + (1.0-profilefrac) * emprralpha[i];
            double b = profilefrac + (1.0-profilefrac) * emprrbeta[i];
            total += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(rr[i]) - b*rr[i];
        }
	}
	return total;
}

void GTRProfileProcess::SampleRR()	{
	for (int i=0; i<GetNrr(); i++)	{
		// rr[i] = LG_RR[i];
		rr[i] = rnd::GetRandom().sExpo();
	}
}

void GTRProfileProcess::PriorSampleRR()	{
    if (profilefrac == 1.0) {
        for (int i=0; i<GetNrr(); i++)	{
            // rr[i] = LG_RR[i];
            rr[i] = rnd::GetRandom().sExpo();
        }
    }
    else    {
        for (int i=0; i<GetNrr(); i++)	{
            double a = profilefrac + (1.0-profilefrac) * emprralpha[i];
            double b = profilefrac + (1.0-profilefrac) * emprrbeta[i];
            rr[i] = rnd::GetRandom().Gamma(a,b);
        }
    }
}

void GTRProfileProcess::SetRR(string type)	{

	rrtype = type;
	if (type != "None")	{
		fixrr = true;
		if ((type == "WAG") || (type == "wag"))	{
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
				rr[i]= WAG_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg6") || (type == "CG6"))	{
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
				rr[i]= CG6RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg10") || (type == "CG10"))	{
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
				rr[i]= CG10RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg20") || (type == "CG20"))	{
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
				rr[i]= CG20RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg30") || (type == "CG30"))	{
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
				rr[i]= CG30RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg40") || (type == "CG40"))	{
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
				rr[i]= CG40RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg50") || (type == "CG50"))	{
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
				rr[i]= CG50RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg60") || (type == "CG60"))	{
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
				rr[i]= CG60RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "LG") || (type == "lg"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : lg only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= LG_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "JTT") || (type == "jtt"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : jtt only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= JTT_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "mtREV") || (type == "mtrev"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : mtrev only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtREV_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "mtZOA") || (type == "mtzoa"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : mtzoa only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtZOA_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "mtART") || (type == "mtart"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : mtart only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtART_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
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
					rr[rrindex(permut[k],permut[l],GetDim())] = tmp;
					
				}
			}
		}	
	}
}


