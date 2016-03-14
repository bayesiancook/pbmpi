#include "PartitionedGTRSubstitutionProcess.h"

//-------------------------------------------------------------------------
//	* elementary computations on conditional likelihood vectors
//	(CPU level 1)
//-------------------------------------------------------------------------

// set the vector uniformly to 1
void PartitionedGTRSubstitutionProcess::Reset(double*** t, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmp = t[i][j];
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					(*tmp++) = 1.0;
					// tmp[k] = 1.0;
				}
				*tmp = 0;
				tmp -= nstate;
				// tmp[GetNstate(i)] = 0;
			}
		}
	}
}

// initialize the vector according to the data observed at a given leaf of the tree (contained in const int* state)
// state[i] == -1 means 'missing data'. in that case, conditional likelihoods are all 1
void PartitionedGTRSubstitutionProcess::Initialize(double*** t, const int* state, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmp = t[i][j];
				int nstate = GetNstate(i);
				tmp[nstate] = 0;
				if (state[i] == -1)	{
					for (int k=0; k<nstate; k++)	{
						(*tmp++) = 1.0;
						// tmp[k] = 1.0;
					}
					tmp -= nstate;
				}
				else	{
					for (int k=0; k<nstate; k++)	{
						(*tmp++) = 0;
						// tmp[k] = 0;
					}
					tmp -= nstate;
					tmp[state[i]] = 1.0;
				}
			}
		}
	}
}

// multiply two conditional likelihood vectors, term by term
void PartitionedGTRSubstitutionProcess::Multiply(double*** from, double*** to, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmpfrom = from[i][j];
				double* tmpto = to[i][j];
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					(*tmpto++) *= (*tmpfrom++);
					// tmpto[k] *= tmpfrom[k];
				}
				*tmpto += *tmpfrom;
				tmpto -= nstate;
				tmpfrom -= nstate;
				// tmpto[GetNstate(i)] += tmpfrom[GetNstate(i)];
			}
		}
	}
}

// multiply a conditional likelihood vector by the (possibly site-specific) stationary probabilities of the process
void PartitionedGTRSubstitutionProcess::MultiplyByStationaries(double*** to, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		{
			const double* stat = GetStationary(i);
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmpto = to[i][j];
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{
						(*tmpto++) *= (*stat++);
						// tmpto[k] *= stat[k];
					}
					tmpto -= nstate;
					stat -= nstate;
				}
			}
		}
	}
}

// to avoid numerical errors: all entries for a given site and a given rate
// are divided by the largest among them
// and the residual is stored in the last entry of the vector
void PartitionedGTRSubstitutionProcess::Offset(double*** t, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmp = t[i][j];
				double max = 0;
				for (int k=0; k<GetNstate(i); k++)	{
					if (tmp[k] <0)	{
						cerr << "error in pruning: negative prob : " << tmp[k] << "\n";
						exit(1);
						tmp[k] = 0;
					}
					if (max < tmp[k])	{
						max = tmp[k];
					}
				}
				if (max == 0)	{
					max = 1e-12;
					/*
					cerr << "error in pruning: null likelihood\n";
					exit(1);
					*/
				}
				for (int k=0; k<GetNstate(i); k++)	{
					tmp[k] /= max;
				}
				tmp[GetNstate(i)] += log(max);
			}
		}
	}
}

//-------------------------------------------------------------------------
//	* likelihood computation (last step, once the recursion has proceeded throughout the entire tree)
//	(CPU level 2)
//-------------------------------------------------------------------------


double PartitionedGTRSubstitutionProcess::ComputeLikelihood(double*** aux, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		if (condalloc)	{
			int j = ratealloc[i];
			double* t = aux[i][j];
			double tot = 0;
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				tot += (*t++);
				// tot += t[k];
			}
			if (tot == 0)	{
				// dirty !
				tot = 1e-12;
				/*
				cerr << "pruning : 0 \n";
				for (int k=0; k<GetNstate(i); k++)	{
					cerr << t[k] << '\n';
				}
				exit(1);
				*/
			}
			sitelogL[i] = log(tot) + (*t);
			// sitelogL[i] = log(tot) + t[GetNstate(i)];
			t -= nstate;
		}
		else	{
			double max = 0;
			double* logl = condsitelogL[i];
			for (int j=0; j<GetNrate(i); j++)	{
				double* t = aux[i][j];
				double tot = 0;
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					tot += (*t++);
					// tot += t[k];
				}
				if (tot == 0)	{
					// dirty !
					tot = 1e-12;
					/*
					cerr << "pruning : 0 \n";
					for (int k=0; k<GetNstate(i); k++)	{
						cerr << t[k] << '\n';
					}
					exit(1);
					*/
				}
				logl[j] = log(tot) + (*t);
				// logl[j] = log(tot) + t[GetNstate(i)];
				t -= nstate;
				if ((!j) || (max < logl[j]))	{
					max = logl[j];
				}
			}
			double total = 0;
			for (int j=0; j<GetNrate(i); j++)	{
				total += GetRateWeight(i,j) * exp(logl[j] - max);
			}
			sitelogL[i] = log(total) + max;
		}
	}

	logL = 0;
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		if(!sitemask[i-sitemin])
			logL += sitelogL[i];
	}
	return logL;
}

void PartitionedGTRSubstitutionProcess::Propagate(double*** from, double*** to, double time, bool condalloc)	{

	// propchrono.Start();
	int i,j,k,l,offset;
	double length,max,maxup;
	const int nstate = GetMatrix(sitemin)->GetNstate();
	// double* bigaux = new double[(sitemax - sitemin) * GetNrate(0) * nstate];
	double* aux = new double[GetNsite() * GetNrate(0) * nstate];
	for (int i=sitemin; i<sitemax; i++)	{
		if(!sitemask[i-sitemin])
		{
			SubMatrix* matrix = GetMatrix(i);
			double** eigenvect = matrix->GetEigenVect();
			double** inveigenvect = matrix->GetInvEigenVect();
			double* eigenval = matrix->GetEigenVal();
			for(j=0; j<GetNrate(i); j++)	{
				if ((!condalloc) || (ratealloc[i] == j))	{
					double* up = from[i][j];
					double* down = to[i][j];
					//SubMatrix* matrix = GetMatrix(i);
					length = time * GetRate(i,j);

					//double** eigenvect = matrix->GetEigenVect();
					//double** inveigenvect= matrix->GetInvEigenVect();
					//double* eigenval = matrix->GetEigenVal();

					// substitution matrix Q = P L P^{-1} where L is diagonal (eigenvalues) and P is the eigenvector matrix
					// we need to compute
					// down = exp(length * Q) . up
					// which we express as
					// down = P ( exp(length * L) . (P^{-1} . up) )

					// thus we successively do the following matrix.vector products

					// P^{-1} . up  -> aux
					// exp(length * L) . aux  -> aux 	(where exp(length*L) is diagonal, so this is linear)
					// P . aux -> down

					/*
					int nstate = GetNstate();
					double* aux = new double[nstate];
					*/
					//double* aux = bigaux + nstate * (i*GetNrate(0)  + j);
					offset = nstate*(i*GetNrate(0) + j);
					// P^{-1} . up  -> aux
					//double* tmpaux = aux;
					for(k=0; k<nstate; k++)	{
						//(*tmpaux++) = 0;
						aux[offset+k] = 0.0;
					}
					//tmpaux -= nstate;
					//double* tmpup = up;
					for(k=0; k<nstate; k++)	{
						//double* tmpinveigen = inveigenvect[i];
						for(l=0; l<nstate; l++)	{
							//(*tmpaux) += (*tmpinveigen++) * (*tmpup++);
							aux[offset+k] += inveigenvect[k][l] * up[l];
						}
						//tmpaux++;
						//tmpup -= nstate;
					}
					//tmpaux -= nstate;

					// exp(length * L) . aux  -> aux
					//double* tmpval = eigenval;
					for(k=0; k<nstate; k++)	{
						//(*tmpaux++) *= exp(length * (*tmpval++));
						aux[offset+k] *= exp(length * eigenval[k]);
					}
					//tmpaux -= nstate;
					//tmpval -= nstate;

					// P . aux -> down
					//double* tmpdown = down;
					for(k=0; k<nstate; k++)	{
						//(*tmpdown++) = 0;
						down[k] = 0.0;
					}
					//tmpdown -= nstate;

					for(k=0; k<nstate; k++)	{
						//double* tmpeigen = eigenvect[i];
						for(l=0; l<nstate; l++)	{
							//(*tmpdown) += (*tmpeigen++) * (*tmpaux++);
							down[k] += eigenvect[k][l] * aux[offset+l];
						}
						//tmpdown++;
						//tmpaux -= nstate;
					}
					//tmpdown -= nstate;

					// exit in case of numerical errors
					for(k=0; k<nstate; k++)	{
						if (isnan(down[k]))	{
							cerr << "error in back prop\n";
							for(l=0; l<nstate; l++)	{
								cerr << up[l] << '\t' << down[l] << '\t' << matrix->Stationary(l) << '\n';
							}
							exit(1);
						}
					}
					maxup = 0.0;
					for(k=0; k<nstate; k++)	{
						if (up[k] < 0.0)	{
							cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
							exit(1);
						}
						if (maxup < up[k])	{
							maxup = up[k];
						}
					}
					max = 0.0;
					for(k=0; k<nstate; k++)	{
						if (down[k] < 0.0)	{
							infprobcount++;
							down[k] = 0.0;
						}
						if (max < down[k])	{
							max = down[k];
						}
					}
					/*if (maxup == 0.0)	{
						cerr << "error in backward propagate: null up array\n";
						for(l=0; l<nstate; l++)	{
							cerr << matrix->Stationary(l) << '\n';
						}
						cerr << time << '\t' << length << '\n';
						exit(1);
					}
					if (max == 0.0)	{
						cerr << "error in backward propagate: null array\n";
						for(k=0; k<nstate; k++)	{
							cerr << up[k] << '\t' << down[k] << '\n';
						}
						cerr << length << '\n';
						cerr << '\n';
						exit(1);
					}*/


					// this is the offset (in log)
					down[nstate] = up[nstate];
				}
			}
		}
	}

	delete[] aux;
	// propchrono.Stop();
}
