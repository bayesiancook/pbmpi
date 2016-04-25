
#include "phylo.h"

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	double cutoff = atof(argv[2]);

	string alphabet;
	is >> alphabet;
	int Nstate = alphabet.size();
	cerr << "nstate : " << Nstate << '\n';

	double stat[Nstate];
	double mat[Nstate][Nstate];
	for (int i=0; i<Nstate; i++)	{
		mat[i][i] = 0;
	}

	for (int i=0; i<Nstate; i++)	{
		is >> stat[i];
	}

	for (int i=0; i<Nstate; i++)	{
		for (int j=i+1; j<Nstate; j++)	{
			double tmp;
			is >> tmp;
			mat[i][j] = tmp;
			// mat[i][j] = tmp * stat[i]*stat[j];
			mat[j][i] = mat[i][j];
		}
	}
	
	int Nrr = Nstate * (Nstate - 1) / 2;
	
	int pair[Nrr][2];
	double pairscore[Nrr];

	int k = 0;
	int keep = 0;
	while (k<Nrr)	{

		double max = 0;
		int imax = 0;
		int jmax = 0;
		for (int i=0; i<Nstate; i++)	{
			for (int j=i+1; j<Nstate; j++)	{
				if (max < mat[i][j])	{
					max = mat[i][j];
					imax = i;
					jmax = j;
				}
			}
		}
		pair[k][0] = imax;
		pair[k][1] = jmax;
		pairscore[k] = max;
		mat[imax][jmax] = 0;
		if (max > cutoff)	{
			keep++;
		}
		k++;
	}

	cerr << keep << '\n';

	int cluster[keep+Nstate][Nstate];
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			cluster[i][j] = 0;
		}
		cluster[i][i] = 1;
		cout << alphabet[i] << '\n';
	}
		
	for (int i=0; i<keep; i++)	{
		for (int j=0; j<Nstate; j++)	{
			cluster[i+Nstate][j] = 0;
		}
		cluster[i+Nstate][pair[i][0]] = 1;
		cluster[i+Nstate][pair[i][1]] = 1;
		cout << alphabet[pair[i][0]] << alphabet[pair[i][1]] << '\n';
	}

	int K = keep + Nstate;
	double dist[K][K];
	for (int i=0; i<K; i++)	{
		dist[i][i] = 0;
	}
	for (int i=0; i<K; i++)	{
		for (int j=i+1; j<K; j++)	{
			double max = 0;
			for (int k=0; k<Nstate; k++)	{
				if (cluster[i][k])	{
					for (int l=0; l<Nstate; l++)	{
						if ((!cluster[j][k]) && cluster[j][l])	{
							if (max < mat[k][l])	{
								max = mat[k][l];
							}
						}
					}
				}
			}
			dist[i][j] = max;
			dist[j][i] = max;
		}
	}

	int left[K];
	for (int i=0; i<K; i++)	{
		left[i] = 1;
	}

	double max = 1;
	int L = K;
	while ((L > 1) && (max > 0))	{
		max = 0;
		int imax = 0;
		int jmax = 0;
		for (int i=0; i<K; i++)	{
			if (left[i])	{
				for (int j=0; j<K; j++)	{
					if (left[j])	{
						if (max < dist[i][j])	{
							max = dist[i][j];
							imax = i;
							jmax = j;
						}
					}
				}
			}
		}

		int tmp[Nstate];
		for (int i=0; i<Nstate; i++)	{
			tmp[i] = (cluster[imax][i] || cluster[jmax][i]);
		}
		int same = -1;
		for (int j=0; j<K; j++)	{
			if (left[j])	{
				int eq = 1;
				for (int i=0; i<Nstate; i++)	{
					if (tmp[i] != cluster[j][i])	{
						eq = 0;
					}
				}
				if (eq)	{
					same = j;
				}
			}
		}
		if (same != -1)	{
			if (same != imax)	{
				left[imax] = 0;
				K--;
			}
			if (same != jmax)	{
				left[jmax] = 0;
				K--;
			}
		}
		else	{
			for (int i=0; i<Nstate; i++)	{
				cluster[imax][i] = tmp[i];
				if (cluster[imax][i])	{
					cout << alphabet[i];
				}
			}
			cout << '\n';
			for (int j=0; j<K; j++)	{
				if ((j!= imax) && (j != jmax) && left[j])	{
					if (dist[imax][j] < dist[jmax][j])	{
						dist[imax][j] = dist[jmax][j];
						dist[j][imax] = dist[imax][j];
					}
				}
			}
			left[jmax] = 0;
			L--;
		}

		/*
		for (int i=0; i<K; i++)	{
			for (int j=i+1; j<K; j++)	{
				if (left[i] && left[j])	{
					int eq = 1;
					for (int k=0; k<Nstate; k++)	{
						if (cluster[i][k] != cluster[j][k])	{
							eq = 0;
						}
					}
					if (eq)	{
						left[j] = 0;
						L--;
					}
				}
			}
		}
		*/
	}
}

