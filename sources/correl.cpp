
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <stdio.h>
#include "correlation.h"
using namespace std;

void compareChainsConvergence(string outname, Correlation ** corr, int nchain, double& disc, double& overlap, double& effsize)
{

	int nbp = 0;
	for(int chain=0; chain<nchain; chain++)	{
		int n = corr[chain]->getNbParameter();
		if (!chain)	{
			nbp = n;
		}
		else	{
			if (nbp != n)	{
			      cerr << "ERROR in compChain compareChainsConvergence(), the two chain have not the same dimensionality\n";
			      cerr.flush();
			      exit(0);
			}
		}
	}

  // int nbyes=0;
  // double tmp;
  double meaneffsize[nbp];
  double absDiff[nbp];
  double meanse[nbp];
  double discrepancy[nbp];

 double maxCIoverlap[nbp];

  int nbrecci[20];

  for(int j=0;j<20;j++)
    nbrecci[j]=0;
  double mean[nchain], se[nchain];

  ofstream os_conv((outname + ".contdiff").c_str());
  os_conv << "name                effsize\trel_diff\n";
  os_conv << '\n';

disc = 0;
overlap = 100;
effsize = corr[0]->getNbSample();

  for(int i=0;i<nbp;i++)
    {
	
	for(int chain=0; chain<nchain; chain++)	{
	      mean[chain]=corr[chain]->getMean(i);
	      se[chain] = corr[chain]->getVariance(i);
	}
      absDiff[i]=0;
	for(int chain1=0; chain1<nchain; chain1++)	{
		for(int chain2=chain1+1; chain2<nchain; chain2++)	{
			double tmp = fabs(mean[chain2] - mean[chain1]);
			if (absDiff[i] < tmp)	{
				absDiff[i] = tmp;
			}
		}
	}

	double max = 0;
	for(int chain=0; chain<nchain; chain++)	{
		if (max < mean[chain])	{
			max = mean[chain];
		}
	}

	meanse[i] = 0;
	for(int chain=0; chain<nchain; chain++)	{
		if (se[chain] > 1e-12)	{
			meanse[i] += sqrt(se[chain]);
		}
	}
	meanse[i] /= nchain;
	if (meanse[i] > 1e-12)	{
		discrepancy[i] = absDiff[i] / meanse[i];
	}
	else	{
		discrepancy[i] = 0;
	}

	meaneffsize[i] = 0;
	for(int chain=0; chain<nchain; chain++)	{
		meaneffsize[i] += corr[chain]->getEffectiveSize(i);
	}
	meaneffsize[i] /= nchain;

	  double ciinf[nchain], cisup[nchain];
	for(int chain=0; chain<nchain; chain++)	{
	      ciinf[chain]=corr[chain]->getInfCI(i);
	      cisup[chain]=corr[chain]->getSupCI(i);
	}

	maxCIoverlap[i] = 0;
	for(int chain1=0; chain1<nchain; chain1++)	{
		for(int chain2=chain1+1; chain2<nchain; chain2++)	{
			
			double CIoverlap, inf, sup, medi, meds;
			
			double ciinf1 = ciinf[chain1];
			double ciinf2 = ciinf[chain2];
			double cisup1 = cisup[chain1];
			double cisup2 = cisup[chain2];

		      if(ciinf1<ciinf2)
			{
			  inf=ciinf1;
			  medi=ciinf2;
			} else {
			  inf=ciinf2;
			  medi=ciinf1;
			}
		      if(cisup1<cisup2)
			{
			  sup=cisup2;
			  meds=cisup1;
			} else {
			  sup=cisup1;
			  meds=cisup2;
			}
	      if(medi>meds)
		CIoverlap=0; // aucun recouvrement      
	      else
		{
		  if(sup-inf != 0)
		    CIoverlap=(meds-medi)*100/(sup-inf);
		  else 
		    CIoverlap=100;
		}
	       for(int j=0;j<20;j++)
		{
		  if(CIoverlap<=5*(j+1))
		    {
		      nbrecci[j]++;
		      break;
		    }
		}

		if (maxCIoverlap[i] < CIoverlap)	{
			maxCIoverlap[i] = CIoverlap;
		}
	}
	}
	      os_conv << corr[0]->getParameterName(i);
	      for (unsigned int k=corr[0]->getParameterName(i).length(); k<20; k++) os_conv << ' ';
	      os_conv << (int) meaneffsize[i] << '\t' << '\t'  << discrepancy[i] << "\n";

	if (disc < discrepancy[i])	{
		disc = discrepancy[i];
	}
	if (effsize > meaneffsize[i])	{
		effsize = meaneffsize[i];
	}
	if (overlap > maxCIoverlap[i])	{
		overlap = maxCIoverlap[i];
	}
    } 
  os_conv.close();
  string cat = "cat " + outname + ".contdiff";
  system(cat.c_str());
}


int SamCompare(int nchain, int burnin, int stop, string* ChainName, double& disc, double& overlap, double& effsize, string outname)	{

	// cerr << "burnin : " << burnin << '\n';
	// cerr << "stop : " << stop << '\n';
  Correlation** corr = new Correlation*[nchain];
  for (int chain=0; chain<nchain; chain++)	{
	corr[chain] = new Correlation();
	corr[chain]->getParameters(ChainName[chain],burnin,stop);
	corr[chain]->computeCovariance();
	corr[chain]->computeEffectiveSize();
	corr[chain]->computeCI(95);
	// corr[chain]->outputResult();
  }

      compareChainsConvergence(outname,corr,nchain,disc,overlap,effsize);

  for (int chain=0; chain<nchain; chain++)	{
	delete corr[chain];
  }
  delete[] corr;

  return 1;
}
