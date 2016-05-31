
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "correlation.h"

Correlation::Correlation(double ci)
{
  burnin=0;
  nbsample=0;
  nbparameter=0;
  if(ci==-1)
    {
      CI=defaultCI;
    } else {
      CI=ci;
    }
  isSort=No;
  parameterName=NULL;
  parameters=NULL;
  sortparameters=NULL;
  covparam=NULL;
  covnorm=NULL;
  meanparam=NULL;
  isConstant=NULL;
  effectiveSize=NULL;
  variance=NULL;
  CIbuffer=NULL;
}

void Correlation::init()
{

  parameterName=new string[nbparameter];
  parameters=new double*[nbparameter];
  sortparameters=new double*[nbparameter];
  covparam=new double*[nbparameter];
  covnorm=new double*[nbparameter];
  meanparam=new double[nbparameter];
  isConstant= new Switch[nbparameter];
  effectiveSize=new double[nbparameter];
  variance=new double[nbparameter];
  CIbuffer=new double*[nbparameter];
  for(int i=0;i<nbparameter;i++)
    {
      CIbuffer[i]=new double[2];
      parameters[i]=NULL;
      sortparameters[i]=NULL;
      covparam[i]=NULL;
      covnorm[i]=NULL;
    }
  weight=NULL;
}

Correlation::~Correlation()
{
    for(int i=0;i<nbparameter;i++)
      {
	if(parameters[i]!=NULL)
	  delete[] parameters[i];
	if(sortparameters[i]!=NULL)
	  delete[] sortparameters[i];
	if(CIbuffer[i]!=NULL)
	  delete[] CIbuffer[i];
	if(covparam[i]!=NULL)
	  delete[] covparam[i];
	if(covnorm[i]!=NULL)
	  delete[] covnorm[i];
      }
    delete[] CIbuffer;
    delete[] covparam;
    delete[] parameters;
    delete[] sortparameters;
    delete[] meanparam;
    if(weight!=NULL)
      delete[] weight;
    delete[] effectiveSize;
    delete[] isConstant;
    delete[] variance;
    delete[] covnorm;
    delete[] parameterName;
}

void Correlation::reccurquicksort(double *to, double *buf, int deb, int fin)
{
  if(fin-deb>1)
    {
      int last=fin;
      int first=deb;
      double pivot=to[deb];
      for(int i=deb+1;i<=fin;i++)
	{
	  if(to[i]>pivot)
	    {
	      buf[last]=to[i];
	      last--;
	    } else {
	      buf[first]=to[i];
	      first++;
	    }
	}
      if(first==last)
	buf[first]=pivot;
      else {
	cerr << "Quiksort: erreur, ou mettre le pivot ???\n";
	cerr.flush();
	exit(0);
      }
      for(int i=deb;i<=fin;i++)
	{
	  to[i]=buf[i];
	}
      if(deb!=first)
	reccurquicksort(to,buf, deb,first-1);
      if(fin!=first)
	reccurquicksort(to,buf, first+1,fin);
    } else {
      if(fin!=deb)
	{
	  double pivot;
	  if(to[fin] < to[deb])
	    {
	      pivot=to[fin];
	      to[fin]=to[deb];
	      to[deb]=pivot;
	    }
	}
    }
}

void Correlation::quicksort(double *from, double *to, int s)
{
  for(int i=0;i<s;i++)
    to[i]=from[i];
  double *tmp=new double[s];
  reccurquicksort(to,tmp,0,s-1);
  delete[] tmp;
}

void Correlation::createWeight()
{
  weight=new double[nbsample/2];
}

void Correlation::createCovarianceBuffer()
{
   if(nbsample!=0)
    { 
      for(int i=0;i<nbparameter;i++)
	{
	  covnorm[i]=new double[nbsample/2];
	  for(int j=0;j<nbsample/2;j++)
	    covnorm[i][j]=0;
	  covparam[i]=new double[nbsample/2];
	  for(int j=0;j<nbsample/2;j++)
	    covparam[i][j]=0;
	}
    } else {
      cerr << "ERROR: in Correlation::createCovarianceBuffer, chain size is 0, exit\n";
      cerr.flush();
      exit(0);
    }
}

void Correlation::createParameterBuffer()
{
  if(nbsample!=0)
    {
      for(int i=0;i<nbparameter;i++)
	{
	  parameters[i]=new double[nbsample];
	  sortparameters[i]=new double[nbsample];	  
	}
    } else {
      cerr << "ERROR: in Correlation::createParameterBuffer, chain size is 0, exit\n";
      cerr.flush();
      exit(0);
    }
}


void Correlation::getParameters(string filename, int start, int stop)
{
  chainName=filename;
  double tmp;
  string strtmp;
  try {

    if (start == -1)	{
	burnin = stop / 5;
    }
    else	{
	    burnin=start;
    }
    nbsample = stop-burnin;
    if(nbsample<=0)
      {
	cerr << "ERROR: in Correlation::getParameters, asking sampling from point " << start << " to point "<< stop <<", exiting\n";
	cerr.flush();
	throw(0);
      }
    cerr << filename << '\t' << " burnin : " << burnin << '\t' << "sample size : " << nbsample << '\n';
    ifstream *is=new ifstream;
    is->open(filename.c_str(), ifstream::in);

    char line[100000];
    string strline;
    is->getline(line,100000);
    strline=line;
    istringstream iss1(strline);

    nbparameter=0;

	/*ALL*/
    iss1 >> strtmp;
    iss1 >> strtmp;
    iss1 >> strtmp;
	/*ALL*/

    strtmp="null";

    while(iss1.good())
      {
	iss1 >> strtmp;
	nbparameter++;
      }
    init();
    istringstream iss2(strline);

	/*ALL*/
    iss2 >> strtmp;
    iss2 >> strtmp;
    iss2 >> strtmp;
	/*ALL*/

    for(int j=0;j<nbparameter;j++)
      iss2 >> parameterName[j];

    createParameterBuffer();
    createWeight();

    for(int i=0;i<stop;i++)
      {

	/*ALL*/
	*is >> tmp;
	*is >> tmp;
	*is >> tmp;
	/*ALL*/

	for(int j=0;j<nbparameter;j++)
	  {
	    *is >> tmp;

	    if(i>=burnin)
	      {
		parameters[j][i-burnin]=tmp;
	      }
	  }
      }
    delete is;
  } catch(...) {
    cerr << "ERROR while reading " << filename << "\n";
    cerr.flush();
    exit(0);
  }
}

void Correlation::sortParameters()
{
  for(int i=0;i < nbparameter; i++)
    {
      quicksort(parameters[i],sortparameters[i],nbsample);
    }
  isSort=Yes;
}

void Correlation::getCI(double *distrib,int s, double ci, double *inf, double* sup)
{
  int iinf, isup;
  if(ci==100)
    {
      iinf=0;
      isup=s-1;
    } else {
      iinf=(int)(0.5*(1-0.01*ci)*s);
      isup=(int)(s-(0.5*(1-0.01*ci)*s));
      isup--;
    }
  *inf=distrib[iinf];
  *sup=distrib[isup];
}

void Correlation::computeCI(double ci)
{
  if(isSort==No)
    sortParameters();
  if(ci==-1)
    {
      ci=defaultCI;
    }
  for(int i=0;i<nbparameter;i++)
    {
      getCI(sortparameters[i],nbsample,ci,&CIbuffer[i][0],&CIbuffer[i][1]);
    }
}

void Correlation::computeMean()
{

  for(int j=0;j<nbparameter;j++)
    meanparam[j]=0;
  for(int j=0;j<nbparameter;j++)
    {
      isConstant[j]=Yes;
      for(int i=0;i<nbsample;i++)
	{
	  if(parameters[j][i]!=parameters[j][0])
	    isConstant[j]=No;
	  meanparam[j]+=parameters[j][i];
	}
    }
  for(int j=0;j<nbparameter;j++)
    {
      if( isConstant[j]==No)
	meanparam[j]/=(double)nbsample;
      else
	meanparam[j]=parameters[j][0];
    }
}

void Correlation::computeCovariance()
{
  if(nbsample==0)
    {
      cerr << "ERROR: in Correlation::computeCovariance, chain size is 0, exit\n";
      cerr.flush();
      exit(0);
    } else {
      createCovarianceBuffer();
      computeMean();

      for(int t=0;t<nbsample/2;t++)
	{
	  for(int j=0;j<nbparameter;j++)
	    {
	      for(int i=0;i<nbsample-t;i++)
		{
		  covparam[j][t]+=(parameters[j][i]-meanparam[j])*(parameters[j][i+t]-meanparam[j]);
		}
	      if(isConstant[j]==Yes)
		covparam[j][t]=0;
	    }
	  if(t==0)
	    {
	      for(int j=0;j<nbparameter;j++)
		{
		  variance[j]=covparam[j][0]/(double)(nbsample-1);
		}
	    } 
	  for(int j=0;j<nbparameter;j++)
	    {
	      covparam[j][t]/=(double)(nbsample);
	      // normalisation
	      if(covparam[j][0]!=0)
		covnorm[j][t]=covparam[j][t]/covparam[j][0];
	      else
		covnorm[j][t]=0;
	    }
	} // fin for t

    }
}



void Correlation::computeWeight()
{
    if(nbsample==0)
    {
      cerr << "ERROR: in Correlation::computeWeight, chain size is 0, exit\n";
      cerr.flush();
      exit(0);
    } else {
      for(int i=0;i<nbsample/2;i++)
	{
	  weight[i]=0.5*(1+cos((2*Pi*i)/(double)(nbsample)));
	}
    }
    
}

void Correlation::computeEffectiveSize()
{
  if(nbsample==0)
    {
      cerr << "ERROR: in Correlation::computeEffectiveSize, chain size is 0, exit\n";
      cerr.flush();
      exit(0);
    } else {
      computeWeight();
      double sum;
      for(int i=0;i<nbparameter;i++)
	{
	  sum=1;//covparam[i][0];*covparam[i][0];
	  for(int j=1;j<nbsample/2;j++)
	    {
	      sum+=2*weight[j]*covnorm[i][j];
	    }
	  effectiveSize[i]=nbsample/sum;
	  if(sum<1)
	    effectiveSize[i]=nbsample;
	} 
    }
}

void Correlation::outputResult()
{
  string output = chainName;
  ofstream os(output.c_str());
  os << "\nparam\t\t\tmean\tSD\tSE \teff size\t " << CI <<"%CI\n";
  os << std::setprecision(4) << std::setiosflags(std::ios::fixed | std::ios::showpoint) << /*0.123456789 <<*/ "\n";
  string space;
  for(int i=0;i<nbparameter;i++)
    {
	space = "  ";
      os << i << " " << parameterName[i] << space << "\t" << meanparam[i] << "\t" << sqrt(variance[i]) << "\t" << sqrt(variance[i]/effectiveSize[i]) << "\t" << effectiveSize[i] << "\t" << CIbuffer[i][0]<< "\t" << CIbuffer[i][1]<<"\n";
    }
  os.close();
}

int Correlation::getNbParameter()
{
  return nbparameter;
}

string Correlation::getParameterName(int n)
{
  if(n<nbparameter)
    return parameterName[n];
  else
    return "NULL";
}

double Correlation::getVariance(int n)
{
  return  variance[n];
}

double Correlation::getEffectiveSize(int n)
{
  return effectiveSize[n];
}

double Correlation::getMean(int n)
{
  return meanparam[n];
}

double Correlation::getInfCI(int n)
{
  return CIbuffer[n][0];
}

double Correlation::getSupCI(int n)
{
  return CIbuffer[n][1];
}

double  Correlation::getNbSample()
{
  return nbsample;
}
