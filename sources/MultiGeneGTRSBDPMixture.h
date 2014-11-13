
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENEGTR_H
#define MULTIGENEGTR_H

#include "ExpoConjugateGTRGeneGammaPhyloProcess.h"
#include "ExpoConjugateGTRSBDPProfileProcess.h"
#include "MultiGeneMixture.h"

class MultiGeneGTRSBDPMixture : public virtual ExpoConjugateGTRSBDPProfileProcess, public virtual MultiGeneMixture	{

	public:

	MultiGeneGTRSBDPMixture(string indatafile, string treefile, string inname, int innratecat, int fixtopo, int kappaprior, double gibbsfactor, int dc, int me, int np);
	MultiGeneGTRSBDPMixture(istream& is, int me, int np);

	void Trace(ostream& os) {

		UpdateOccupancyNumbers();

		os << ((int) (chronototal.GetTime() / 1000));
		os << '\t' << ((int) (propchrono.GetTime() / 1000));
		if (chronototal.GetTime())	{
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
		}
		else	{
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanLength();
		os << '\t' << GetMeanAlpha();
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}
		os << '\n';
		os.flush();
	}

	ExpoConjugateGTRGeneGammaPhyloProcess* GetGTRProcess(int gene)	{
		ExpoConjugateGTRGeneGammaPhyloProcess* tmp = dynamic_cast<ExpoConjugateGTRGeneGammaPhyloProcess*> (process[gene]);
		if (! tmp)	{
			cerr << "error : null gene process\n";
			exit(1);
		}
		return tmp;
	}

	void TraceHeader(ostream& os) {
		os << "#time\ttimeinmix\tprop\tlnL\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		// os << "\tkappa\tallocent";
		os << '\n'; 
	}

	void ToStream(ostream& os)	{
		ExpoConjugateGTRSBDPProfileProcess::ToStream(os);
		GlobalToStream(os);
	}

	void FromStream(istream& is)	{
		ExpoConjugateGTRSBDPProfileProcess::FromStream(is);
		GlobalFromStream(is);
	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		GlobalGeneMove();

		propchrono.Start();
		ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		propchrono.Stop();

		LengthRelRateMove(1,10);
		LengthRelRateMove(0.1,10);
		LengthRelRateMove(0.01,10);

		GlobalUpdateParameters();
		GlobalUnfold();

		chronototal.Stop();

		// Trace(cerr);
	}

	double* GetEmpiricalFreq()	{
		cerr << "error : in get empirical freq\n";
		exit(1);
		return 0;
		// return empfreq;
	}

	protected:

        virtual void SlaveExecute(MESSAGE);

	virtual void SlaveMixMove()	{
		ExpoConjugateGTRSBDPProfileProcess::SlaveMixMove();
	}

	double LengthRelRateMove(double tuning, int nrep);
	void SlaveLengthFactorMove();

	void Create(int inNsite, int Nstate, int nratecat, int dc);
	void CreateSuffStat();
	void DeleteSuffStat();

	void GlobalUpdateParameters();
	void SlaveUpdateParameters();

	void GlobalUpdateSiteProfileSuffStat();
	void SlaveUpdateSiteProfileSuffStat();

	void GlobalUpdateRRSuffStat();
	void SlaveUpdateRRSuffStat();
	void UpdateRRSuffStat() {
		cerr << "error: in update site rr suff stat of global mixture\n";
		exit(1);
	}

	int* GetSiteProfileSuffStatCount(int site)	{
		return siteprofilesuffstatcount[site];
	}

	double* GetSiteProfileSuffStatBeta(int site)	{
		return siteprofilesuffstatbeta[site];
	}

	int* allocsiteprofilesuffstatcount;
	double* allocsiteprofilesuffstatbeta;
	int** siteprofilesuffstatcount;
	double** siteprofilesuffstatbeta;

	int* tmpallocsiteprofilesuffstatcount;
	double* tmpallocsiteprofilesuffstatbeta;
	int** tmpsiteprofilesuffstatcount;
	double** tmpsiteprofilesuffstatbeta;
	int* tmprrsuffstatcount;
	double* tmprrsuffstatbeta;

};

#endif
