
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENEPOISSON_H
#define MULTIGENEPOISSON_H

#include "PoissonGeneGammaPhyloProcess.h"
#include "PoissonSBDPProfileProcess.h"
#include "MultiGeneMixture.h"

class MultiGenePoissonSBDPMixture : public virtual PoissonSBDPProfileProcess, public virtual MultiGeneMixture	{

	public:

	MultiGenePoissonSBDPMixture(string indatafile, string treefile, string inname, int innratecat, int fixtopo, int kappaprior, double gibbsfactor, int dc, int me, int np);
	MultiGenePoissonSBDPMixture(istream& is, int me, int np);

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
		os << '\n';
		os.flush();
	}

	PoissonGeneGammaPhyloProcess* GetPoissonProcess(int gene)	{
		PoissonGeneGammaPhyloProcess* tmp = dynamic_cast<PoissonGeneGammaPhyloProcess*> (process[gene]);
		if (! tmp)	{
			cerr << "error : null gene process\n";
			exit(1);
		}
		return tmp;
	}

	void TraceHeader(ostream& os) {
		os << "#time\ttimeinmix\tprop\tlnL\tlength\talpha\tNmode\tstatent\tstatalpha";
		os << '\n'; 
	}

	void ToStream(ostream& os)	{
		PoissonSBDPProfileProcess::ToStream(os);
		GlobalToStream(os);
	}

	void FromStream(istream& is)	{
		PoissonSBDPProfileProcess::FromStream(is);
		GlobalFromStream(is);
	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		GlobalGeneMove();

		propchrono.Start();
		PoissonSBDPProfileProcess::Move(1,1,10);
		propchrono.Stop();

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

	// virtual double GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual void SlaveMixMove()	{
		PoissonSBDPProfileProcess::SlaveMixMove();
	}

	void Create(int inNsite, int Nstate, int nratecat, int dc);
	void CreateSuffStat();
	void DeleteSuffStat();

	void GlobalUpdateParameters();
	void SlaveUpdateParameters();

	void GlobalUpdateSiteProfileSuffStat();
	void SlaveUpdateSiteProfileSuffStat();

	int* GetSiteProfileSuffStatCount(int site)	{
		return siteprofilesuffstatcount[site];
	}

	int* allocsiteprofilesuffstatcount;
	int** siteprofilesuffstatcount;

	int* tmpallocsiteprofilesuffstatcount;
	int** tmpsiteprofilesuffstatcount;
	int* tmprrsuffstatcount;

};

#endif
