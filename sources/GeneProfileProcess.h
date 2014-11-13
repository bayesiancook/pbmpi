
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENEPROFILE_H
#define GENEPROFILE_H

#include "PhyloProcess.h"

class GeneProfileProcess : public virtual ProfileProcess {

	public:

	GeneProfileProcess() : hub(0), offset(0) {}
	virtual ~GeneProfileProcess() {}

	virtual double* GetProfile(int site)	{
		return hub->GetProfile(site + offset);
	}

	virtual StateSpace* GetStateSpace()	{
		return hub->GetStateSpace();
	}

	virtual double LogProfilePrior()	{
		return hub->LogProfilePrior();
	}

	virtual void SampleProfile()	{
		// cerr << "in gene sample profile\n";
		// exit(1);
	}

	virtual double GetMeanStationaryEntropy() {return 0;}
	virtual double GetSiteStationaryEntropy(int site) {return 0;}

	virtual void GlobalUpdateParameters() {}
	virtual void SlaveUpdateParameters() {}

	virtual double ProfileSuffStatLogProb() {}

	virtual void SetMixtureParameters() = 0;

	protected:

	virtual void Create(int nsite, int dim, ProfileProcess* inhub, int inoffset)	{
		ProfileProcess::Create(nsite,dim);
		if (! hub)	{
			hub = inhub;
			offset = inoffset;
		}
		else	{
			cerr << "hub already specified\n";
			exit(1);
		}
	}

	virtual void Delete()	{
		ProfileProcess::Delete();
	}

	// UpdateParameters
	// Send back suff stats

	ProfileProcess* hub;
	int offset;

};


class GenePhyloProcess : public virtual GeneProfileProcess, public virtual PhyloProcess	{


	public:

	GenePhyloProcess() {}
	virtual ~GenePhyloProcess() {}

	void SetGibbsFactor(double inf)	{
		gibbsfactor = inf;
	}

	double GetGibbsFactor()	{
		return gibbsfactor;
	}

	virtual void GlobalUpdateSiteProfileSuffStat() {}
	virtual void SlaveUpdateSiteProfileSuffStat() {}

	virtual StateSpace* GetStateSpace()	{
		return GeneProfileProcess::GetStateSpace();
	}

	void TraceHeader(ostream& os)	{
	}

	void Trace(ostream& os)	{
	}

	void SaveTree(string name)	{
		ofstream os((name + datafile + ".treelist").c_str(),ios_base::app);
		GetLengthTree()->ToStream(os);
	}

	void ToStreamHeader(ostream& os)	{
	}

	virtual void ToStream(ostream& os) {}
	virtual void FromStream(istream& is) {}

	virtual void Delete()	{
		PhyloProcess::Delete();
		GeneProfileProcess::Delete();
	}

	double gibbsfactor;

};
#endif
