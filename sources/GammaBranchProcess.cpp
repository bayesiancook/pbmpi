
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GammaBranchProcess.h"
#include "Random.h"

void GammaBranchProcess::Create(Tree* intree, double inalpha, double inbeta)  {
    BranchProcess::Create(intree);
    branchalpha = inalpha;
    branchbeta = inbeta;
    // RecursiveSampleLength(GetRoot());

    if (! branchempalpha)   {
        branchempalpha = new double[GetNbranch()];
        branchempbeta = new double[GetNbranch()];
        for (int j=0; j<GetNbranch(); j++)  {
            branchempalpha[j] = 1.0;
            branchempbeta[j] = 1.0;
        }
    }
}

void GammaBranchProcess::Delete()   {
    if (branchempalpha) {
        delete[] branchempalpha;
        delete[] branchempbeta;
        branchempalpha = 0;
    }
}

void GammaBranchProcess::ToStream(ostream& os)	{

	SetNamesFromLengths();
	tree->ToStream(os);
	os << branchalpha << '\n';
	os << branchbeta << '\n';
	/*
	for (int j=0; j<GetNbranch(); j++)	{
		os << blarray[j] << '\t';
	}
	os << '\n';
	*/
}

void GammaBranchProcess::ToStreamWithLengths(ostream& os, const Link* from)	{

	if (from->isLeaf())	{
		os << from->GetNode()->GetName();
	}
	else	{
		os << "(";
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStreamWithLengths(os, link->Out());
			if (link->Next() != from)	{
				os << ",";
			}
		}
		os << ")";
	}
	if (! from->isRoot())	{
		os << ":" << blarray[from->GetBranch()->GetIndex()];
	}
}


void GammaBranchProcess::FromStream(istream& is)	{

	tree->ReadFromStream(is);
	tree->RegisterWith(tree->GetTaxonSet());
	SetLengthsFromNames();
	branchalpha = -1;
	branchbeta = -1;
	is >> branchalpha;
	is >> branchbeta;
	/*
	for (int j=0; j<GetNbranch(); j++)	{
		is >> blarray[j];
	}
	*/
}
	
double GammaBranchProcess::LogBranchLengthPrior(const Branch* branch)	{
	int index = branch->GetIndex();
    double a = lengthfrac*branchalpha + (1-lengthfrac)*branchempalpha[index];
    double b = lengthfrac*branchbeta + (1-lengthfrac)*branchempbeta[index];
	return a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(blarray[index]) - b*blarray[index];
	// return branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha-1) * log(blarray[index]) - branchbeta * blarray[index];
}

void GammaBranchProcess::SampleLength(const Branch* branch)	{
	int index = branch->GetIndex();
    double a = lengthfrac*branchalpha + (1-lengthfrac)*branchempalpha[index];
    double b = lengthfrac*branchbeta + (1-lengthfrac)*branchempbeta[index];
	blarray[index] = rnd::GetRandom().Gamma(a,b);
	// blarray[index] = rnd::GetRandom().Gamma(branchalpha,branchbeta);
}
	
void GammaBranchProcess::SampleLength()	{
	cerr << "sample length\n";
	exit(1);
	branchalpha = rnd::GetRandom().sExpo();
	branchbeta = rnd::GetRandom().sExpo();
	blarray[0] = 0;
    RecursiveSampleLength(GetRoot());
}

void GammaBranchProcess::PriorSampleLength()    {
    // fixed to 1.0 in all models
	// branchalpha = rnd::GetRandom().sExpo();
	if (betaprior == 1)	{
        double u = 8 * (rnd::GetRandom().Uniform() - 4.0);
        branchbeta = exp( u * log(10.0) );
    }
    else    {
        branchbeta = 10 * rnd::GetRandom().sExpo();
    }
	blarray[0] = 0;
    RecursiveSampleLength(GetRoot());
}

double GammaBranchProcess::LogHyperPrior()	{
	double total = -branchalpha;
	if (betaprior == 1)	{
		total -= log(branchbeta);
		if ((branchbeta < 1e-4) || (branchbeta > 1e4))	{
			total = -std::numeric_limits<double>::infinity();
			// total -= 1.0 / 0;
		}
	}
	else	{
		total -= 0.1 * branchbeta;
	}
	return total;
}

double GammaBranchProcess::Move(double tuning, int nrep)	{
	double total = MoveLength();
	total += MoveBranchBeta(tuning,nrep);
	// cerr << "in bl move : " << GetTotalLength() << '\n';
	return total;
}

double GammaBranchProcess::MoveBranchBeta(double tuning, int nrep)	{
	int Naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogLengthPrior() - LengthSuffStatLogProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		branchbeta *= e;
		deltalogprob += LogHyperPrior() + LogLengthPrior() + LengthSuffStatLogProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			Naccepted ++;
		}
		else	{
			branchbeta /= e;
		}
	}
	return ((double) Naccepted) / nrep;
}

double GammaBranchProcess::MoveLength()	{

	GlobalUpdateBranchLengthSuffStat();
	for (int i=1; i<GetNbranch(); i++)	{
        double a = lengthfrac*branchalpha + (1-lengthfrac)*branchempalpha[i];
        double b = lengthfrac*branchbeta + (1-lengthfrac)*branchempbeta[i];
		blarray[i] = rnd::GetRandom().Gamma(a + GetBranchLengthSuffStatCount(i), b + GetBranchLengthSuffStatBeta(i));
		// blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
	}
	return 1.0;
}

double GammaBranchProcess::NonMPIMoveLength()	{

	UpdateBranchLengthSuffStat();
	for (int i=1; i<GetNbranch(); i++)	{
        double a = lengthfrac*branchalpha + (1-lengthfrac)*branchempalpha[i];
        double b = lengthfrac*branchbeta + (1-lengthfrac)*branchempbeta[i];
		blarray[i] = rnd::GetRandom().Gamma(a + GetBranchLengthSuffStatCount(i), b + GetBranchLengthSuffStatBeta(i));
		// blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
	}
	return 1.0;
}

double GammaBranchProcess::NonMPIMove(double tuning, int nrep)	{
	NonMPIMoveLength();
	MoveBranchBeta(tuning,nrep);
    return 1.0;
}

