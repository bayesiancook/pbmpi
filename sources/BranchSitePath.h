
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef SITEPATH_H
#define SITEPATH_H

#include <string>
#include <map>
#include "Tree.h"
#include "SubMatrix.h"

class Plink	{

	friend class BranchSitePath;

	public:

				Plink();
				Plink(int instate, double inrel_time);
				~Plink();
	Plink* 			Prev();
	Plink* 			Next();

	bool 			IsFirst();
	bool 			IsLast();

	void 			Splice();
	void 			Insert(Plink* link);

	void 			SetState(int instate);
	int 			GetState();

	void 			SetRelativeTime(double inrel_time);
	double 			GetRelativeTime();

	private:

	Plink* next;
	Plink* prev;
	
	int state;
	double rel_time;

};

class BranchSitePath  {

	public:

				// special constructor: for Poisson processes (for which the detailed path does not matter,
				// only the total number of substitutions and the final state)
				BranchSitePath(int incount, int instate);

				BranchSitePath(int instate = -1);
	virtual			~BranchSitePath();

	Plink* 			Init();
	Plink*  		Last();
	int 			GetNsub();
	int 			GetInitState();
	int 			GetFinalState();
	double 			GetRelativeTime(Plink* link) {return link->GetRelativeTime();}

	void 			Reset(int state);
	void 			Append(int instate, double reltimelength);

	void 			Print(ostream& os)	{
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			os << state << "  :  " << GetRelativeTime(link);
			link = link->Next();
		}
		os << "\t:::\t" << GetNsub();
		os << '\n';
	}


	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	void AddRateSuffStat(int& count, double& beta, double factor, const double* rr, const double* stat, int nstate)	{
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			double tmp = 0;
			for (int k=0; k<nstate; k++)	{
				if (k != state)	{
					tmp += rr[rrindex(state,k,nstate)] * stat[k];
				}
			}
			beta += GetRelativeTime(link) * factor * tmp;
			if (link != last)	{
				count ++;
			}
			link = link->Next();
		}
	}

	void AddProfileSuffStat(int* count, double* beta, double factor, const double* rr, int nstate)	{
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			for (int k=0; k<nstate; k++)	{
				if (k!=state)	{
					beta[k] += GetRelativeTime(link) * factor * rr[rrindex(state,k,nstate)];
				}
			}
			if (link != last)	{
				count[link->Next()->GetState()] ++;
			}
			link = link->Next();
		}
	}

	void AddRRSuffStat(int* count, double* beta, double factor, const double* stat, int nstate)	{
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			for (int k=0; k<nstate; k++)	{
				if (k!=state)	{
					beta[rrindex(state,k,nstate)] += GetRelativeTime(link) * factor * stat[k];
				}
			}
			if (link != last)	{
				count[rrindex(state,link->Next()->GetState(),nstate)] ++;
			}
			link = link->Next();
		}
	}

	void AddGeneralPathRateSuffStat(int& count, double& beta, double factor, SubMatrix* mat)	{
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			beta -= GetRelativeTime(link) * factor * (*mat)(state,state);
			if (link != last)	{
				count ++;
			}
			link = link->Next();
		}
	}

	void AddGeneralPathSuffStat(map<pair<int,int>,int>& paircount, map<int,double>& waitingtime, double factor)	{
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			waitingtime[state] += GetRelativeTime(link) * factor;
			if (link != last)	{
				paircount[pair<int,int>(state,link->Next()->GetState())]++;
			}
			link = link->Next();
		}
	}

	Plink* init;
	Plink* last;
	int nsub;

};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//	* Plink
//-------------------------------------------------------------------------

inline Plink::Plink() : next(0), prev(0), state(0), rel_time(0)  {}
inline Plink::Plink(int instate, double inrel_time) : next(0), prev(0), state(instate), rel_time(inrel_time){}
inline Plink::~Plink()	{Splice();}

inline Plink* Plink::Prev() {return prev;}
inline Plink* Plink::Next() {return next;}

inline bool Plink::IsFirst()	{return prev == 0;}
inline bool Plink::IsLast() {return next == 0;}

inline void Plink::Insert(Plink* link)	{
	link->next = next;
	if (next)	{
		next->prev = link;
	}
	link->prev = this;
	next = link;
}

inline void Plink::Splice()	{
	if (prev)	{
		prev->next = next;
	}
	if (next)	{
		next->prev = prev;
	}
	prev = next = 0;
}

inline void Plink::SetState(int instate) {state = instate;}
inline void Plink::SetRelativeTime(double inrel_time) {
	if (std::isnan(rel_time))	{
		cerr << "in Plink::SetRelativeTime: setting to nan\n";
		exit(1);
	}
	rel_time = inrel_time;
}
inline double Plink::GetRelativeTime() {return rel_time;}
inline int Plink::GetState() {return state;}

//-------------------------------------------------------------------------
//	* BranchSitePath
//-------------------------------------------------------------------------

inline Plink* BranchSitePath::Init() {return init;}
inline Plink*  BranchSitePath::Last() {return last;}

inline void BranchSitePath::Append(int instate, double reltimelength)	{
	last->SetRelativeTime(reltimelength);
	Plink* link = new Plink(instate,0);
	last->Insert(link);
	last = link;
	nsub++;
}	

inline int BranchSitePath::GetNsub()	{
	return nsub;
}

inline int BranchSitePath::GetInitState() {return init->GetState();}
inline int BranchSitePath::GetFinalState() {return last->GetState();}


inline BranchSitePath::BranchSitePath(int incount, int infinalstate)	{
	init = last = new Plink;
	init->SetState(infinalstate);
	nsub = incount;
}

inline BranchSitePath::BranchSitePath(int instate)	{
	init = last = new Plink;
	init->SetState(instate);
	nsub = 0;
}

inline BranchSitePath::~BranchSitePath()	{
	Reset(0);
	delete init;
}

inline void BranchSitePath::Reset(int state)	{
	Plink* link=last;
	while (link!=init)	{
		Plink* prev = link->Prev();
		delete link;
		link = prev;
	}
	nsub = 0;
	init->SetState(state);
	init->SetRelativeTime(0);
	last = init;
}

#endif 
