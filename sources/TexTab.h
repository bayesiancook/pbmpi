
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef TEXTAB_H
#define TEXTAB_H

#include <sstream>
#include <list>

/*	
string textabentry(double mean, double min, double var, int sign = 2, double pp = 0.5, double c1 = 0.95, double c2 = 0.975)	{

	double pow = 1;
	for (int i=0l i<sign)	{
		pow *= 10;
	}

	int mean1 = (int) mean;
	int mean2 = (int) (pow * (mean - mean1) + 0.5);
	
	int min1 = (int) mean;
	int min2 = (int) (pow * (min - min1) + 0.5);
	
	int max1 = (int) max;
	int max2 = (int) (pow * (max - max1) + 0.5);
	

	ostringstream s;
	s << mean1 << "&" << mean2;
	if ((pp > c2) || (pp < (1-c2)))	{
		s << "${}^{**}$";
	}
	else if ((pp > c1) || (pp < (1-c1)))	{
		s << "${}^{*}$";
	}
	s << "\t&\t" <<;
	s << min1 << "&" << min2;
	s << "\t&\t" <<;
	s << max1 << "&" << max2;

	return s.str();
	
}
*/


inline string textabentry(list<double>& l, bool withmean = true, bool withsign = true, bool withci = true, int sign_digit = 2, double c = 0.95, double c1 = 0.95, double c2 = 0.975)	{

	double mean = 0;
	double pp = 0;

	for (list<double>::iterator i=l.begin(); i!=l.end(); i++)	{
		mean += *i;
		if (*i > 0)	{
			pp++;
		}
	}
	mean /= l.size();
	pp /= l.size();

	l.sort();
	int n = ((int) (((double) l.size()) * (1-c)));
	list<double>::const_iterator i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	double min = *i;
	n = ((int) (((double) l.size()) * c));
	i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	double max = *i;

	double pow = 1;
	for (int i=0; i<sign_digit; i++)	{
		pow *= 10;
	}

	ostringstream s;

	if (withmean)	{
		int meansign = 1;
		if (mean)	{
			meansign = (int) (mean / fabs(mean));
		}
		mean = fabs(mean);
		int mean1 = (int) mean;
		int mean2 = (int) (pow * (mean - mean1) + 0.5);
		if (mean2 == pow)	{
			mean2 = 0;
			mean1++;
		}

		if (meansign == -1)	{
			s << "-";
		}
		s << mean1 << "&";
		// s.precision(sign_digit);
		s << setw(sign_digit) << setfill('0') << mean2;
		// s.precision(0);

		if (withsign)	{
			if ((pp > c2) || (pp < (1-c2)))	{
				s << "${}^{**}$";
			}
			else if ((pp > c1) || (pp < (1-c1)))	{
				s << "${}^{*}$";
			}
		}
	}
	
	if (withci)	{
		s << "&\t";

		int minsign = 1;
		if (min)	{
			minsign = (int) (min / fabs(min));
		}
		min = fabs(min);

		int maxsign = 1;
		if (max)	{
			maxsign = (int) (max / fabs(max));
		}
		max = fabs(max);

		int min1 = (int) min;
		int min2 = (int) (pow * (min - min1) + 0.5);
		if (min2 == pow)	{
			min2 = 0;
			min1 ++;
		}
		
		int max1 = (int) max;
		int max2 = (int) (pow * (max - max1) + 0.5);
		if (max2 == pow)	{
			max2 = 0;
			max1 ++;
		}

		if (minsign == -1)	{
			s << "-";
		}
		s << min1 << "&";
		// s.precision(sign_digit);
		s << setw(sign_digit) << setfill('0') << min2;
		// s << min2;
		// s.precision(0);
		s << "&\t";

		if (maxsign == -1)	{
			s << "-";
		}
		s << max1 << "&";
		// s.precision(sign_digit);
		s << setw(sign_digit) << setfill('0') << max2;
		// s << max2;
		// s.precision(0);
		s << "&\t";
	}

	return s.str();
	
}

inline void printCI(list<double>& l, ostream& os, double c = 0.95)	{

	double mean = 0;
	double pp = 0;

	for (list<double>::iterator i=l.begin(); i!=l.end(); i++)	{
		mean += *i;
		if (*i > 0)	{
			pp++;
		}
	}
	mean /= l.size();
	pp /= l.size();

	l.sort();
	int n = ((int) (((double) l.size()) * (1-c)));
	list<double>::const_iterator i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	double min = *i;
	n = ((int) (((double) l.size()) * c));
	i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	double max = *i;

	os << mean << " (" << min << "," << max << ")";
}


#endif

