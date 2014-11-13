
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef CHRONO_H
#define CHRONO_H

#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <cstdlib>

using namespace std;

class Chrono	{

	public:

	Chrono(){
		Reset();
	};
	~Chrono() {};
	void Reset();
	void Start();
	void Stop();

	void ToStream(ostream& os);
	void FromStream(istream& is);

	int operator++();

	double GetTime();
	double GetTimePerCount();
	int GetCount();

	private:

	// this is in milli seconds
	double sec1;
	double sec2;
	double milli1;
	double milli2;

	double TotalTime;
	int N;

}
;


#endif		
