
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#pragma once

#include <chrono>
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
	double GetTime();

	int operator++();
	int GetCount();

	void FromStream(std::istream& is);
    void ToStream(std::ostream& os) const;

	private:

    // in milliseconds
    std::chrono::steady_clock::time_point time1;
    std::chrono::steady_clock::time_point time2;
    double TotalTime;
	int N;
};

