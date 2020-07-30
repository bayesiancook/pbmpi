
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "Chrono.h"

void Chrono::Reset()	{
	TotalTime = 0;
	N = 0;
}

void Chrono::Start()	{
    time1 = std::chrono::steady_clock::now();
}

void Chrono::Stop()	{
    time2 = std::chrono::steady_clock::now();
    TotalTime += std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
}

void Chrono::ToStream(ostream& os)	const {
	os << TotalTime << '\t' << N << '\n';
}

void Chrono::FromStream(istream& is)	{
	is >> TotalTime >> N;
}

int Chrono::operator++()	{
	return N++;
}

double Chrono::GetTime()	{
    return TotalTime;
}

