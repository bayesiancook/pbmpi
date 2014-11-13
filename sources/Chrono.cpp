
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
	// clock_gettime(CLOCK_REALTIME, &nano1);

	struct timeval tod;
	gettimeofday(&tod,NULL);
	sec1 = ((double) tod.tv_sec);
	milli1 = ((double) tod.tv_usec) / 1000;
}

void Chrono::Stop()	{
	/*
	clock_gettime(CLOCK_REALTIME, &nano2);
	double t1 = ((double) (nano2.tv_sec))  - ((double) (nano1.tv_sec));
	double t2 = (nano2.tv_nsec - nano1.tv_nsec) * (1e-9);
	double duration = t1 + t2;
	*/

	struct timeval tod;
	gettimeofday(&tod,NULL);
	sec2 = ((double) tod.tv_sec);
	milli2 = ((double) tod.tv_usec) / 1000;
	double duration = 1000*(sec2 - sec1) + milli2 - milli1;

	TotalTime += duration;
}

void Chrono::ToStream(ostream& os)	{
	os << TotalTime << '\t' << N << '\n';
}

void Chrono::FromStream(istream& is)	{
	is >> TotalTime >> N;
}

int Chrono::operator++()	{
	return N++;
}

double Chrono::GetTime()	{
	double tmp = ((double) (((long) (1000 * TotalTime)) / 1000));
	if (tmp < 0)	{
		cerr << "error : negative time : " << TotalTime << '\t' << 1000 * TotalTime << '\t' << (int) (1000 * TotalTime) << '\n';
		exit(1);
	}
}

double Chrono::GetTimePerCount()	{
	return TotalTime / N;
}

int Chrono::GetCount()	{
	return N;
}

