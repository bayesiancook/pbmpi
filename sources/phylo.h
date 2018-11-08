
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include <sstream>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

using namespace std;

typedef int Boolean;

	enum	Switch		{No = 0, Yes = 1};

	enum 	Prior 		{Flat = 0, PowerLaw = 1, GammaInv = 2, Dirichlet = 3 , MultiGamma = 4, Exponential = 5, OnOff = 6, GammaDistributed = 7, Dirac = 8, GammaMixture = 9, AlphaBased = 10, BetaBased = 11, Trapezoidal = 12, DirichletProcess = 13, Normal = 14, Cauchy=15};

	const int CheckLevel = 0;
	const int verbose  = 0;

#include "StringStreamUtils.h"

#include "TaxaParameters.h"
#include "Bipartition.h"
#include "BipartitionList.h"
#include "TreeList.h"

#include "PolyNode.h"
#include "PBTree.h"
#include "AmphiNode.h"
#include "Consensus.h"

