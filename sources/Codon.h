
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/



class AAMutSelProfileProcess : public virtual GeneralPathSuffStatMatrixProfileProcess	{

	// y mettre les variables globales (taux de mutation essentiellement)

	// s'inspirer de GTRProfileProcess et GeneralPathSuffStatGTRProfileProcess
};

class AAMutSelDPProfileProcess : public virtual AAMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixDPProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRDPProfileProcess
	
};


class AAMutSelDPSubstitutionProcess : public virtual AAMutSelDPProfileProcess, public virtual UniformRateProcess, public virtual GeneralPathSuffStatMatrixSubstitutionProcess {

	// s'inspirer de GeneralPathSuffStatGTRSubstitutionProcess
	// et GeneralPathSuffStatRASCATGTRSubstitutionProcess


};


class AAMutSelDPPhyloProcess : public virtual AAMutSelDPSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess	{

	// s'inspirer de GeneralPathSuffStatGTRPhyloProcess
	// et GeneralPathSuffStatRASCATGTRPhyloProcess

};

// enfin, le PhyloProcess ainsi construit peut etre instancie dans le main.cpp


