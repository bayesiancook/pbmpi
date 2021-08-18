
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef __PARALLELH
#define __PARALLELH 

#include "mpi.h"

const int TAG1 = 91;

enum MESSAGE {KILL,SCAN,UPDATE_RATE,UPDATE_RRATE,UPDATE_BLENGTH,UPDATE_SRATE,UPDATE_SPROFILE,PARAMETER_DIFFUSION,UNFOLD,COLLAPSE,LIKELIHOOD,RESET,RESETALL,MULTIPLY,SMULTIPLY,INITIALIZE,PROPAGATE,PROPOSE,RESTORE,UPDATE,DETACH,ATTACH,NNI,KNIT,BRANCHPROPAGATE,ROOT,REALLOC_MOVE,PROFILE_MOVE,MIX_MOVE,REALLOC_DONE,GIVEMEMORE,BCAST_TREE,UNCLAMP,SETDATA,SETNODESTATES,CVSCORE,SETTESTDATA,GENE_MOVE,SAMPLE,LENGTH,ALPHA,SAVETREES, LENGTHFACTOR, FROMSTREAM, TOSTREAM, SITELOGL, STEPPINGSITELOGL, FULLSITELOGL, RESTOREDATA, WRITE_MAPPING,NONSYNMAPPING,COUNTMAPPING,SITERATE,SIMULATE,SETRATEPRIOR,SETPROFILEPRIOR,SETROOTPRIOR,STATEPOSTPROBS,SITELOGLCUTOFF,SITELOGCV, PREPARESTEPPING, SETSTEPPINGFRAC, EMPIRICALFRAC, EMPIRICALPRIOR, CREATESITE, DELETESITE};

#endif


