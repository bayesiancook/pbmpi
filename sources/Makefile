CC=mpic++
CPPFLAGS= -Wall -O3 -std=c++11
LDFLAGS= -O3
SRCS=  TaxonSet.cpp Tree.cpp Random.cpp SequenceAlignment.cpp CodonSequenceAlignment.cpp \
	StateSpace.cpp CodonStateSpace.cpp ZippedSequenceAlignment.cpp SubMatrix.cpp \
	GTRSubMatrix.cpp CodonSubMatrix.cpp linalg.cpp Chrono.cpp BranchProcess.cpp \
	GammaBranchProcess.cpp RateProcess.cpp DGamRateProcess.cpp ProfileProcess.cpp \
	OneProfileProcess.cpp MatrixProfileProcess.cpp MatrixOneProfileProcess.cpp \
	GTRProfileProcess.cpp ExpoConjugateGTRProfileProcess.cpp \
	ExpoConjugateGTRDPProfileProcess.cpp ExpoConjugateGTRSBDPProfileProcess.cpp \
	ExpoConjugateGTRFiniteProfileProcess.cpp \
	MixtureProfileProcess.cpp \
	PoissonMixtureProfileProcess.cpp MatrixMixtureProfileProcess.cpp \
	GTRMixtureProfileProcess.cpp ExpoConjugateGTRMixtureProfileProcess.cpp \
	GeneralPathSuffStatMatrixMixtureProfileProcess.cpp \
	DPProfileProcess.cpp SBDPProfileProcess.cpp \
	PoissonDPProfileProcess.cpp PoissonSBDPProfileProcess.cpp MatrixDPProfileProcess.cpp \
	MatrixSBDPProfileProcess.cpp FiniteProfileProcess.cpp PoissonFiniteProfileProcess.cpp \
	MatrixFiniteProfileProcess.cpp CodonMutSelProfileProcess.cpp \
	AACodonMutSelProfileProcess.cpp SubstitutionProcess.cpp Propagate.cpp \
	PoissonSubstitutionProcess.cpp MatrixSubstitutionProcess.cpp \
	ExpoConjugateGTRSubstitutionProcess.cpp GeneralPathSuffStatMatrixSubstitutionProcess.cpp \
	PhyloProcess.cpp PoissonPhyloProcess.cpp MatrixPhyloProcess.cpp Stepping.cpp \
	ExpoConjugateGTRPhyloProcess.cpp GeneralPathSuffStatMatrixPhyloProcess.cpp \
	RASCATGammaPhyloProcess.cpp RASCATSBDPGammaPhyloProcess.cpp \
	RASCATFiniteGammaPhyloProcess.cpp RASCATGTRFiniteGammaPhyloProcess.cpp \
	RASCATGTRDPGammaPhyloProcess.cpp RASCATGTRSBDPGammaPhyloProcess.cpp \
	AACodonMutSelFinitePhyloProcess.cpp CodonMutSelFinitePhyloProcess.cpp\
	CodonMutSelSBDPPhyloProcess.cpp \
	AACodonMutSelSBDPPhyloProcess.cpp \
	Bipartition.cpp BipartitionList.cpp Consensus.cpp TaxaParameters.cpp PBTree.cpp TreeList.cpp PolyNode.cpp correl.cpp correlation.cpp NNI.cpp


OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= pb_mpi readpb_mpi tracecomp bpcomp cvrep
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

$(PROGSDIR)/sitecv : SiteCV.o $(OBJS)
	$(CC) SiteCV.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/cvrep : CVrep.o $(OBJS)
	$(CC) CVrep.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/jackknife: JackKnife.o $(OBJS)
	$(CC) JackKnife.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/readpb_mpi: ReadPB.o $(OBJS)
	$(CC) ReadPB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/pb_mpi: PB.o $(OBJS)
	$(CC) PB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/tracecomp: TraceComp.o $(OBJS)
	$(CC) TraceComp.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/bpcomp: BPCompare.o $(OBJS)
	$(CC) BPCompare.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

