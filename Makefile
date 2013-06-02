.PHONY: all clean test
SAMTOOLDIR=../samtools-0.1.18/
SEQANDIR=seqan_library/
BOOST=/Users/usuyama/bin/boost_1_45_0/
BOOSTLIB=${BOOST}stage/lib/

CXX= g++
CPPFLAGS= -g -DNDEBUG -D_IOLIB=2 -DMINREADS=2 -DHAPMUC# -DMMTEST #for mutation_model_test
CXXFLAGS= -g -I$(SAMTOOLDIR) -I$(SEQANDIR) -I${BOOST} -I./ -Wno-deprecated  -O3
LDFLAGS= -L$(SAMTOOLDIR) -lbam -lz -L${BOOSTLIB} -lboost_program_options
SRCSDINDEL=HapMuC.cpp Utils.cpp MutationModel.cpp HapBlock.cpp HaplotypeDistribution.cpp ObservationModelFB.cpp Haps.cpp Haps2.cpp EMBasic.cpp EMfor2.cpp MutationCall.cpp
OBJSDINDEL=$(SRCSDINDEL:%.cpp=%.o)
all: hapmuc
hapmuc:$(OBJSDINDEL) 
	$(CXX) -o $@ $(CXXFLAGS) $(DINDELFLAGS) $(OBJSDINDEL) $(LDFLAGS)

OBJSTEST=Utils.o MutationModel.o HapBlock.o HaplotypeDistribution.o ObservationModelFB.o Haps.o Haps2.o EMBasic.o EMfor2.o MutationCall.o mutation_model_test.o
mutation_model_test: $(OBJSTEST)
	$(CXX) -o mutation_model_test $(CXXFLAGS) $(OBJSTEST) $(LDFLAGS)

clean:
	rm -f $(OBJSDINDEL) $(OBJSCOMPAREVARIANTS)  $(OBJSMAKEGLF)

test: hapmuc
	@sh hapmuc_testcases/test.sh 1
	@sh hapmuc_testcases/test.sh 2
	@sh hapmuc_testcases/test.sh 3
	@sh hapmuc_testcases/test.sh 4
	@sh hapmuc_testcases/test.sh 5
