.PHONY: all clean test
SAMTOOLDIR=libs/samtools-0.1.19/
SEQANDIR=libs/seqan_library/
BOOST=libs/boost_1_45_0_subset/

CXX= g++
CPPFLAGS= -D_IOLIB=2 -DMINREADS=2 #-DLOGDEBUG -DDEBUGREADS #-DNDEBUG -DMMTEST
CXXFLAGS= -I$(SAMTOOLDIR) -I$(SEQANDIR) -I${BOOST} -I./ -Wno-deprecated -O3
LDFLAGS= -L$(SAMTOOLDIR) -lbam -lz -lpthread
SRCSHAPMUC=HapMuC.cpp Utils.cpp MutationModel.cpp HapBlock.cpp HaplotypeDistribution.cpp ObservationModelFB.cpp Haps.cpp Haps2.cpp EMBasic.cpp MutationCall.cpp
OBJSHAPMUC=$(SRCSHAPMUC:%.cpp=%.o)

hapmuc:$(OBJSHAPMUC) EMBasic.hpp Fasta.hpp HapBlock.hpp HapMuC.hpp Haplotype.hpp HaplotypeDistribution.hpp Haps.hpp Haps2.hpp MLAlignment.hpp MutationCall.hpp MutationModel.hpp MyBam.hpp ObservationModel.hpp ObservationModelFB.hpp ObservationModelSeqAn.hpp OutputData.hpp Read.hpp ReadIndelErrorModel.hpp StringHash.hpp Utils.hpp Variant.hpp VariantFile.hpp cmdline.h log.h
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) $(OBJSHAPMUC) $(LDFLAGS)

dependencies:
	cd libs; make

clean:
	rm -f $(OBJSHAPMUC)
