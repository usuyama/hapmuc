SAMTOOLDIR=../samtools-0.1.18/
SEQANDIR=seqan_library/
BOOST=/Users/usuyama/bin/boost_1_45_0/
BOOSTLIB=${BOOST}stage/lib/

CXX= g++
CPPFLAGS= -DNDEBUG -D_IOLIB=2 -DMINREADS=2 -DDINDEL
CXXFLAGS= -I$(SAMTOOLDIR) -I$(SEQANDIR) -I${BOOST} -I./ -Wno-deprecated  -O3
LDFLAGS= -L$(SAMTOOLDIR) -lbam -lz -L${BOOSTLIB} -lboost_program_options
SRCSDINDEL=DInDel.cpp Utils.cpp MutationModel.cpp HapBlock.cpp HaplotypeDistribution.cpp ObservationModelFB.cpp GetCandidates.cpp Faster.cpp Haps.cpp Haps2.cpp EMBasic.cpp EMfor2.cpp MutationCall.cpp
OBJSDINDEL=$(SRCSDINDEL:%.cpp=%.o)
all: dindel
dindel:$(OBJSDINDEL) Utils.hpp Read.hpp DInDel.hpp HapBlock.hpp Haplotype.hpp HaplotypeDistribution.hpp MyBam.hpp GetCandidates.hpp Variant.hpp EMBasic.hpp Fasta.hpp OutputData.hpp MLAlignment.hpp ObservationModelSeqAn.hpp VariantFile.hpp ReadIndelErrorModel.hpp Library.hpp Faster.hpp EMfor2.hpp Haps.hpp Haps2.hpp MutationCall.hpp
	$(CXX) -o $@ $(CXXFLAGS) $(DINDELFLAGS) $(OBJSDINDEL) $(LDFLAGS)

test: dindel
	./dindel --opl --doPooled --bamFile /Users/usuyama/data/test/test.bam --ref /Users/usuyama/data/hg19/hg19.fasta --varFile sample.realign_windows.1.txt --libFile sample.dindel_output.libraries.txt --outputFile sample.dindel_stage2_output_windows.1 > log

mc: dindel
	./dindel --analysis mutationCall --tumorBamFiles testTumorBams --normalBamFiles testNormalBams  --ref /Users/usuyama/data/hg19/hg19.fasta --varFile sample.realign_windows.1.txt --libFile sample.dindel_output.libraries.txt --outputFile mutationCall.txt > logmc

mc2: dindel
	./dindel --analysis mutationCall --tumorBamFiles testTumorBams2 --normalBamFiles testNormalBams2  --ref /Users/usuyama/data/hg19/hg19.fasta --varFile window2.txt --libFile sample.dindel_output.libraries.txt --outputFile mutationCall2.txt > logmc2

clean:
	rm -f $(OBJSDINDEL) $(OBJSCOMPAREVARIANTS)  $(OBJSMAKEGLF)
