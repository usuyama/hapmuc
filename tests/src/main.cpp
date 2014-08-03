#include "gtest/gtest.h"

#include <iostream>
#include <list>

#include "IRead.h"
#include "SamRead.h"
#include "log.h"
#include "Parameters.h"
#include "BamReader.h"
#include "HaplotypeBuilder.h"
#include "ProfileHMM.h"
#include "Alignment.h"
#include "CandidateWindow.h"
#include "MutationModel.hpp"
#include "MutationCaller.h"
#include "PairedRead.h"
#include "ProfileHMMUtils.h"

using namespace std;
Parameters params = Parameters();
string data_path; // path to data folder.
faidx_t *fai = NULL; // reference genome

// The fixture for testing class.
class HapMuCTest : public ::testing::Test {
protected:
    // You can remove any or all of the following functions if its body
    // is empty.
    
    HapMuCTest() {
        // You can do set-up work for each test here.
        params.refFileName = data_path + "/random_ref.fasta";
        fai = fai_load(params.refFileName.c_str());
    }
    
    virtual ~HapMuCTest() {
        // You can do clean-up work that doesn't throw exceptions here.
        fai_destroy(fai);
    }
    
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:
    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }
    
    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }
    
    // Objects declared here can be used by all tests in the test case for Foo.
};

// readPair should map to h1
void mapToH1(const Haplotype &h0, const Haplotype &h1, const IRead *read) {
    if (typeid(*read) != typeid(PairedRead)) {
        LOG(logWARNING) << "should be PairedRead" << endl;
        return;
    }
    const PairedRead *readPair = dynamic_cast<const PairedRead *>(read);
    //readPair->first->printVariants();
    //readPair->second->printVariants();
    {
        ProfileHMM pHMM0 = ProfileHMMUtils::genProfileHMM(h0, readPair->first);
        Alignment a0 = pHMM0.viterbi();
        ProfileHMM pHMM1 = ProfileHMMUtils::genProfileHMM(h1, readPair->first);
        Alignment a1 = pHMM1.viterbi();
        a0.print();
        a1.print();
        EXPECT_GT(a1.likelihood, a0.likelihood) << "first " << read->getSeqName();
    }
    {
        ProfileHMM pHMM0 = ProfileHMMUtils::genProfileHMM(h0, readPair->second);
        Alignment a0 = pHMM0.viterbi();
        ProfileHMM pHMM1 = ProfileHMMUtils::genProfileHMM(h1, readPair->second);
        Alignment a1 = pHMM1.viterbi();
        EXPECT_EQ(a0.likelihood, a1.likelihood) << "second " << read->getSeqName();
    }
    {
        vector<Haplotype> haps;
        haps.push_back(h0);
        haps.push_back(h1);

        double lik0 = ProfileHMMUtils::calc(Haplotype::getAllVariants(haps), h0, read);
        double lik1 = ProfileHMMUtils::calc(Haplotype::getAllVariants(haps), h1, read);
        EXPECT_GT(lik1, lik0) << "paired " << read->getSeqName();
    }
}

// readPair should map to h2
void mapToH2(const Haplotype &h0, const Haplotype &h1, const Haplotype &h2, const IRead *read) {
    if (typeid(*read) != typeid(PairedRead)) {
        LOG(logWARNING) << "should be PairedRead" << endl;
        return;
    }
     const PairedRead *readPair = dynamic_cast<const PairedRead *>(read);
    {
        ProfileHMM pHMM0 = ProfileHMMUtils::genProfileHMM(h0, readPair->first);
        Alignment a0 = pHMM0.viterbi();
        ProfileHMM pHMM1 = ProfileHMMUtils::genProfileHMM(h1, readPair->first);
        Alignment a1 = pHMM1.viterbi();
        ProfileHMM pHMM2 = ProfileHMMUtils::genProfileHMM(h2, readPair->first);
        Alignment a2 = pHMM1.viterbi();
        EXPECT_GT(a2.likelihood, a1.likelihood) << "first " << read->getSeqName();
        EXPECT_GT(a2.likelihood, a0.likelihood) << "first " << read->getSeqName();
    }
    {
        ProfileHMM pHMM0 = ProfileHMMUtils::genProfileHMM(h0, readPair->second);
        Alignment a0 = pHMM0.viterbi();
        ProfileHMM pHMM1 = ProfileHMMUtils::genProfileHMM(h1, readPair->second);
        Alignment a1 = pHMM1.viterbi();
        ProfileHMM pHMM2 = ProfileHMMUtils::genProfileHMM(h2, readPair->second);
        Alignment a2 = pHMM1.viterbi();
        EXPECT_EQ(a2.likelihood, a1.likelihood) << "second " << read->getSeqName();
        EXPECT_EQ(a1.likelihood, a0.likelihood) << "second " << read->getSeqName();
    }
    {
        vector<Haplotype> haps;
        haps.push_back(h0);
        haps.push_back(h1);
        haps.push_back(h2);
        double lik0 = ProfileHMMUtils::calc(Haplotype::getAllVariants(haps), h0, read);
        double lik1 = ProfileHMMUtils::calc(Haplotype::getAllVariants(haps), h1, read);
        double lik2 = ProfileHMMUtils::calc(Haplotype::getAllVariants(haps), h2, read);
        EXPECT_GT(lik2, lik0) << "paired " << read->getSeqName();
        EXPECT_GT(lik2, lik1) << "paired " << read->getSeqName();
    }
}

IRead *findBySeqName(string name, const vector<IRead *> &reads) {
    for (int i = 0; i < reads.size(); i++) {
        if (reads[i]->getSeqName() == name) {
            return reads[i];
        }
    }
    throw string("sequence with the given name was not found");
}

// av should be included in the alignment of r given h.
void checkIfVariantIsIncluded(Haplotype h, IRead *read, Variant &av) {
    if (typeid(*read) != typeid(PairedRead)) {
        LOG(logWARNING) << "should be PairedRead" << endl;
        return;
    }
    const PairedRead *readPair = dynamic_cast<const PairedRead *>(read);
    ProfileHMM pHMM = ProfileHMMUtils::genProfileHMM(h, readPair->first);
    Alignment a = pHMM.viterbi();
    vector<Variant> variants = a.getAllVariants();
    vector<Variant>::iterator it = variants.begin();
    bool found = false;
    while (it != variants.end()) {
        if (av.originalString == it->originalString && av.startInGenome == it->startInGenome) {
            found = true;
        }
        it++;
    }
    EXPECT_EQ(true, found);
}

#include "sim_test.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        std::cerr << "please specify the path of data folder" << std::endl;
        exit(1);
    }
    data_path = string(argv[1]);
    return RUN_ALL_TESTS();
}

