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


TEST_F(HapMuCTest, HaplotypeBuilderBasic) {
    vector<Variant> variants = vector<Variant>();
    int pos = 1436;
    string chr = "chrA";
    variants.push_back(Variant("T=>G", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h0 = hb.build(pos, vector<Variant>(), chr, 1400, 1500);
    EXPECT_EQ("CAGGCGCGATGTGATAGGACATTTATATCCCGCGCGTCCTAGGTCTTGAATAATACGCAGCTGGGCCCAAAGCTGTCTGATCTATGTTGGATTCAGAACT", h0.seq);
    Haplotype h1 = hb.build(pos, variants, chr, 1400, 1500);
    EXPECT_EQ("CAGGCGCGATGTGATAGGACATTTATATCCCGCGCGGCCTAGGTCTTGAATAATACGCAGCTGGGCCCAAAGCTGTCTGATCTATGTTGGATTCAGAACT", h1.seq);
}

TEST_F(HapMuCTest, HaplotypeBuilderDeletion) {
    int pos = 1826; // to 1829
    string chr = "chrA";
    vector<Variant> variants = vector<Variant>();
    variants.push_back(Variant("-GAC", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h1 = hb.build(pos, variants, chr, 1800, 1850);
    EXPECT_EQ("CGCACTGTCTATAAGCGCATATAGAAGGACCGTTATCCTCATTCAGG", h1.seq);
}

TEST_F(HapMuCTest, HaplotypeBuilderInsertion) {
    int pos = 1328; // zero-based
    string chr = "chrB";
    vector<Variant> variants = vector<Variant>();
    variants.push_back(Variant("+GTC", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h1 = hb.build(pos, variants, chr, 1300, 1350);
    EXPECT_EQ("TGAGCAGGCCGTCTCGCCGTATCAGCAGGTCTCGCAGTAACTGAACTTCCATT", h1.seq);
}

TEST_F(HapMuCTest, ProfileHMM) {
    MyBam bam = MyBam(data_path + "/1/tumor.bam");
    BamReader tumorBamReader = BamReader(&bam, fai);
    int pos = 1436;
    string chr = "chrA";
    vector<IRead *> tumorReadPairs = tumorBamReader.getPairedReads(chr, pos, pos - 1000, pos + 1000);
    vector<Variant> variants = vector<Variant>();
    variants.push_back(Variant("T=>G", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h0 = hb.build(pos, vector<Variant>(), chr, pos - 800, pos + 800 + 100);
    Haplotype h1 = hb.build(pos, variants, chr, pos - 800, pos + 800 + 100);
    mapToH1(h0, h1, findBySeqName("33", tumorReadPairs));
    mapToH1(h1, h0, findBySeqName("49", tumorReadPairs));
    mapToH1(h0, h1, findBySeqName("25", tumorReadPairs));
    mapToH1(h1, h0, findBySeqName("42", tumorReadPairs));
    checkIfVariantIsIncluded(h0, findBySeqName("33", tumorReadPairs), variants[0]);
    checkIfVariantIsIncluded(h0, findBySeqName("25", tumorReadPairs), variants[0]);
    // check variants
}

TEST_F(HapMuCTest, SoftClipping) {
    MyBam bam = MyBam(data_path + "/1/tumor.bam");
    BamReader tumorBamReader = BamReader(&bam, fai);
    int pos = 1436;
    string chr = "chrA";
    vector<IRead *> tumorReadPairs = tumorBamReader.getPairedReads(chr, pos, pos - 1000, pos + 1000);
    vector<Variant> variants = vector<Variant>();
    variants.push_back(Variant("T=>G", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h0 = hb.build(pos, vector<Variant>(), chr, pos - 800, pos + 800 + 100);
    Haplotype h1 = hb.build(pos, variants, chr, pos - 800, pos + 800 + 100);
    IRead *softClippedIRead = findBySeqName("8", tumorReadPairs);
    const PairedRead *softClippedReadPair = dynamic_cast<const PairedRead *>(softClippedIRead);
    EXPECT_EQ(1368, softClippedReadPair->first->leftMostPos); // NOTE: 1-based positions in sam/bam files
    mapToH1(h1, h0, softClippedIRead);
}

TEST_F(HapMuCTest, ProfileHMMDeletion) {
    MyBam bam = MyBam(data_path + "/del/tumor.bam");
    BamReader tumorBamReader = BamReader(&bam, fai);
    int pos = 1826;
    string chr = "chrA";
    vector<IRead *> tumorReadPairs = tumorBamReader.getPairedReads(chr, pos, pos - 1000, pos + 1000);
    vector<Variant> variants = vector<Variant>();
    variants.push_back(Variant("-GAC", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h0 = hb.build(pos, vector<Variant>(), chr, pos - 800, pos + 800 + 100);
    Haplotype h1 = hb.build(pos, variants, chr, pos - 800, pos + 800 + 100);
    mapToH1(h0, h1, findBySeqName("18", tumorReadPairs));
    mapToH1(h1, h0, findBySeqName("16", tumorReadPairs));
    mapToH1(h0, h1, findBySeqName("10", tumorReadPairs));
    mapToH1(h1, h0, findBySeqName("6", tumorReadPairs));
    checkIfVariantIsIncluded(h0, findBySeqName("18", tumorReadPairs), variants[0]);
    checkIfVariantIsIncluded(h0, findBySeqName("10", tumorReadPairs), variants[0]);
}

TEST_F(HapMuCTest, ProfileHMMInsertion) {
    MyBam bam = MyBam(data_path + "/ins/tumor.bam");
    BamReader tumorBamReader = BamReader(&bam, fai);
    int pos = 1328;
    string chr = "chrB";
    vector<IRead *> tumorReadPairs = tumorBamReader.getPairedReads(chr, pos, pos - 1000, pos + 1000);
    vector<Variant> variants = vector<Variant>();
    variants.push_back(Variant("+GTC", pos));
    HaplotypeBuilder hb = HaplotypeBuilder(fai);
    Haplotype h1 = hb.build(pos, variants, chr, pos - 800, pos + 800 + 100);
    Haplotype h0 = hb.build(pos, vector<Variant>(), chr, pos - 800, pos + 800 + 100);
    mapToH1(h0, h1, findBySeqName("18", tumorReadPairs));
    mapToH1(h1, h0, findBySeqName("36", tumorReadPairs));
    mapToH1(h0, h1, findBySeqName("11", tumorReadPairs));
    mapToH1(h1, h0, findBySeqName("17", tumorReadPairs));
    checkIfVariantIsIncluded(h0, findBySeqName("18", tumorReadPairs), variants[0]);
    checkIfVariantIsIncluded(h0, findBySeqName("11", tumorReadPairs), variants[0]);
}


TEST_F(HapMuCTest, Total1) {
    MyBam tumorBam = MyBam(data_path + "/1/tumor.bam");
    BamReader tumorBamReader = BamReader(&tumorBam, fai);
    MyBam normalBam = MyBam(data_path + "/1/normal.bam");
    BamReader normalBamReader = BamReader(&normalBam, fai);
    MutationCaller caller = MutationCaller(tumorBamReader, normalBamReader, params, fai);
    ifstream varfile((data_path + "/1/windows").c_str());
    string line;
    vector<string> results;
    vector<MutationCallResult> mcr_results;
    results.push_back(MutationCallResult::getHeader());
    while (getline(varfile, line)) {
        try {
            CandidateWindow cw = CandidateWindow(line);
            MutationCallResult result = caller.call(cw);
            mcr_results.push_back(result);
            results.push_back(result.getOutput());
        } catch (std::string &s) {
            LOG(logERROR) << s << std::endl;
        }
    }
    for (int i = 0; i < results.size(); i++) {
        LOG(logDEBUG) << results[i] << endl;
    }
    EXPECT_EQ(results.size(), 4);
    EXPECT_GT(mcr_results[0].getBayesFactor(), mcr_results[0].getBayesFactorBasic());
    EXPECT_GT(mcr_results[0].getBayesFactorBasic(), 1);
    EXPECT_LT(mcr_results[1].getBayesFactor(), mcr_results[1].getBayesFactorBasic());
    EXPECT_GT(mcr_results[1].getBayesFactorBasic(), 1);
    EXPECT_LT(mcr_results[1].getBayesFactor(), 1);
    EXPECT_FALSE(mcr_results[2].hasBayesFactor);
}


TEST_F(HapMuCTest, Total2) {
    MyBam tumorBam = MyBam(data_path + "/two/tumor.bam");
    BamReader tumorBamReader = BamReader(&tumorBam, fai);
    MyBam normalBam = MyBam(data_path +  "/two/normal.bam");
    BamReader normalBamReader = BamReader(&normalBam, fai);
    MutationCaller caller = MutationCaller(tumorBamReader, normalBamReader, params, fai);
    ifstream varfile((data_path + "/two/windows").c_str());
    string line;
    getline(varfile, line);
    CandidateWindow cw = CandidateWindow(line);
    MutationCallResult result = caller.call(cw);
    EXPECT_GT(result.getBayesFactorBasic(), 1);
    EXPECT_FALSE(result.hasBayesFactor);
    getline(varfile, line);
    CandidateWindow cw2 = CandidateWindow(line);
    MutationCallResult result2 = caller.call(cw2);
    EXPECT_LT(result2.getBayesFactor(), 1);
    EXPECT_GT(result2.getBayesFactorBasic(), 1);
}

TEST_F(HapMuCTest, TotalIns) {
    MyBam tumorBam = MyBam(data_path + "/ins/tumor.bam");
    BamReader tumorBamReader = BamReader(&tumorBam, fai);
    MyBam normalBam = MyBam(data_path + "/ins/normal.bam");
    BamReader normalBamReader = BamReader(&normalBam, fai);
    MutationCaller caller = MutationCaller(tumorBamReader, normalBamReader, params, fai);
    ifstream varfile((data_path + "/ins/windows").c_str());
    string line;
    vector<string> results;
    vector<MutationCallResult> mcr_results;
    results.push_back(MutationCallResult::getHeader());
    while (getline(varfile, line)) {
        CandidateWindow cw = CandidateWindow(line);
        MutationCallResult result = caller.call(cw);
        mcr_results.push_back(result);
        results.push_back(result.getOutput());
    }
    for (int i = 0; i < results.size(); i++) {
        LOG(logDEBUG) << results[i] << endl;
    }
    EXPECT_EQ(results.size(), 2);
    EXPECT_GT(mcr_results[0].getBayesFactorBasic(), 1);
}

TEST_F(HapMuCTest, TotalDel) {
    MyBam tumorBam = MyBam(data_path + "/ins/tumor.bam");
    BamReader tumorBamReader = BamReader(&tumorBam, fai);
    MyBam normalBam = MyBam(data_path + "/ins/normal.bam");
    BamReader normalBamReader = BamReader(&normalBam, fai);
    MutationCaller caller = MutationCaller(tumorBamReader, normalBamReader, params, fai);
    ifstream varfile((data_path + "/ins/windows").c_str());
    string line;
    vector<string> results;
    vector<MutationCallResult> mcr_results;
    results.push_back(MutationCallResult::getHeader());
    while (getline(varfile, line)) {
        CandidateWindow cw = CandidateWindow(line);
        MutationCallResult result = caller.call(cw);
        mcr_results.push_back(result);
        results.push_back(result.getOutput());
    }
    for (int i = 0; i < results.size(); i++) {
        LOG(logDEBUG) << results[i] << endl;
    }
    EXPECT_EQ(results.size(), 2);
    EXPECT_GT(mcr_results[0].getBayesFactorBasic(), 1);
}

TEST_F(HapMuCTest, TotalSNPs) {
	MyBam tumorBam = MyBam(data_path + "/snps/tumor.bam");
	BamReader tumorBamReader = BamReader(&tumorBam, fai);
	MyBam normalBam = MyBam(data_path + "/snps/normal.bam");
	BamReader normalBamReader = BamReader(&normalBam, fai);
	MutationCaller caller = MutationCaller(tumorBamReader, normalBamReader, params, fai);
	ifstream varfile((data_path + "/snps/windows").c_str());
	string line;
	vector<string> results;
	vector<MutationCallResult> mcr_results;
	results.push_back(MutationCallResult::getHeader());
	while (getline(varfile, line)) {
		CandidateWindow cw = CandidateWindow(line);
		MutationCallResult result = caller.call(cw);
		mcr_results.push_back(result);
		results.push_back(result.getOutput());
	}
	for (int i = 0; i < results.size(); i++) {
		LOG(logDEBUG) << results[i] << endl;
	}
	EXPECT_EQ(results.size(), 2);
	EXPECT_LT(mcr_results[0].getBayesFactor(), mcr_results[0].getBayesFactorBasic());
	EXPECT_GT(mcr_results[0].getBayesFactorBasic(), 1);
}
