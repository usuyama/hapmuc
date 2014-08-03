//
//  Parameters.cpp
//  ReadTest
//
//  Created by 臼山 直人 on 1/13/14.
//  Copyright (c) 2014 Univ. of Tokyo. All rights reserved.
//

#include <cstdlib>
#include "Parameters.h"
#include "cmdline.h"

Parameters::Parameters() {
    setDefaultValues();
}

Parameters::Parameters(int argc, const char *argv[]) {
    setDefaultValues();
    getFromCommandLineArguments(argc, argv);
}

void Parameters::setDefaultValues() {
    maxReads = 20000;
    minReads = 4;
    averageMapQualThreshold = 0.0;
    mapQualThreshold = 0.9683772233983162;
    //baseQualThreshold = 0.9;
    maxReadLength = 100;
    minDistanceGermlineIndel = 25;
    showHapAlignments = false;
//    showReads = false;
    priorIndel = 1.0 / 10000;
    priorSNP = 1.0 / 1000.0;
//    filterReadAux = std::string("");
    quiet = true;
    singleReads = false;
    maxInsertSize = 800;
    minDistanceSNP = 3;
    maxIndelLength = 25;
    priorDelByError = 5e-4;
    hap3_params.clear_all();
    hap2_params.clear_all();
    parseHyperParameters(hap3_params.mut_a0, "1.0,1.0,1.0");
    parseHyperParameters(hap3_params.mut_b0, "0.1,10.0");
    parseHyperParameters(hap3_params.mut_c0, "1.0,1.0");
    parseHyperParameters(hap3_params.err_a0, "1.0,1.0");
    parseHyperParameters(hap3_params.err_b0, "1.0,10.0");
    parseHyperParameters(hap3_params.err_c0, "1.0,1.0");
    hap2_params = hap3_params;
    EMtol = 1e-4;
    indelPosessionFreqThreshold = 1.0;
    softClipPosessionFreqThreshold = 1.0;
    withoutBayesFactor = false;
}

inline void Parameters::parseHyperParameters(std::vector<double> &vec, std::string s) {
    std::vector<std::string> params = Utils::split(s, ",");
    vec.clear();
    for (int i = 0; i < params.size(); i++) {
        vec.push_back(atof(params[i].c_str()));
    }
}

void Parameters::getFromCommandLineArguments(int argc, const char *argv[]) {
    cmdline::parser a;
    a.add<std::string>("ref", 'f', "fasta reference sequence (should be indexed with .fai file)", true, "");
    a.add<std::string>("tumor", 'a', "tumor bam file", true, "");
    a.add<std::string>("normal", 'b', "normal bam file", true, "");
    a.add<std::string>("out", 'o', "file-prefix for output results", true, "");
    a.add<std::string>("windows", 'w', "file with candidate variants to be tested.", true, "");

    a.add("quiet", '\0', "quiet output");

    a.add<double>("priorSNP", '\0', "prior probability of a SNP site", false, 1.0 / 1000.0);
    a.add<double>("priorIndel", '\0', "prior probability of a detected indel not being a sequencing error", false, 1.0 / 10000.0);
    a.add<int>("maxReads", '\0', "maximum number of reads in a window", false, 20000);
    a.add<int>("minReads", '\0', "minimum number of reads", false, 4);
    a.add<double>("averageMapQualThreshold", '\0', "lower limit of average read mapping quality for considering haplotypes", false, 0.0);
    a.add<double>("mapQualThreshold", '\0', "lower limit of read mapping quality", false, 0.9683772233983162);
    //a.add<double>("baseQualThreshold", '\0', "lower limit for base quality", false, 0.9);
    a.add<int>("maxReadLength", '\0', "maximum length of reads", false, 100);

//    a.add<std::string>("filterReadAux", '\0', "match string for exclusion of reads based on auxilary information", false, "");

    a.add<std::string>("mutationModelAlpha", 'h', "hyper parameter for tumor haplotype frequencies in mutation model.", false, "1.0,1.0,1.0");
    a.add<std::string>("mutationModelBeta", 'i', "hyper parameter for error rates in mutation model.", false, "0.1,10.0");
    a.add<std::string>("mutationModelGamma", 'j', "hyper parameter for normal haplotype frequencies in mutation model.", false, "1.0,1.0");

    a.add<std::string>("errorModelAlpha", 'k', "hyper parameter for tumor haplotype frequencies in error model.", false, "1.0,1.0");
    a.add<std::string>("errorModelBeta", 'l', "hyper parameter for error rates in error model.", false, "1.0,10.0");
    a.add<std::string>("errorModelGamma", 'm', "hyper parameter for normal haplotype frequencies in error model.", false, "1.0,1.0");

    a.add<double>("pError", '\0', "probability of a read indel", false, 5e-4);
    a.add<double>("pMut", '\0', "probability of a mutation in the read", false, 1e-5);
    a.add<int>("maxIndelLength", '\0', "maximum length of a _sequencing error_ indel in read", false, 25);

    a.add<int>("maxInsertSize", '\0', "maximum insert size of paired-end reads", false, 800);
    
    a.add<int>("minDistanceSNP", '\0', "minimum distance to the heterozygous SNP", false, 3);
    a.add<int>("minDistanceGermlineIndel", '\0', "minimum distance to the heterozygous indel", false, 25);
//    a.add("showHapAlignments", '\0', "show for each haplotype which reads map to it");
//    a.add("showReads", '\0', "show reads");
    a.add("singleReads", '\0', "Do not use paired-end information");
    a.add<double>("indelPosessionFreqThreshold", '\0', "Threshold for frequency of reads that have indels", false, 1.0);
    a.add<double>("softClipPosessionFreqThreshold", '\0', "Threshold for frequency of reads that have indels", false, 1.0);

    a.add("withoutBayesFactor", '\0', "Do not calculate Bayes factor (debug purpose)");
    
    a.parse_check(argc, argv);

    maxReads = a.get<int>("maxReads");
    minReads = a.get<int>("minReads");
    averageMapQualThreshold = a.get<double>("averageMapQualThreshold");
    mapQualThreshold = a.get<double>("mapQualThreshold");
    //baseQualThreshold = a.get<double>("baseQualThreshold");
    maxReadLength = a.get<int>("maxReadLength");
    priorSNP = a.get<double>("priorSNP");
    priorIndel = a.get<double>("priorIndel");
    refFileName = a.get<std::string>("ref");
    outFilePrefix = a.get<std::string>("out");
    tumorBam = a.get<std::string>("tumor");
    normalBam = a.get<std::string>("normal");
    windowFile = a.get<std::string>("windows");
    pError = a.get<double>("pError");
    pMut = a.get<double>("pMut");
    maxIndelLength = a.get<int>("maxIndelLength");
    minDistanceSNP = a.get<int>("minDistanceSNP");
    minDistanceGermlineIndel = a.get<int>("minDistanceGermlineIndel");
//    showReads = a.exist("showReads") ? true : false;
    quiet = a.exist("quiet") ? true : false;
    maxInsertSize = a.get<int>("maxInsertSize");
    singleReads = a.exist("singleReads") ? true : false;
//    showHapAlignments = a.exist("showHapAlignments") ? true : false;
//    filterReadAux = a.get<std::string>("filterReadAux");
    parseHyperParameters(hap3_params.mut_a0, a.get<std::string>("mutationModelAlpha"));
    parseHyperParameters(hap3_params.mut_b0, a.get<std::string>("mutationModelBeta"));
    parseHyperParameters(hap3_params.mut_c0, a.get<std::string>("mutationModelGamma"));
    parseHyperParameters(hap3_params.err_a0, a.get<std::string>("errorModelAlpha"));
    parseHyperParameters(hap3_params.err_b0, a.get<std::string>("errorModelBeta"));
    parseHyperParameters(hap3_params.err_c0, a.get<std::string>("errorModelGamma"));
    hap2_params = hap3_params;
    softClipPosessionFreqThreshold = a.get<double>("softClipPosessionFreqThreshold");
    indelPosessionFreqThreshold = a.get<double>("indelPosessionFreqThreshold");
    if (mapQualThreshold > 1.0) {
        throw std::string("mapQualThreshold should be 0.0 ~ 1.0. ex: 30 -> 0.999");
    }
    withoutBayesFactor = a.exist("withoutBayesFactor") ? true : false;
}
