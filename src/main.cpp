/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <iostream>
#include <fstream>

#include "faidx.h"
#include "Parameters.h"
#include "MyBam.hpp"
#include "BamReader.h"
#include "MutationCaller.h"
#include "CandidateWindow.h"
#include "MutationCallResult.h"

int main(int argc, const char *argv[])
{
    //
    // Parse parameters
    //
    Parameters parameters = Parameters(argc, argv);

    //
    // Open log file
    //
    FILE* logFile = fopen((parameters.outFilePrefix + ".log").c_str(), "w");
    Output2FILE::Stream() = logFile;
    LOG(logINFO) << "HapMuC ver1.1" << std::endl;

    //
    // Open reference genome
    //
    faidx_t *fai = NULL;
    fai = fai_load(parameters.refFileName.c_str());
    if (!fai) {
        LOG(logERROR) << "Cannot open reference sequence file." << std::endl;
        exit(1);
    }

    //
    // Prepare bam files
    //
    MyBam tumorBam = MyBam(parameters.tumorBam);
    MyBam normalBam = MyBam(parameters.normalBam);
    BamReader tumorBamReader = BamReader(&tumorBam, fai, parameters.maxReads);
    BamReader normalBamReader = BamReader(&normalBam, fai, parameters.maxReads);

    //
    // Prepare main algorithm
    //
    MutationCaller mutationCaller = MutationCaller(tumorBamReader, normalBamReader,
                                                   parameters, fai);

    //
    // Open input file and output file
    //
    std::ifstream inputWindowStream(parameters.windowFile.c_str());
    std::ofstream outStream((parameters.outFilePrefix + ".calls.txt").c_str());

    //
    // Output column names
    //
    outStream << MutationCallResult::getHeader() << std::endl;

    std::string line;
    while (getline(inputWindowStream, line)) {
        //
        // Parse a candidate window
        //
        CandidateWindow window = CandidateWindow(line);

        LOG(logINFO) << "************************* target: " <<
        window.info.chr << " " << window.info.start <<
        " " << window.target << " *************************" << std::endl;

        try {
            //
            // Evaluate whether the candidate somatic mutation exists or not
            // by calculating Bayes factor
            //
            MutationCallResult result = mutationCaller.call(window);
            
            //
            // Output the result
            //
            outStream << result.getOutput() << std::endl;
            LOG(logDEBUG) << result.getOutput() << std::endl;
            
        } catch (std::string &s) {
            LOG(logERROR) << "something unexpected happened. exit." << std::endl;
            LOG(logERROR) << s << std::endl;
            exit(1);
        }
    }

    outStream.close();
    inputWindowStream.close();
    fai_destroy(fai);

    return 0;
}

