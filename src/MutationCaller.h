/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef __MutationCall__
#define __MutationCall__

#include "faidx.h"
#include "BamReader.h"
#include "Parameters.h"
#include "MutationCallResult.h"
#include "CandidateWindow.h"
#include "HaplotypeBuilder.h"

class MutationCaller {
public:
    HaplotypeBuilder hapBuilder;
    BamReader &tumorBamReader;
    BamReader &normalBamReader;
    Parameters &params;
    const faidx_t *fai;
    MutationCaller(BamReader &tumorBamReader, BamReader &normalBamReader, Parameters &params, const faidx_t *fai);
    MutationCallResult call(const CandidateWindow &candidateWindow);
};

#endif /* defined(__MutationCall__) */
