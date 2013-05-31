//
//  MutationCall.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/10/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//
#ifndef HAPS_HPP_
#define HAPS_HPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include "foreach.hpp"
#include "bam.h"
#include "Haplotype.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
#include "Utils.hpp"
#include "faidx.h"
#include "GetCandidates.hpp"
#include "ObservationModelSeqAn.hpp"
#include "VariantFile.hpp"
#include "Faster.hpp"
#include <ext/hash_map>
#include <exception>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/algorithm/string.hpp>

#include <iostream>
namespace Haps {
    
    bool alignHaplotypes(vector<Haplotype> & haps,  uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants, Parameters params, string refSeqForAlign);
    void selectHaplotypesAndReads(vector<Haplotype> & haps, vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, OutputData & glfData, int & normal_count, int & tumor_count, Parameters params);
    void filterHaplotypes(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks,  vector<int> & filtered, map<pair<int, AlignedVariant>, VariantCoverage> & varCoverage, bool doFilter, Parameters params);
    bool getHaplotypes(vector<Haplotype> &haps, const vector<Read> & reads,uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign);
    void computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap, Parameters params);
    void filterReads(vector<Read> & mergedReads, vector<vector<MLAlignment> > & mergedLiks, int & normal_count, int & tumor_count);
    
    
}
#endif