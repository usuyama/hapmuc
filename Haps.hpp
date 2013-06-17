//
//  MutationCall.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/10/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//
#ifndef HAPS_HPP_
#define HAPS_HPP_

#include <vector>
#include <map>
#include "Variant.hpp"
#include "HapMuC.hpp" //for params
#include "Read.hpp"
#include "Haplotype.hpp"
#include "MLAlignment.hpp"

namespace Haps {
    
    bool alignHaplotypes(vector<Haplotype> & haps,  uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants, Parameters params, string refSeqForAlign);
    void selectHaplotypesAndReads(vector<Haplotype> & haps, vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, int & normal_count, int & tumor_count, Parameters params);
    void filterHaplotypes(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks,  vector<int> & filtered, map<pair<int, AlignedVariant>, VariantCoverage> & varCoverage, bool doFilter, Parameters params);
    bool getHaplotypes(vector<Haplotype> &haps, const vector<Read> & reads,uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign);
    void computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap, Parameters params);
    void filterReads(vector<Read> & mergedReads, vector<vector<MLAlignment> > & mergedLiks, int & normal_count, int & tumor_count);
    
}
#endif
