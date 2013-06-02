//
//  EMfor2.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/12/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//
#ifndef EM2_HPP_
#define EM2_HPP_
#include <iostream>
namespace EMfor2 {
    void estimate_basic(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors, uint32_t candPos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & variantPosteriors, Parameters params, double a0, double b0);
    
}
#endif