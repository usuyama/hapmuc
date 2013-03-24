//
//  EMBasic.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/12/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//
#ifndef EMBasic_HPP_
#define EMBasic_HPP_
#include <iostream>
namespace EMBasic {
    void estimate(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, OutputData & glfData, int index, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & variantPosteriors, double a0, string program, Parameters params);
}
#endif