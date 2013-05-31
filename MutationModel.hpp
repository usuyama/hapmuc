//
//  EMfor2.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/12/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//
#ifndef MM_HPP_
#define MM_HPP_
#include "Read.hpp"
#include "HapMuC.hpp"
#include "VariantFile.hpp"
#include "OutputData.hpp"
namespace MutationModel {
    void estimate(const vector<Haplotype> & haps, const vector<Read> & tumor_reads, const vector<Read> & normal_reads, const vector<double> & rlt, const vector<double> & rln, vector<double> & tumorHapFreqs, vector<double> & normalHapFreqs, vector <HapEstResult > & tumorPosteriors, vector <HapEstResult > & normalPosteriors, uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, int index, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & tumorVariantPosteriors, map<AlignedVariant, double> & normalVariantPosteriors, Parameters params, string est_time);
}
#endif
