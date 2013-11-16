/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#ifndef MM_HPP_
#define MM_HPP_
#include "Read.hpp"
#include "HapMuC.hpp"
#include "VariantFile.hpp"
#include "OutputData.hpp"
namespace MutationModel {
    void estimate(const vector<Haplotype> & haps, const vector<Read> & tumor_reads, const vector<Read> & normal_reads, const vector<double> & rlt, const vector<double> & rln, vector<double> & tumorHapFreqs, vector<double> & normalHapFreqs, vector <HapEstResult > & tumorPosteriors, vector <HapEstResult > & normalPosteriors, uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, Hap3Param params, string est_time, string log_prefix);
}
#endif
