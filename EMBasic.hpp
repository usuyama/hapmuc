/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#ifndef EMBasic_HPP_
#define EMBasic_HPP_
#include <iostream>
namespace EMBasic {
    void estimate(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & variantPosteriors, double a0, string program, Parameters params);
}
#endif