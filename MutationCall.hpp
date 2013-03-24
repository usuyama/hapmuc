//
//  MutationCall.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/10/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//

#ifndef MUTATION_CALL_HPP_
#define MUTATION_CALL_HPP_

#include <iostream>
namespace MutationCall {
    
    void computeBayesFactors(int index, const vector<Read> & normalReads, const vector<Read> & tumorReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants,  OutputData & oData, OutputData & glfData, Parameters params, string refSeq, string refSeqForAlign, vector<AlignedVariant> & close_somatic, vector<AlignedVariant> & close_germline);
    
    void outputLowerBounds(vector<Haplotype> &haps, string fname, uint32_t leftPos, uint32_t rightPos, double bf, lower_bound_t& normal_lb, lower_bound_t &tumor_lb, lower_bound_t &merged_lb, vector<double> normalHapFreqs, vector<double> tumorHapFreqs, vector<double> mergedHapFreqs);
    
    double calc_hap2_bf(const vector<Read> & normalReads, const vector<Read> & tumorReads, const vector<Read> & mergedReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, OutputData & glfData, Parameters params, string refSeq, string refSeqForAlign, double a0, double b0, double free_a0);
}
#endif