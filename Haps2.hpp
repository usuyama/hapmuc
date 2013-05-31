//
//  MutationCall.hpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/10/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//
#ifndef HAPS2_HPP_
#define HAPS2_HPP_

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
#include "HapMuC.hpp"
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
namespace Haps2 {
    void filter_reads(vector<Read> & reads, vector<vector<MLAlignment> > & liks);
    bool alignHaplotypes(vector<Haplotype> & haps,  uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants, Parameters params, string refSeq);
    string get_hap(string refseq, int left, int right, vector<AlignedVariant> variants);
    vector<vector<AlignedVariant > > get_comb_vars(vector<vector<AlignedVariant> > current, vector<AlignedVariant> variants) ;
    vector<string> getAdditionalCombSeq(string refSeq, int left, int right, vector<AlignedVariant> variants, vector<vector<AlignedVariant> > current, Parameters params);
    vector<string> getAllCombSeq(string refSeq, int left, int right, vector<AlignedVariant> variants, Parameters params);
    vector<Haplotype> getAllHaps(string refSeq, string refSeqForAlign, uint32_t & leftPos, uint32_t & rightPos, vector<AlignedVariant> variants, Parameters params);
    bool getHaplotypes(vector<Haplotype> &haps, const vector<Read> & normal_reads, const vector<Read> & tumor_reads, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, vector<AlignedVariant> & close_somatic, vector<AlignedVariant>& close_germline, vector<vector<MLAlignment> > & normal_liks, vector<vector<MLAlignment> > & tumor_liks, OutputData & glfData, int index);
    bool getHaplotypesBasic(vector<Haplotype> &haps, const vector<Read> & normal_reads, const vector<Read> & tumor_reads, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, vector<vector<MLAlignment> > & normal_liks, vector<vector<MLAlignment> > & tumor_liks, OutputData & glfData, int index);
}
#endif