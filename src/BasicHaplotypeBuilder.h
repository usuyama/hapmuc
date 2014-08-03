#ifndef __haplotype_builder__BasicHaplotypeBuilder__
#define __haplotype_builder__BasicHaplotypeBuilder__

#include <string>
#include <vector>

#include "faidx.h"

#include "Variant.h"
#include "Haplotype.hpp"

class BasicHaplotypeBuilder {
   const faidx_t *fai;
public:
    std::string getRefSeq(std::string chr, int lpos, int rpos);
    std::string getHap(const std::string &refseq, int left, int right,
                       const std::vector<Variant> &variants);
    BasicHaplotypeBuilder(const faidx_t *fai);
    Haplotype build(const int targetPos, const std::vector<Variant> &variants, std::string chr, int left, int right);
};

#endif /* defined(__haplotype_builder__BasicHaplotypeBuilder__) */
