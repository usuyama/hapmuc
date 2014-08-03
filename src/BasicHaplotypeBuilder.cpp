#include "BasicHaplotypeBuilder.h"

#include <algorithm>
#include <sstream>

#include "Haplotype.hpp"
#include "log.h"

BasicHaplotypeBuilder::BasicHaplotypeBuilder(const faidx_t *_fai) : fai(_fai) {};

// get referense sequence
// NOTE: 1-based positions
std::string BasicHaplotypeBuilder::getRefSeq(std::string chr, int lpos, int rpos) {
    if (!fai) {
        throw std::string("FAI error.");
    }
    std::stringstream ss;
    ss << chr << ":" << lpos << "-" << rpos;
    char *ref;
    int len;
    ref = fai_fetch(fai, ss.str().c_str(), &len);
    if (len == 0) {
        throw std::string("faidx error: len==0");
    }
    std::string res(ref);
    free(ref);
    std::transform(res.begin(), res.end(), res.begin(), ::toupper);
    return res;
}

std::string BasicHaplotypeBuilder::getHap(const std::string &refSeq, int left, int right, const std::vector<Variant> &variants) {
    // refSeq            :  1 2 3 4 5 6
    // 0-based positions : 0 1 2 3 4 5 6
    std::vector<int> positionOfIndexInHap(refSeq.size() + 1);
    for (int i = 0; i < positionOfIndexInHap.size(); i++) {
        positionOfIndexInHap[i] = left + i;
    }
    std::string seq = refSeq;
    for (int i = 0; i < variants.size(); i++) {
        Variant var = variants[i];
        std::vector<int>::iterator it = find(positionOfIndexInHap.begin(), positionOfIndexInHap.end(), var.startInGenome);
        if (it != positionOfIndexInHap.end()) {
            int idx = (int) distance(positionOfIndexInHap.begin(), it);
            if (var.type == Variant::DEL) {
                std::string delseq = seq.substr(idx, var.length);
                if (delseq != var.originalString.substr(1, var.length)) {
                    LOG(logERROR) << "the deleted sequence differ from the reference" << std::endl;
                    throw std::string("error_making_haplotype");
                }
                seq.erase(idx, var.length);
                positionOfIndexInHap.erase(it, it + var.length);
            } else if (var.type == Variant::INS) {
                seq.insert(idx, var.seq);
                positionOfIndexInHap.insert(it, var.seq.length(), var.startInGenome);
            } else if (var.type == Variant::SNP) {
                char ref = var.seq[0];
                char nuc = var.seq[3];
                if (seq[idx] == ref && seq[idx] != nuc) {
                    seq[idx] = nuc;
                } else {
                    LOG(logWARNING) << "variant matched to the reference. is it OK?" << std::endl;
                }
            }
        } else {
            throw std::string("error_making_haplotype");
        }
    }
    return seq;
}

Haplotype BasicHaplotypeBuilder::build(const int targetPos, const std::vector<Variant> &variants,
                                  std::string chr, int left, int right) {
    // check overlaps between variants
    for (int i = 0; i < variants.size(); i++) {
        for (int j = i + 1;j < variants.size(); j++) {
            if (variants[i].isOverlap(variants[j])) {
                throw std::string("variants_overlapped_during_build");
            }
        }
    }
    
    std::string refSeq = getRefSeq(chr, left + 1, right);
    std::string seq = getHap(refSeq, left, right, variants);
    Haplotype hap;
    hap.chr = chr;
    hap.leftPos = left;
    hap.rightPos = right;
    hap.seq = seq;
    hap.variants = variants;
    hap.targetPos = targetPos;
    hap.hasIndel = false;
    for (int i = 0; i < variants.size(); i++) {
        if (variants[i].isIndel()) {
            hap.hasIndel = true;
        }
    }
    return hap;
}

