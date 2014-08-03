#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_

#include <vector>
#include <string>
#include "Variant.h"

class Haplotype {
public:
    int leftPos, rightPos; // 0-based positions
    std::string seq;
    std::vector<Variant> variants;
    int targetPos; // index for target variant (in genome position)
    bool hasIndel;
    std::string chr;
    
    Haplotype() {};

    Haplotype(const Haplotype &h) {
        chr = h.chr;
        seq = h.seq;
        leftPos = h.leftPos;
        rightPos = h.rightPos;
        hasIndel = h.hasIndel;
        variants = h.variants;
        targetPos = h.targetPos;
    }

    Haplotype &operator=(const Haplotype &h) {
        if (&h != this) {
            chr = h.chr;
            seq = h.seq;
            leftPos = h.leftPos;
            rightPos = h.rightPos;
            hasIndel = h.hasIndel;
            variants = h.variants;
            targetPos = h.targetPos;
        }
        return *this;
    }

    static std::vector<Variant> getAllVariants(const std::vector<Haplotype> &haps) {
        std::vector<Variant> allVariants;
        for (int i = 0; i < haps.size(); i++) {
            allVariants.insert(allVariants.end(), haps[i].variants.begin(), haps[i].variants.end());
        }
        return allVariants;
    }

    // 0-based positions
    int getStartPosInGenome() const {
        return leftPos;
    }

    int getEndPosInGenome() const {
        return rightPos;
    }

    friend std::ostream& operator<<(std::ostream &stream, const Haplotype &hap) {
        stream << "Haplotype " << hap.chr << ":" << hap.leftPos << "-" << hap.rightPos << " " << hap.targetPos << std::endl;
        stream << hap.seq << std::endl << "variants: ";
        for (int i = 0;i < hap.variants.size();i++) {
            stream << hap.variants[i] << " ";
        }
        return stream;
    }
};

#endif /*HAPLOTYPE_HPP_*/

