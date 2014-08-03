#ifndef __haplotype_builder__VariantUtils__
#define __haplotype_builder__VariantUtils__

#include <vector>
#include <sstream>
#include <algorithm>

#include "Variant.h"
#include "VariantFromSAM.h"
#include "log.h"

class VariantUtils {
public:
    static inline Variant variantFromSAMToVariant(const VariantFromSAM &_v) {
        Variant v;
        v.startInGenome = _v.left;
        v.endInGenome = _v.right;
        v.length = _v.length;
        v.ref = _v.ref;
        v.obs = _v.obs;
        std::stringstream ss;
        if (_v.type == VariantFromSAM::MISMATCH) {
            ss << _v.ref << "=>" << _v.obs;
            v.type = Variant::SNP;
        } else if (_v.type == VariantFromSAM::INS) {
            ss << "+" << _v.obs;
            v.type = Variant::INS;
        } else if (_v.type == VariantFromSAM::DEL) {
            ss << "-" << _v.ref;
            v.type = Variant::DEL;
        } else {
            throw std::string("fromVariantFromSAM_not_supported");
        }
        v.seq = ss.str();
        v.originalString = v.seq;
        return v;
    }
    
    static inline VariantFromSAM variantToVariantFromSAM(const Variant &_v) {
        VariantFromSAM v;
        v.left = _v.startInGenome;
        v.right = _v.endInGenome;
        v.length = _v.length;
        v.ref = _v.ref;
        v.obs = _v.obs;
        if (_v.type == Variant::SNP) {
            v.type = VariantFromSAM::MISMATCH;
        } else if (_v.type == Variant::INS) {
            v.type = VariantFromSAM::INS;
        } else if (_v.type == Variant::DEL) {
            v.type = VariantFromSAM::DEL;
        }
        return v;
    }
    
    static inline std::vector<Variant> variantFromSAMVecToVariantVec(const std::vector<VariantFromSAM> &variants) {
        std::vector<Variant> newVariants;
        std::transform(variants.begin(), variants.end(), back_inserter(newVariants), variantFromSAMToVariant);
        return newVariants;
    }
    
    static inline std::vector<VariantFromSAM> variantVecToVariantFromSAMVec(const std::vector<Variant> &variants) {
        std::vector<VariantFromSAM> newVariants;
        std::transform(variants.begin(), variants.end(), back_inserter(newVariants), variantToVariantFromSAM);
        return newVariants;
    }
    
    static inline int getMinDistance(Variant &target, std::vector<Variant> &closeVariants) {
        int min = 10e5;
        for (std::vector<Variant>::iterator it = closeVariants.begin(); it != closeVariants.end(); it++) {
            int distance = std::abs(target.startInGenome - it->startInGenome);
            if (distance < min) {
                min = distance;
            }
        }
        return min;
    }
    
    static inline const Variant &findClosestVariant(const Variant &target, const std::vector<Variant> &closeVariants) {
        int min = 10e5;
        const Variant *closest;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin(); it != closeVariants.end(); it++) {
            int distance = std::abs(target.startInGenome - it->startInGenome);
            if (distance < min) {
                min = distance;
                closest = &(*it);
            }
        }
        return *closest;
    }
    
    static inline const Variant &findFurthestVariant(const Variant &target, const std::vector<Variant> &closeVariants) {
        int max = 0;
        const Variant *furthest;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin(); it != closeVariants.end(); it++) {
            int distance = std::abs(target.startInGenome - it->startInGenome);
            if (distance > max) {
                max = distance;
                furthest = &(*it);
            }
        }
        return *furthest;
    }
    
    static std::string getSymbols(std::vector<Variant> &variants) {
        if (variants.empty()) {
            return "-";
        } else if (variants.size() == 1) {
            return variants[0].getSymbol();
        }
        std::string out = variants[0].getSymbol();
        for (int i = 1;i < variants.size();i++) {
            out += "," + variants[i].getSymbol();
        }
        return out;
    }
    
    static inline std::pair<int, int> getRange(const Variant &target, const std::vector<Variant> &closeVariants) {
        int left = target.startInGenome, right = target.endInGenome;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin();it != closeVariants.end();it++) {
            if (left > it->startInGenome) {
                left = it->startInGenome;
            }
            if (right < it->endInGenome) {
                right = it->endInGenome;
            }
        }
        return std::make_pair(left, right);
    }
    
    static inline std::vector<Variant> excludeFurthestVariant(const Variant &target, const std::vector<Variant> &closeVariants) {
        if (closeVariants.size() <= 1) {
            return std::vector<Variant>();
        }
        const Variant &furthest = findFurthestVariant(target, closeVariants);
        std::vector<Variant> variants;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin();it != closeVariants.end();it++) {
            if (furthest.startInGenome != it->startInGenome) {
                variants.push_back(*it);
            } else {
                LOG(logWARNING) << "exclude: " << *it << std::endl;
            }
        }
        return variants;
    }
    
};

#endif /* defined(__haplotype_builder__VariantUtils__) */
