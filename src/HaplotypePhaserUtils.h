#ifndef haplotype_builder_HaplotypePhaserUtils_h
#define haplotype_builder_HaplotypePhaserUtils_h

#include <vector>

#include "haplotype_phaser/Read.h"
#include "IRead.h"
#include "VariantFromSAM.h"

class HaplotypePhaserUtils {
public:
    static inline haplotype_phaser::Read IReadToRead(const std::vector<VariantFromSAM> &variants, const IRead *iread) {
        std::string line = "";
        for (std::vector<VariantFromSAM>::const_iterator it = variants.begin(); it != variants.end(); it++) {
            if (iread->coverPosition(it->left) && iread->coverPosition(it->right)) {
                if (iread->hasVariant(*it)) {
                    line += "1";
                } else {
                    line += "0";
                }
            } else {
                line += "-";
            }
        }
        return haplotype_phaser::Read::fromString(line, iread->getSeqName());
    }
    
    static inline std::vector<haplotype_phaser::Read> generateIReadVecFromReadVec(const std::vector<VariantFromSAM> &variants,
                                                                        const std::vector<IRead *> &ireads) {
        std::vector<haplotype_phaser::Read> hpReads;
        for (std::vector<IRead *>::const_iterator it = ireads.begin(); it != ireads.end(); it++) {
            haplotype_phaser::Read r = IReadToRead(variants, *it);
            if (r.isHaplotypeInformative()) {
                hpReads.push_back(r);
            }
        }
        return hpReads;
    }
    
    static inline std::vector<std::vector<VariantFromSAM> > HapToVariants(const std::string &hap,
                                                                          const std::vector<VariantFromSAM> &variants) {
        std::vector<VariantFromSAM> variants1, variants2;
        for (int i = 0;i < hap.length();i++) {
            if (hap[i] == '1') {
                variants1.push_back(variants[i]);
            } else {
                variants2.push_back(variants[i]);
            }
        }
        std::vector<std::vector<VariantFromSAM> > results;
        results.push_back(variants1);
        results.push_back(variants2);
        return results;
    }
    
};

#endif
