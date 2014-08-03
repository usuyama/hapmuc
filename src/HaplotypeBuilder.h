/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#ifndef __HaplotypeBuilder__
#define __HaplotypeBuilder__

#include <vector>
#include <string>

#include "faidx.h"

#include "Haplotype.hpp"
#include "Variant.h"
#include "BasicHaplotypeBuilder.h"
#include "IRead.h"
#include "Parameters.h"
#include "MutationCallResult.h"

class HaplotypeBuilder : public BasicHaplotypeBuilder {
    const faidx_t *fai;
public:
    int windowSize;
    std::vector<Haplotype> buildHaplotypes(const Variant &targetVariant, const std::vector<Variant> &closeVariants,
                                           std::string chr,
                                           const std::vector<IRead *> &tumorReads,
                                           const std::vector<IRead *> &normalReads);
    std::vector<Haplotype> constructNormalHaplotypes(const Variant &targetVariant,
                                                     const std::vector<Variant> &closeVariants,
                                                     std::string chr,
                                                     const std::vector<IRead *> &tumorReads,
                                                     const std::vector<IRead *> &normalReads,
                                                     int *_MECScore, int *_MECNumReads);
    std::vector<Haplotype> addSomaticAndErrorHaplotypes(const std::vector<Haplotype> &normalHaps,
                                                        const Variant &targetVariant,
                                                        const std::vector<Variant> &closeVariants,
                                                        std::string chr,
                                                        const std::vector<IRead *> &tumorReads,
                                                        const std::vector<IRead *> &normalReads);
    static std::vector<Variant> checkCloseVariants(const Variant &targetVariant,
                                            const std::vector<Variant> &closeVariants,
                                            std::string chr,
                                            const std::vector<IRead *> &tumorReads,
                                            const std::vector<IRead *> &normalReads,
                                            Parameters &params, MutationCallResult &result);
    HaplotypeBuilder(const faidx_t *fai);
    static void checkAndSort(std::vector<Haplotype> *haps,
                             std::vector<std::vector<double> > *normalLiks,
                             std::vector<std::vector<double> > *tumorLiks,
                             const std::vector<IRead *> &normalReads,
                             const std::vector<IRead *> &tumorReads,
                             const Variant &targetVariant,
                             int *numHaplotypeInfomativeReads);
    static std::vector<std::vector<double> > convertLiks2basic(const std::vector<std::vector<double> > &liks);
    static std::vector<Haplotype> convertHaps2basic(const std::vector<Haplotype> &haps);
    std::vector<Haplotype> prepareHaplotypesForBasic(std::string chr,const Variant &targetVariant);
};

#endif /* defined(__HaplotypeBuilder__) */
