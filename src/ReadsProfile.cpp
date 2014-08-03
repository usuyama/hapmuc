#include "ReadsProfile.h"

#include <algorithm>

#include "ReadUtils.h"
#include "VariantUtils.h"

ReadsProfile::Profile ReadsProfile::makeData(const std::vector<IRead *> &reads) {
    Profile d;
    d.numReads = (int) reads.size();
    d.averageMapQuality = 0.0;
    int softClipCount = 0, indelCount = 0;
    for (std::vector<IRead *>::const_iterator it = reads.begin(); it != reads.end(); it++) {
        d.averageMapQuality += (*it)->getMapQualityPrimary();
        if ((*it)->hasSoftClipPrimary()) {
            softClipCount++;
        }
        if ((*it)->hasIndelPrimary()) {
            indelCount++;
        }
    }
    d.averageMapQuality /= d.numReads;
    d.indelPosesssionFreq = (double)indelCount / d.numReads;
    d.softClipPosessionFreq = (double)softClipCount / d.numReads;
    return d;
}

ReadsProfile::ReadsProfile(const std::vector<IRead *> &reads, const Variant &variant) {
    all = makeData(reads);
    VariantFromSAM vsam = VariantUtils::variantToVariantFromSAM(variant);
    // with variant
    std::vector<IRead *> readsWithVariant = ReadUtils::hasVariant(reads, vsam);
    withVariant = makeData(readsWithVariant);
    // without variant
    std::vector<IRead *> readsWithoutVariant;
    readsWithoutVariant.reserve(all.numReads - withVariant.numReads);
    for (std::vector<IRead *>::const_iterator it = reads.begin();it != reads.end();it++) {
        if (std::find(readsWithVariant.begin(), readsWithVariant.end(), *it) == readsWithVariant.end()) {
            readsWithoutVariant.push_back(*it);
        }
    }
    withoutVariant = makeData(readsWithoutVariant);
}
