/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <string>
#include <vector>
#include <algorithm> 
#include <numeric>

#include "faidx.h"

#include "log.h"
#include "HaplotypeBuilder.h"
#include "Parameters.h"
#include "Haplotype.hpp"
#include "haplotype_phaser/Read.h"
#include "haplotype_phaser/OptimalDPPhaser.h"
#include "HaplotypePhaserUtils.h"
#include "VariantUtils.h"
#include "Utils.h"
#include "ReadsProfile.h"
#include "VariantFromSAM.h"

HaplotypeBuilder::HaplotypeBuilder(const faidx_t *fai_) : BasicHaplotypeBuilder(fai_), fai(fai_) {
    windowSize = 2000;
};

std::vector<Haplotype> HaplotypeBuilder::buildHaplotypes(const Variant &targetVariant,
                                                         const std::vector<Variant> &closeVariants,
                                                         std::string chr,
                                                         const std::vector<IRead *> &tumorReads,
                                                         const std::vector<IRead *> &normalReads) {
    std::vector<Haplotype> haps;

    // setting left & right pos
    int left = targetVariant.startInGenome - windowSize / 2;
    left = left > 0 ? left : 0;
    int right = targetVariant.startInGenome + windowSize / 2;

   // if (closeVariants.empty()) {
        std::vector<Variant> variants;
        variants.push_back(targetVariant);
        Haplotype h01 = BasicHaplotypeBuilder::build(targetVariant.startInGenome, std::vector<Variant>(), chr, left, right);
        Haplotype h23 = BasicHaplotypeBuilder::build(targetVariant.startInGenome, variants, chr, left, right);
        haps.push_back(h01);
        haps.push_back(h01);
        haps.push_back(h23);
        haps.push_back(h23);
//    }

    return haps;
}

std::vector<Haplotype> HaplotypeBuilder::constructNormalHaplotypes(const Variant &targetVariant,
                                                                   const std::vector<Variant> &closeVariants,
                                                                   std::string chr,
                                                                   const std::vector<IRead *> &tumorReads,
                                                                   const std::vector<IRead *> &normalReads,
                                                                   int *_MECScore, int *_MECNumReads) {
    std::vector<Haplotype> haps;
    
    // setting left & right pos
    int pos = targetVariant.startInGenome;
    int left = pos - windowSize / 2;
    left = left > 0 ? left : 0;
    int right = pos + windowSize / 2;
    
    if (closeVariants.empty()) {
        LOG(logINFO) << "no close germline variants" << std::endl;
        Haplotype h01 = BasicHaplotypeBuilder::build(pos, std::vector<Variant>(), chr, left, right);
        haps.push_back(h01);
        haps.push_back(h01);
    } else if (closeVariants.size() == 1) {
        LOG(logINFO) << "only one close germline variant" << std::endl;
        Haplotype h0 = BasicHaplotypeBuilder::build(pos, std::vector<Variant>(), chr, left, right);
        Haplotype h1 = BasicHaplotypeBuilder::build(pos, closeVariants, chr, left, right);
        haps.push_back(h0);
        haps.push_back(h1);
    } else {
        LOG(logINFO) << "infer local normal haplotypes by a haplotype-phasing algorithm" << std::endl;
        const int maxCloseVariants = 20;
        if (closeVariants.size() > maxCloseVariants) {
            LOG(logWARNING) << "too many variants for inference of local normal haplotypes by haplotype-phasing algorithm" << std::endl;
            std::vector<Variant> newCloseVariants;
            for (int i = 0;i < maxCloseVariants; i++) {
                newCloseVariants.push_back(closeVariants[i]);
            }
            return constructNormalHaplotypes(targetVariant, newCloseVariants, chr, tumorReads, normalReads, _MECScore, _MECNumReads);
        }
        // prepare inputs
        std::vector<VariantFromSAM> variants = VariantUtils::variantVecToVariantFromSAMVec(closeVariants);
        std::vector<haplotype_phaser::Read> hpReads;
    
        while (true) {
            std::vector<haplotype_phaser::Read> hpReads_normal = HaplotypePhaserUtils::generateIReadVecFromReadVec(variants, normalReads);
            std::vector<haplotype_phaser::Read> hpReads_tumor = HaplotypePhaserUtils::generateIReadVecFromReadVec(variants, tumorReads);
            hpReads.insert(hpReads.end(), hpReads_normal.begin(), hpReads_normal.end());
            hpReads.insert(hpReads.end(), hpReads_tumor.begin(), hpReads_tumor.end());

            // check
             std::vector<bool> atLeastOne(variants.size(), false);
            for (int i = 0; i < hpReads.size(); i++) {
                std::string &seq = hpReads[i].seq;
                for (int si = 0; si < seq.length(); si++) {
                    if (seq[si] != '-') {
                        atLeastOne[si] = true;
                    }
                }
            }

            std::vector<VariantFromSAM> newVariants;
            for (int x = 0; x < variants.size(); x++) {
                if (atLeastOne[x]) {
                    newVariants.push_back(variants[x]);
                } else {
                    LOG(logWARNING) << variants[x] << " was not covered" << std::endl;
                }
            }

            if (newVariants.size() != variants.size()) {
                hpReads.clear();
                variants = newVariants;
            } else {
                break;
            }
        }

        // check if there are enough haplotype-informative reads
        if (!hpReads.empty()) {
            std::vector<bool> bridged_intervals(variants.size() - 1, false);
            for (int i = 0; i < hpReads.size(); i++) {
                std::string &seq = hpReads[i].seq;
                // LOG(logDEBUG) << seq << std::endl;
                int first_index = (int)variants.size(), last_index = -1;
                for (int si = 0; si < seq.length(); si++) {
                    if (seq[si] != '-') {
                        if (first_index == (int)variants.size()) {
                            first_index = si;
                        }
                        last_index = si;
                    }
                }
                //LOG(logDEBUG) << first_index << "~" << last_index << std::endl;
                for (int x = first_index; x < last_index; x++) {
                    bridged_intervals[x] = true;
                }
            }

            LOG(logDEBUG) << "check bridged intervals" << std::endl;
            for (int i = 0; i < bridged_intervals.size(); i++) {
                // LOGP(logDEBUG) << bridged_intervals[i] << " ";
                if (bridged_intervals[i] == false) {
                    LOG(logWARNING) << "There is a gap between close variants for haplotype-phasing" << std::endl;
                    std::vector<bool> bridged_intervals(variants.size() - 1, false);
                    for (int i = 0; i < hpReads.size(); i++) {
                        std::string &seq = hpReads[i].seq;
                        LOG(logDEBUG) << seq << std::endl;
                    }
                    hpReads.clear();
                    break;
                }
            }
        }
        
        LOG(logDEBUG) << "read vec for haplotype-phasing prepared" << std::endl;
        if (hpReads.empty()) {
            LOG(logWARNING) << "no haplotype-informative reads for constructing"
            << "normal haplotypes: " << hpReads.size() << std::endl;
            const Variant &closestVariant = VariantUtils::findClosestVariant(targetVariant, closeVariants);
            LOG(logWARNING) << "only use the nearlest variant: " << closestVariant << std::endl;
            std::vector<Variant> newCloseVariants;
            newCloseVariants.push_back(closestVariant);
            return constructNormalHaplotypes(targetVariant, newCloseVariants, chr, tumorReads, normalReads, _MECScore, _MECNumReads);
        }
        
        // do haplotype phasing
        haplotype_phaser::OptimalDPPhaser phaser(hpReads);
        std::string hap;
        int MEC;
        phaser.inferHaplotype(&hap, &MEC);
        LOG(logDEBUG) << phaser.toString(hap, MEC);

        if (MEC / (float)hpReads.size() > 1.0) {
            LOG(logERROR) << "very low-quality haplotype-phasing " << MEC << " / " << hpReads.size() << std::endl;
            const Variant &closestVariant = VariantUtils::findClosestVariant(targetVariant, closeVariants);
            LOG(logWARNING) << "only use the nearlest variant: " << closestVariant << std::endl;
            std::vector<Variant> newCloseVariants;
            newCloseVariants.push_back(closestVariant);
            return constructNormalHaplotypes(targetVariant, newCloseVariants, chr, tumorReads, normalReads, _MECScore, _MECNumReads);
        }

        // set variables for the result
        *_MECScore = MEC;
        *_MECNumReads = (int)hpReads.size();

        // process output
        std::vector<std::vector<VariantFromSAM> > vs = HaplotypePhaserUtils::HapToVariants(hap, variants);
        std::vector<Variant> vars0 = VariantUtils::variantFromSAMVecToVariantVec(vs[0]);
        std::vector<Variant> vars1 = VariantUtils::variantFromSAMVecToVariantVec(vs[1]);
        
        //generate haps
        Haplotype h0 = BasicHaplotypeBuilder::build(pos, vars0, chr, left, right);
        Haplotype h1 = BasicHaplotypeBuilder::build(pos, vars1, chr, left, right);
        haps.push_back(h0);
        haps.push_back(h1);
    }
    
    LOG(logDEBUG) << "normal haplotypes:" << std::endl;
    for (std::vector<Haplotype>::iterator it = haps.begin(); it != haps.end(); it++) {
        for (std::vector<Variant>::iterator vit = it->variants.begin(); vit != it->variants.end(); vit++) {
            LOGP(logDEBUG) << *vit << " ";
        }
        LOGP(logDEBUG) <<std::endl;
    }
         
    return haps;
}

std::vector<Haplotype> HaplotypeBuilder::prepareHaplotypesForBasic(std::string chr,const Variant &targetVariant) {
    std::vector<Haplotype> haps;
    int pos = targetVariant.startInGenome;
    int left = pos - windowSize / 2;
    left = left > 0 ? left : 0;
    int right = pos + windowSize / 2;
    Haplotype h0 = BasicHaplotypeBuilder::build(pos, std::vector<Variant>(), chr, left, right);
    haps.push_back(h0);
    std::vector<Variant> variants;
    variants.push_back(targetVariant);
    Haplotype h2 = BasicHaplotypeBuilder::build(pos, variants, chr, left, right);
    haps.push_back(h2);
    return haps;
}

std::vector<Haplotype> HaplotypeBuilder::addSomaticAndErrorHaplotypes(const std::vector<Haplotype> &normalHaps,
                                                                      const Variant &targetVariant,
                                                                      const std::vector<Variant> &closeVariants,
                                                                      std::string chr,
                                                                      const std::vector<IRead *> &tumorReads,
                                                                      const std::vector<IRead *> &normalReads) {
    std::vector<Haplotype> haps = normalHaps;
    
    // setting left & right pos
    int pos = targetVariant.startInGenome;
    int left = pos - windowSize / 2;
    left = left > 0 ? left : 0;
    int right = pos + windowSize / 2;
    
    // generate somatic and error haplotypes from h0 and h1
    for (std::vector<Haplotype>::const_iterator it = normalHaps.begin(); it != normalHaps.end();it++) {
        std::vector<Variant> variants = it->variants;
        variants.push_back(targetVariant);
        Haplotype newHaplotype = BasicHaplotypeBuilder::build(pos, variants, chr, left, right);
        haps.push_back(newHaplotype);
    }
    
    return haps;
}

static inline bool checkMapQuality(ReadsProfile::Profile &profile, double threshold) {
    return profile.numReads == 0 || profile.averageMapQuality >= threshold;
}

static inline bool checkIndelPosession(ReadsProfile::Profile &profile, double threshold) {
    return profile.numReads == 0 || profile.indelPosesssionFreq <= threshold;
}

static inline bool checkSoftClipPosession(ReadsProfile::Profile &profile, double threshold) {
    return profile.numReads == 0 || profile.softClipPosessionFreq <= threshold;
}

std::vector<Variant> HaplotypeBuilder::checkCloseVariants(const Variant &targetVariant,
                                                          const std::vector<Variant> &closeVariants,
                                                          std::string chr,
                                                          const std::vector<IRead *> &tumorReads,
                                                          const std::vector<IRead *> &normalReads,
                                                          Parameters &params, MutationCallResult &result) {
    //
    // check reads that cover the target position
    //
    ReadsProfile tumorReadsProfile(tumorReads, targetVariant);
    ReadsProfile normalReadsProfile(normalReads, targetVariant);
    LOG(logINFO) << "profiles of reads, which cover the candidate position:" << std::endl;
    LOGP(logINFO) << "tumor reads: " << tumorReadsProfile << std::endl;
    LOGP(logINFO) << "normal reads: " << normalReadsProfile << std::endl;

    std::string errorStatus = "";

    {
        result.avgMapQualTumorObs = tumorReadsProfile.withVariant.averageMapQuality;
        result.avgMapQualTumorRef = tumorReadsProfile.withoutVariant.averageMapQuality;
        result.avgMapQualNormalObs = normalReadsProfile.withVariant.averageMapQuality;
        result.avgMapQualNormalRef = normalReadsProfile.withoutVariant.averageMapQuality;

        double &th = params.averageMapQualThreshold;
        if (!(checkMapQuality(tumorReadsProfile.withVariant, th) &&
              checkMapQuality(tumorReadsProfile.withoutVariant, th) &&
              checkMapQuality(normalReadsProfile.withVariant, th) &&
              checkMapQuality(normalReadsProfile.withoutVariant, th))) {
            LOG(logWARNING) << "the average mapping quality is lower than the specified threshold (";
            LOGP(logWARNING) << th << ")" << std::endl;
            LOG(logWARNING) << "give up considering close germline variants." << std::endl;
            errorStatus = "low_mapping_quality";
        }
    }

    {
        result.softClipFreqTumorObs = tumorReadsProfile.withVariant.softClipPosessionFreq;
        result.softClipFreqTumorRef = tumorReadsProfile.withoutVariant.softClipPosessionFreq;
        result.softClipFreqNormalObs = normalReadsProfile.withVariant.softClipPosessionFreq;
        result.softClipFreqNormalRef = normalReadsProfile.withoutVariant.softClipPosessionFreq;

        double &th = params.softClipPosessionFreqThreshold;
        if (errorStatus == "" &&
            !(checkSoftClipPosession(tumorReadsProfile.withVariant, th) &&
              checkSoftClipPosession(tumorReadsProfile.withoutVariant, th) &&
              checkSoftClipPosession(normalReadsProfile.withVariant, th) &&
              checkSoftClipPosession(normalReadsProfile.withoutVariant, th))) {
            LOG(logWARNING) << "SoftClip posession rate is higher than the threshold (";
            LOGP(logWARNING) << th << ")" << std::endl;
            LOG(logWARNING) << "give up considering close germline variants." << std::endl;
            errorStatus = "too_many_softclips_nearby";
        }
    }

    {
        result.indelFreqTumorObs = tumorReadsProfile.withVariant.indelPosesssionFreq;
        result.indelFreqTumorRef = tumorReadsProfile.withoutVariant.indelPosesssionFreq;
        result.indelFreqNormalObs = normalReadsProfile.withVariant.indelPosesssionFreq;
        result.indelFreqNormalRef = normalReadsProfile.withoutVariant.indelPosesssionFreq;

        double &th = params.indelPosessionFreqThreshold;
        if (errorStatus == "" &&
            !(checkIndelPosession(tumorReadsProfile.withVariant, th) &&
              checkIndelPosession(tumorReadsProfile.withoutVariant, th) &&
              checkIndelPosession(normalReadsProfile.withVariant, th) &&
              checkIndelPosession(normalReadsProfile.withoutVariant, th))) {
            LOG(logWARNING) << "Indel posession rate is higher than the threshold (";
            LOGP(logWARNING) << th << ")" << std::endl;
            LOG(logWARNING) << "give up considering close germline variants." << std::endl;
            throw std::string("too_many_indels_nearby");
        }
    }

    if (errorStatus != "") {
        throw errorStatus;
    }

    for (std::vector<Variant>::const_iterator it = closeVariants.begin();it != closeVariants.end();it++) {
        if (it->isIndel() && targetVariant.isOverlap(*it, params.minDistanceGermlineIndel)) {
            LOG(logWARNING) << "too close germline indel to the tarfget" << std::endl;
            throw std::string("germline_indel_too_close");
        }
    }

    for (std::vector<Variant>::const_iterator it = closeVariants.begin();it != closeVariants.end();it++) {
        if (targetVariant.isOverlap(*it, params.minDistanceSNP)) {
            LOG(logWARNING) << "too close germline variants to the target" << std::endl;
            throw std::string("germline_variant_too_close");
        }
    }

    for (int i = 0; i < closeVariants.size(); i++) {
        for (int j = i + 1;j < closeVariants.size(); j++) {
            if (closeVariants[i].isOverlap(closeVariants[j], params.minDistanceSNP)) {
                LOG(logWARNING) << "too close germline variants each other" << std::endl;
                throw std::string("germline_variants_overlapped");
            }
        }
    }

    if (closeVariants.empty()) {
        throw std::string("no_close_germline_variants");
    }
    
    return closeVariants;
}

static inline std::string votesToString(std::vector<double> &votes) {
    std::stringstream ss;
    for (int i = 0; i < votes.size(); i++) {
        ss << votes[i] << " ";
    }
    return ss.str();
}

void HaplotypeBuilder::checkAndSort(std::vector<Haplotype> *haps,
                                    std::vector<std::vector<double> > *normalLiks,
                                    std::vector<std::vector<double> > *tumorLiks,
                                    const std::vector<IRead *> &normalReads,
                                    const std::vector<IRead *> &tumorReads,
                                    const Variant &targetVariant,
                                    int *numHaplotypeInfomativeReads) {
    if (haps->size() != 4) {
        throw std::string("checkAndSort is for haps.size() == 4");
    }

    if (normalLiks->size() < 4 || tumorLiks->size() < 4) {
        throw std::string("# of reads/liks too small");
    }

    // TODO: filter low-quality reads
    std::vector<double> flatNormalLiks = std::vector<double>();
    for (int i = 0;i < normalLiks->size();i++) flatNormalLiks.insert(flatNormalLiks.end(), normalLiks->at(i).begin(), normalLiks->at(i).end());
    double maxNormalLiks = *std::max_element(flatNormalLiks.begin(), flatNormalLiks.end());
    LOG(logDEBUG) << "maxNormalLiks: " << maxNormalLiks << std::endl;
    for (int i = 0;i < tumorLiks->size();i++) {
        bool flag = false;
        double max = *std::max_element(tumorLiks->at(i).begin(), tumorLiks->at(i).end());
        if ((max - maxNormalLiks) < -25.0) flag = true;
        if (flag) {
            LOG(logDEBUG) << "filter low-quality reads" << std::endl;
            for (int j = 0;j < 4;j++) {
                LOGP(logDEBUG) << tumorLiks->at(i)[j] << " ";
                tumorLiks->at(i)[j] = max;
            }
            LOGP(logDEBUG) << std::endl;
        }
    }
    for (int i = 0;i < normalLiks->size();i++) {
        bool flag = false;
        double max = *std::max_element(normalLiks->at(i).begin(), normalLiks->at(i).end());
        if ((max - maxNormalLiks) < -25.0) flag = true;
        if (flag) {
            LOG(logDEBUG) << "filter low-quality reads" << std::endl;
            for (int j = 0;j < 4;j++) {
                LOGP(logDEBUG) << normalLiks->at(i)[j] << " ";
                normalLiks->at(i)[j] = max;
            }
            LOGP(logDEBUG) << std::endl;
        }
    }

    VariantFromSAM target = VariantUtils::variantToVariantFromSAM(targetVariant);

    // check the number of haplotype-informative-variant-supporting reads in tumor
    int count = 0;
    for (int i = 0;i < tumorLiks->size();i++) {
        std::vector<double> tmp = tumorLiks->at(i);
        double max = tmp[0];
        int max_index = 0;
        for (int j = 1;j < tmp.size();j++) {
            if (max < tmp[j]) {
                max = tmp[j];
                max_index = j;
            }
        }
        if (std::abs(tmp[2] - tmp[3]) > 10e-2 && (max_index == 2 || max_index == 3)) {
            count++;
        }
    }

    *numHaplotypeInfomativeReads = count;

    LOG(logINFO) << "# of haplotype-informative variant-supporting reads: " << count << std::endl;
    if (count < 2) {
        throw std::string("too_few_haplotype_informative_variant_supporting_reads");
    }

    std::vector<double> normalVotes = Utils::votes(*normalLiks), tumorVotes = Utils::votes(*tumorLiks);
    int normalReadsSize = (int)normalLiks->size();
    int tumorReadsSize = (int)tumorLiks->size();

    if ((normalVotes[0] + normalVotes[1] < normalVotes[2] + normalVotes[3]) ||
        (normalVotes[0] < normalReadsSize / 4.0 && normalVotes[1] < normalReadsSize / 4.0)) {
        LOG(logERROR) << "normalVotes ";
        LOGP(logERROR) << votesToString(normalVotes) << std::endl;
        LOG(logERROR) << "tumorVotes ";
        LOGP(logERROR) << votesToString(tumorVotes) << std::endl;
        throw std::string("umbalanced_normal_haplotypes");
    }

    if (tumorVotes[2] < tumorVotes[3]) {
        LOG(logDEBUG) << "switch h0 <-> h1 and h2 <-> h3" << std::endl;
        std::iter_swap(haps->begin(), haps->begin() + 1);
        std::iter_swap(haps->begin() + 2, haps->begin() + 3);
        for (int i = 0; i < normalLiks->size(); i++) {
            std::vector<double> &tmp = normalLiks->at(i);
            std::swap(tmp[0], tmp[1]);
            std::swap(tmp[2], tmp[3]);
        }
        for (int i = 0; i < tumorLiks->size(); i++) {
            std::vector<double> &tmp = tumorLiks->at(i);
            std::swap(tmp[0], tmp[1]);
            std::swap(tmp[2], tmp[3]);
        }
    }

    ReadsProfile tumorReadsProfile(tumorReads, targetVariant);

    for (int i = 0;i < tumorLiks->size();i++) {
        IRead *read = tumorReads.at(i);
        if (*numHaplotypeInfomativeReads >= 3
            && read->hasSoftClip() == false
            && read->hasIndel() == false
            && read->getMapQuality() > 0.9999
            && checkIndelPosession(tumorReadsProfile.withVariant, 0.2)
            && checkIndelPosession(tumorReadsProfile.withoutVariant, 0.2)
            && checkSoftClipPosession(tumorReadsProfile.withVariant, 0.2)
            && checkSoftClipPosession(tumorReadsProfile.withoutVariant, 0.2)
            && checkMapQuality(tumorReadsProfile.withVariant, 0.99)
            && checkMapQuality(tumorReadsProfile.withoutVariant, 0.99)) {
            continue;
        }

        std::vector<double> &tmp = tumorLiks->at(i);
        double max = tmp[0];
        int max_index = 0;
        for (int j = 1;j < tmp.size();j++) {
            if (max < tmp[j]) {
                max = tmp[j];
                max_index = j;
            }
        }
        if (max_index == 2) {
            LOG(logDEBUG) << "ignoring " << read->getSeqName() << std::endl;
            tmp[3] = tmp[2];
            tmp[1] = tmp[0];
        }
    }

    LOG(logDEBUG) << "tumor liks: " << std::endl;
    for (int i = 0;i < tumorLiks->size();i++) {
        for (int j = 0;j < tumorLiks->at(i).size();j++) {
            LOGP(logDEBUG) << tumorLiks->at(i)[j] << " ";
        }
        LOGP(logDEBUG) << tumorReads[i]->getSeqName() << std::endl;
    }

    LOG(logDEBUG) << "normal liks: " << std::endl;
    for (int i = 0;i < normalLiks->size();i++) {
        for (int j = 0;j < normalLiks->at(i).size();j++) {
            LOGP(logDEBUG) << normalLiks->at(i)[j] << " ";
        }
        LOGP(logDEBUG) << normalReads[i]->getSeqName() << std::endl;
    }

}


// extract h0 and h3
std::vector<std::vector<double> > HaplotypeBuilder::convertLiks2basic(const std::vector<std::vector<double> > &liks) {
    if (liks.size() == 0) {
        throw std::string("liks.size() == 0");
    }
    if (liks[0].size() == 4) {
        std::vector<std::vector<double> > newLiks;
        newLiks.reserve(liks.size());
        for (int i = 0; i < liks.size(); i++) {
            std::vector<double> tmp;
            tmp.push_back(liks[i][0]);
            tmp.push_back(liks[i][0]);
            tmp.push_back(liks[i][2]);
            tmp.push_back(liks[i][2]);
            newLiks.push_back(tmp);
        }
        return newLiks;
    } else if (liks[0].size() == 2) {
        std::vector<std::vector<double> > newLiks;
        newLiks.reserve(liks.size());
        for (int i = 0; i < liks.size(); i++) {
            std::vector<double> tmp;
            tmp.push_back(liks[i][0]);
            tmp.push_back(liks[i][0]);
            tmp.push_back(liks[i][1]);
            tmp.push_back(liks[i][1]);
            newLiks.push_back(tmp);
        }
        return newLiks;
    } else {
        throw std::string("liks[0].size() != 2 or 4");
    }
}

std::vector<Haplotype> HaplotypeBuilder::convertHaps2basic(const std::vector<Haplotype> &haps) {
    if (haps.size() == 4) {
        std::vector<Haplotype> newHaps;
        newHaps.push_back(haps[0]);
        newHaps.push_back(haps[0]);
        newHaps.push_back(haps[2]);
        newHaps.push_back(haps[2]);
        return newHaps;
    } else if (haps.size() == 2) {
        std::vector<Haplotype> newHaps;
        newHaps.push_back(haps[0]);
        newHaps.push_back(haps[0]);
        newHaps.push_back(haps[1]);
        newHaps.push_back(haps[1]);
        return newHaps;
    } else {
        throw std::string("haps.size() != 4 or 2");
    }
}
