/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include "faidx.h"
#include "MutationCaller.h"
#include "BamReader.h"
#include "Parameters.h"
#include "MutationCallResult.h"
#include "CandidateWindow.h"
#include "HaplotypeBuilder.h"
#include "MutationModel.hpp"
#include "ProfileHMMUtils.h"
#include "ReadUtils.h"
#include "ReadsProfile.h"

MutationCaller::MutationCaller(BamReader &_tumorBamReader,
                               BamReader &_normalBamReader,
                               Parameters &_params,
                               const faidx_t *_fai) : params(_params), tumorBamReader(_tumorBamReader),
normalBamReader(_normalBamReader), hapBuilder(_fai), fai(_fai) {
    hapBuilder.windowSize = 2 * (params.maxInsertSize + params.maxReadLength);
}

MutationCallResult MutationCaller::call(const CandidateWindow &window) {
    // prepare result
    MutationCallResult result(window);

    //
    // get target chr/pos from window
    //
    int pos = window.target.startInGenome;
    std::string chr = window.info.chr;
    LOG(logINFO) << chr << ":" << pos << " " << window.target << std::endl;
    LOG(logINFO) << "close variants: ";
    for (std::vector<Variant>::const_iterator it = window.closeVariants.begin(); it != window.closeVariants.end(); it++) {
        LOGP(logINFO) << *it << " ";
    }
    LOGP(logINFO) << std::endl;

    if (window.closeVariants.empty()) {
        result.setStatus("no_close_germline_variants");
    }

    std::vector<IRead *> allTumorReads, allNormalReads;
    try {
        //
        // get reads around the target
        // NOTE: these IRead * must be deleted manually.
        // adjast the window size by maxReads and variants
        //
        std::vector<Variant> closeVariantsForFetchingReads = window.closeVariants;

        if (params.withoutBayesFactor) {
            closeVariantsForFetchingReads.clear();
        }

        while (true) {
            try {
                std::pair<int, int> windowRange = VariantUtils::getRange(window.target,
                                                                         closeVariantsForFetchingReads);

                if (params.singleReads == false) {
                    LOG(logINFO) << "fetching tumor paired-reads..." << std::endl;

                    allTumorReads = tumorBamReader.getPairedReads(chr, pos,
                                                                  windowRange.first,
                                                                  windowRange.second,
                                                                  params.maxInsertSize);
                } else {
                    LOG(logINFO) << "fetching tumor reads without pair information..." << std::endl;
                    allTumorReads = tumorBamReader.getSingleReads(chr, pos,
                                                                  windowRange.first,
                                                                  windowRange.second);
                }
                if (params.singleReads == false) {
                    LOG(logINFO) << "fetching normal paired-reads..." << std::endl;

                    allNormalReads = normalBamReader.getPairedReads(chr, pos,
                                                                    windowRange.first,
                                                                    windowRange.second,
                                                                    params.maxInsertSize);
                } else {
                    LOG(logINFO) << "fetching normal reads without pair information..." << std::endl;
                    allNormalReads = normalBamReader.getSingleReads(chr, pos,
                                                                    windowRange.first,
                                                                    windowRange.second);
                }
                break;
            } catch (std::string &s) {
                if (s == "too_many_reads_in_window") {
                    if (closeVariantsForFetchingReads.empty()) {
                        throw s;
                    } else {
                        LOG(logWARNING) << "excluding furthest germline variants" << std::endl;
                        std::vector<Variant> tmp = VariantUtils::excludeFurthestVariant(window.target, closeVariantsForFetchingReads);
                        std::swap(tmp, closeVariantsForFetchingReads);
                        LOG(logWARNING) << "try fetching again" << std::endl;
                        continue;
                    }
                } else {
                    LOG(logERROR) << "error during fethcing reads: " << s << std::endl;
                    throw s;
                }
            }
        }

        //
        // check reads
        //
        ReadsProfile allTumorReadsProfile(allTumorReads, window.target);
        ReadsProfile allNormalReadsProfile(allNormalReads, window.target);
        LOG(logINFO) << "profiles of all fetched reads:" << std::endl;
        LOGP(logINFO) << "tumor reads: " << allTumorReadsProfile << std::endl;
        LOGP(logINFO) << "normal reads: " << allNormalReadsProfile << std::endl;

        //
        // reads that cover the target position
        //
        std::vector<IRead *> tmpTumorReads = ReadUtils::coverPosition(allTumorReads, pos);
        std::vector<IRead *> tmpNormalReads = ReadUtils::coverPosition(allNormalReads, pos);


        //
        // check close variants before constructing haplotypes
        //
        std::vector<Variant> closeVariants;
        try {
            closeVariants = HaplotypeBuilder::checkCloseVariants(window.target,
                                                                 closeVariantsForFetchingReads,
                                                                 chr, tmpTumorReads,
                                                                 tmpNormalReads, params, result);
        } catch (std::string &s) {
            closeVariants = std::vector<Variant>();
            if (s == "no_close_germline_variants") {
                result.setStatus(s);
            } else if (s == "low_mapping_quality") {
                result.setStatus(s);
            } else if (s == "germline_variant_too_close") {
                result.setStatus(s);
            } else if (s == "germline_variants_overlapped") {
                result.setStatus(s);
            } else if (s == "germline_indel_too_close") {
                result.setStatus(s);
            } else if (s == "too_many_softclips_nearby") {
                result.setStatus(s);
            } else if (s == "too_many_indels_nearby") {
                result.setStatus(s);
            }
        }

        if (params.withoutBayesFactor == true) {
            // do nothing
        } else {
            //
            // construct local normal haplotypes
            //
            LOG(logINFO) << "constructing local normal haplotypes..." << std::endl;
            std::vector<Haplotype> normalHaps;
            try {
                normalHaps = hapBuilder.constructNormalHaplotypes(window.target, closeVariants, chr,
                                                                  allTumorReads, allNormalReads,
                                                                  &result.MECScore, &result.MECnumReads);
            } catch (std::string &s) {
                LOG(logERROR) << "something strange had happpend during custructing normal haplotypes" << std::endl;
                if (s == "error_making_haplotype") {
                    LOG(logERROR) << "try again withoug close germline variants." << std::endl;
                    CandidateWindow new_window = window;
                    new_window.closeVariants.clear();
                    return call(new_window);
                } else {
                    throw s;
                }
            }

            //
            // filter by mapping quality
            //
            std::vector<IRead *> tumorReads = ReadUtils::filterByMapQuality(tmpTumorReads, params.mapQualThreshold);
            std::vector<IRead *> normalReads = ReadUtils::filterByMapQuality(tmpNormalReads, params.mapQualThreshold);
            ReadsProfile filteredTumorReadsProfile(tumorReads, window.target);
            ReadsProfile filteredNormalReadsProfile(normalReads, window.target);
            LOG(logINFO) << "profiles of reads after mapping-quality filter:" << std::endl;
            LOGP(logINFO) << "tumor reads: " << filteredTumorReadsProfile << std::endl;
            LOGP(logINFO) << "normal reads: " << filteredNormalReadsProfile << std::endl;


            if (tumorReads.size() < params.minReads || normalReads.size() < params.minReads) {
                LOG(logERROR) << "too few number of reads: tumor = " << tumorReads.size()
                << ", normal = " << normalReads.size() << std::endl;
                LOG(logERROR) << "maybe, you should consider check mapQualThreshold" << std::endl;
                throw std::string("too_few_reads");
            }

            //
            // add a somatic haplotype and an error haplotype
            //
            LOG(logINFO) << "add somatic & error haplotypes..." << std::endl;
            std::vector<Haplotype> haps;
            try {
                haps = hapBuilder.addSomaticAndErrorHaplotypes(normalHaps,
                                                               window.target, closeVariants, chr,
                                                               allTumorReads, allNormalReads);
            } catch (std::string &s) {
                LOG(logERROR) << "something strange had happpend during custructing somatic and error haplotypes" << std::endl;
                if (s == "error_making_haplotype") {
                    LOG(logERROR) << "try again withoug close germline variants." << std::endl;
                    CandidateWindow new_window = window;
                    new_window.closeVariants.clear();
                    return call(new_window);
                } else {
                    throw s;
                }
            }

            //
            // calculate alignment likelihoods of reads against haplotypes
            //
            LOG(logINFO) << "calculating alignment likelihoods..." << std::endl;
            std::vector<std::vector<double> > tumorLiks = ProfileHMMUtils::calcLikelihoods(haps, tumorReads);
            std::vector<std::vector<double> > normalLiks = ProfileHMMUtils::calcLikelihoods(haps, normalReads);

            //
            // calculate Bayes factor using a heterozygous germline variant nearby
            //
            if (!(haps[0].variants.empty() && haps[1].variants.empty())) {
                try {
                    HaplotypeBuilder::checkAndSort(&haps, &normalLiks, &tumorLiks,
                                                   normalReads, tumorReads,
                                                   window.target,
                                                   &result.numHaplotypeInformativeReads);
                    MutationModel mm = MutationModel(haps, tumorLiks, normalLiks, params.hap3_params, "mutation");
                    double mmr = mm.estimate();
                    MutationModel em = MutationModel(haps, tumorLiks, normalLiks, params.hap3_params, "error");
                    double emr = em.estimate();
                    result.setBayesFactor(mmr - emr);
                    result.usedCloseGermlineVariants1 = haps[0].variants;
                    result.usedCloseGermlineVariants2 = haps[1].variants;
                    LOG(logINFO) << "Bayes factor: " << result.getBayesFactor() << std::endl;
                } catch (const std::exception& ex) {
                    LOG(logERROR) << ex.what() << std::endl;
                } catch (const std::string& s) {
                    if (s == "umbalanced_normal_haplotypes") {
                        result.setStatus(s);
                    } else if (s == "too_few_haplotype_informative_variant_supporting_reads") {
                        result.setStatus(s);
                    }
                    LOG(logERROR) << s << std::endl;
                } catch (...) {
                    LOG(logERROR) << "something wrong during calculation of bayes factor" << std::endl;
                }
            }

            // prepare for Basic HapmucScore
            std::vector<Haplotype> hapsForBasicTmp = hapBuilder.prepareHaplotypesForBasic(chr, window.target);
            std::vector<std::vector<double> > tumorLiksForBasicTmp = ProfileHMMUtils::calcBasicLikelihoods(hapsForBasicTmp, tumorReads);
            std::vector<std::vector<double> > normalLiksForBasicTmp = ProfileHMMUtils::calcBasicLikelihoods(hapsForBasicTmp, normalReads);
            std::vector<std::vector<double> > tumorLiksForBasic = HaplotypeBuilder::convertLiks2basic(tumorLiksForBasicTmp);
            std::vector<std::vector<double> > normalLiksForBasic = HaplotypeBuilder::convertLiks2basic(normalLiksForBasicTmp);
            std::vector<Haplotype> hapsForBasic = HaplotypeBuilder::convertHaps2basic(hapsForBasicTmp);

            LOG(logDEBUG) << "liks for basic" << std::endl;
            for (int r = 0; r < tumorLiksForBasic.size();r++) {
                for (int x = 0; x < tumorLiksForBasic[r].size(); x++) {
                    LOGP(logDEBUG) << tumorLiksForBasic[r][x] << " ";
                }
                LOGP(logDEBUG) << tumorReads[r]->getSeqName() << std::endl;
            }

            LOG(logDEBUG) << "liks for normal" << std::endl;

            for (int r = 0; r < normalLiksForBasic.size();r++) {
                for (int x = 0; x < normalLiksForBasic[r].size(); x++) {
                    LOGP(logDEBUG) << normalLiksForBasic[r][x] << " ";
                }
                LOGP(logDEBUG) << normalReads[r]->getSeqName() << std::endl;
            }

            //
            // calculate Bayes factor without heterozygous germline variants nearby
            //
            try {
                MutationModel mm2 = MutationModel(hapsForBasic, tumorLiksForBasic, normalLiksForBasic, params.hap2_params, "mutation");
                double mmr2 = mm2.estimate();
                MutationModel em2 = MutationModel(hapsForBasic, tumorLiksForBasic, normalLiksForBasic, params.hap2_params, "error");
                double emr2 = em2.estimate();
                result.setBayesFactorBasic(mmr2 - emr2);
                LOG(logINFO) << "Bayes factor (without germline variants): " << result.getBayesFactorBasic() << std::endl;
            } catch (...) {
                LOG(logERROR) << "something wrong during calculation of bayes factor basic" << std::endl;
            }
        }

    } catch (std::string &s) {
        if (s == "too_many_reads_in_window") {
            result.setStatus(s, true);
        } else if (s == "too_few_reads") {
            result.setStatus(s, true);
        }
        LOG(logERROR) << s << std::endl;
        goto FINALLY; // C++ does not have the 'finally' statement.
    }

FINALLY:
    //
    // fin.
    //
    for (std::vector<IRead *>::iterator it = allTumorReads.begin();it != allTumorReads.end();it++) {
        delete *it;
    }
    for (std::vector<IRead *>::iterator it = allNormalReads.begin();it != allNormalReads.end();it++) {
        delete *it;
    }
    return result;
}
