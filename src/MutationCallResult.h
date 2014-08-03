/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef MutationCallResult_h
#define MutationCallResult_h

#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

#include "CandidateWindow.h"
#include "Variant.h"
#include "VariantUtils.h"

class MutationCallResult {
    double bayesFactor, bayesFactorBasic;
    std::string status;
public:
    double avgMapQualTumorObs, avgMapQualTumorRef, avgMapQualNormalObs, avgMapQualNormalRef;
    double softClipFreqTumorObs, softClipFreqTumorRef, softClipFreqNormalObs, softClipFreqNormalRef;
    double indelFreqTumorObs, indelFreqTumorRef, indelFreqNormalObs, indelFreqNormalRef;
    int MECScore, MECnumReads;
    int numHaplotypeInformativeReads;
    bool hasBayesFactor, hasBayesFactorBasic;
    void setBayesFactor(double bf) {
        bayesFactor = bf;
        hasBayesFactor = true;
    }
    
    void setBayesFactorBasic(double bfb) {
        bayesFactorBasic = bfb;
        hasBayesFactorBasic = true;
    }
    
    double getBayesFactor() {
        if (hasBayesFactor) {
            return bayesFactor;
        } else {
            throw std::string("not calculated");
        }
    }
    
    double getBayesFactorBasic() {
        if (hasBayesFactorBasic) {
            return bayesFactorBasic;
        } else {
            throw std::string("not calculated");
        }
    }
    
    void setStatus(std::string s, bool force = false) {
        if (status == "" || force) {
            status = s;
        }
    }
    
    std::string getStatus() {
        if (status == "") {
            return "considered_close_hetero_germline_variants";
        } else {
            return status;
        }
    }
    
    std::vector<Variant> usedCloseGermlineVariants1, usedCloseGermlineVariants2;
    CandidateWindow candidateWindow;
    MutationCallResult(CandidateWindow candidateWindow_) {
        candidateWindow = candidateWindow_;
		hasBayesFactor = false;
		hasBayesFactorBasic = false;
        status = "";
        avgMapQualTumorRef = -1;
        avgMapQualTumorObs = -1;
        avgMapQualNormalRef = -1;
        avgMapQualNormalObs = -1;
        softClipFreqTumorRef = -1;
        softClipFreqTumorObs = -1;
        softClipFreqNormalRef = -1;
        softClipFreqNormalObs = -1;
        indelFreqTumorRef = -1;
        indelFreqTumorObs = -1;
        indelFreqNormalRef = -1;
        indelFreqNormalObs = -1;
        MECScore = -1;
        MECnumReads = -1;
        numHaplotypeInformativeReads = -1;
	}
    
    static std::string getHeader() {
        std::string header = "chr\tstart\tend\tref\tobs\t";
        header += "ref_count_tumor\tobs_count_tumor\tref_count_normal\tobs_count_normal\t";
        header += "missrate_tumor\tstrandrate_tumor\tmissrate_normal\tstrandrate_normal\t";
        header += "ref_bq_tumor\tobs_bq_tumor\tref_bq_normal\tobs_bq_normal\t";
        header += "triallelic_site_check\tindel_cover_check\t";
        header += "fisher\t";
        header += "bf_without_hetero_germline_variants\tbf_with_hetero_germline_variants\thapmuc_score\tstatus\t";
        header += "germline_snp_nearby1\tgermline_snp_nearby2\t";
        header += "distance\t";
        header += "num_haplotype_informative_variant_supporting_reads\t";
        header += "avg_mapping_quality_tumor_ref\tavg_mapping_quality_tumor_obs\tavg_mapping_quality_normal_ref\tavg_mapping_quality_normal_obs\t";
        header += "softclip_freq_tumor_ref\tsoftclip_freq_tumor_obs\tsoftclip_freq_normal_ref\tsoftclip_freq_normal_obs\t";
        header += "indel_freq_tumor_ref\tindel_freq_tumor_obs\tindel_freq_normal_ref\tindel_freq_normal_obs\t";
        header += "MEC_score\tMEC_num_reads\t";

        return header;
    }
    
    template <class X> static inline void push(std::stringstream &ss, X x) {
        ss << x << "\t";
    }

    template <class X> static inline void pushValue(std::stringstream &ss, X x) {
        if (std::isnan(x) || x < 0) {
            ss << "-" << "\t";
        } else {
            if (std::abs(x - 0.0) < 1e-4) {
                ss << 0 << "\t";
            } else if (std::abs(x - 1.0) < 1e-4) {
                ss << 1 << "\t";
            } else {
                ss << x << "\t";
            }
        }
    }
    
    std::string getOutput() {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(4);
        const CandidateWindow::Information &vi = candidateWindow.info;
        push(ss, vi.chr);
        push(ss, vi.start);
        push(ss, vi.end);
        push(ss, vi.ref);
        push(ss, vi.obs);
        push(ss, vi.ref_count_tumor);
        push(ss, vi.obs_count_tumor);
        push(ss, vi.ref_count_normal);
        push(ss, vi.obs_count_normal);
        push(ss, vi.missrate_tumor);
        push(ss, vi.strandrate_tumor);
        push(ss, vi.missrate_normal);
        push(ss, vi.strandrate_normal);
        push(ss, vi.ref_bq_tumor);
        push(ss, vi.obs_bq_tumor);
        push(ss, vi.ref_bq_normal);
        push(ss, vi.obs_bq_normal);
        push(ss, vi.triallelic_site_check);
        push(ss, vi.indel_cover_check);
        push(ss, vi.fisher_score);
        if (hasBayesFactorBasic) {
            push(ss, bayesFactorBasic);
            if (hasBayesFactor) {
                push(ss, bayesFactor); // bf_with_snp
                push(ss, bayesFactor); // hapmuc_score
                push(ss, getStatus());
                push(ss, VariantUtils::getSymbols(usedCloseGermlineVariants1));
                push(ss, VariantUtils::getSymbols(usedCloseGermlineVariants2));
                std::vector<Variant> allVariants = usedCloseGermlineVariants1;
                allVariants.insert(allVariants.end(), usedCloseGermlineVariants2.begin(), usedCloseGermlineVariants2.end());
                push(ss, VariantUtils::getMinDistance(candidateWindow.target, allVariants));
                push(ss, numHaplotypeInformativeReads);
            } else {
                push(ss, "-");
                push(ss, bayesFactorBasic); //hapmuc_score
                push(ss, getStatus());
                push(ss, "-");
                push(ss, "-");
                push(ss, "-");
                push(ss, "-");
            }
        } else {
            // do not have both bf_without_snp or bf_with_snp
            push(ss, "-"); // bf_without_snp
            push(ss, "-"); // bf_with_snp
            push(ss, "-"); // hapmuc_score
            push(ss, getStatus());
            push(ss, "-"); // variants1
            push(ss, "-"); // variants2
            push(ss, "-"); // distance
            push(ss, "-");
        }
        pushValue(ss, avgMapQualTumorRef);
        pushValue(ss, avgMapQualTumorObs);
        pushValue(ss, avgMapQualNormalRef);
        pushValue(ss, avgMapQualNormalObs);
        pushValue(ss, softClipFreqTumorRef);
        pushValue(ss, softClipFreqTumorObs);
        pushValue(ss, softClipFreqNormalRef);
        pushValue(ss, softClipFreqNormalObs);
        pushValue(ss, indelFreqTumorRef);
        pushValue(ss, indelFreqTumorObs);
        pushValue(ss, indelFreqNormalRef);
        pushValue(ss, indelFreqNormalObs);
        pushValue(ss, MECScore);
        if (MECnumReads < 0) {
            ss << "-";
        } else {
            ss << MECnumReads;
        }
        return ss.str();
    }
};

#endif
