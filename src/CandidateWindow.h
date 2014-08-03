/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef __CandidateWindow__
#define __CandidateWindow__

#include <vector>
#include "Variant.h"

class CandidateWindow {
    std::vector<Variant> parseCloseVariants(std::string str);
    std::string getVariantSymbol(std::string ref, std::string obs);
public:
    struct Information {
        std::string chr;
        int start, end;
        std::string ref, obs;
        std::string ref_count_tumor, obs_count_tumor;
        std::string missrate_tumor, strandrate_tumor;
        std::string ref_count_normal, obs_count_normal;
        std::string missrate_normal, strandrate_normal;
        std::string ref_bq_tumor, obs_bq_tumor, ref_bq_normal, obs_bq_normal;
        std::string triallelic_site_check, indel_cover_check;
        std::string fisher_score;
      //  std::string estimated_error_rate;
    };
    Information info;
    Variant target;
    std::vector<Variant> closeVariants;
    CandidateWindow(std::string line);
    CandidateWindow() {
        target = Variant();
        closeVariants = std::vector<Variant>();
    }
    CandidateWindow(const CandidateWindow&); //copy constructor
    CandidateWindow& operator=(const CandidateWindow&);
};

#endif /* defined(__CandidateWindow__) */
