/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <cstdlib>
#include <sstream>
#include <algorithm>
#include "CandidateWindow.h"
#include "Utils.h"

// ex: "1495:G=>C", "1495:G=>C,1410:+CGG", "1495:-AA"
std::vector<Variant> CandidateWindow::parseCloseVariants(std::string str) {
    std::vector<Variant> variants;
    if (str == "") {
        throw 1;
    }
    if (str == "-") {
        return variants;
    }
    std::vector<std::string> var;
    std::vector<std::string> variantsInStr = Utils::split(str, ",");
    variants.reserve(variantsInStr.size());
    for (int i = 0; i < variantsInStr.size(); i++) {
        var = Utils::split(variantsInStr[i], ":");
        Variant v(var[1], atoi(var[0].c_str()));
        variants.push_back(v);
    }
    return variants;
}

inline std::string CandidateWindow::getVariantSymbol(std::string ref, std::string obs) {
    if (ref == "-") {
        return ("+" + obs);
    } else if (obs == "-") {
        return ("-" + ref);
    } else {
        return (ref + "=>" + obs);
    }
}

CandidateWindow::CandidateWindow(std::string line) {
    try {
        std::stringstream linestream(line);
        CandidateWindow::Information &vi = info;
        std::string closeGermlineVariants;
        linestream >> vi.chr >> vi.start >> vi.end >> vi.ref >> vi.obs
        >> vi.ref_count_tumor >> vi.obs_count_tumor >> vi.ref_count_normal
        >> vi.obs_count_normal >> vi.missrate_tumor >> vi.strandrate_tumor
        >> vi.missrate_normal >> vi.strandrate_normal >> vi.ref_bq_tumor
        >> vi.obs_bq_tumor >> vi.ref_bq_normal >> vi.obs_bq_normal;
        if ((int)std::count(line.begin(), line.end(), '\t') == 20) {
            linestream >> vi.triallelic_site_check >> vi.indel_cover_check
                       >> vi.fisher_score >> closeGermlineVariants;
        } else {
            linestream >> vi.fisher_score >> closeGermlineVariants;
        }
        target = Variant(getVariantSymbol(vi.ref, vi.obs), vi.start);
        closeVariants = parseCloseVariants(closeGermlineVariants);
    } catch (...) {
        throw std::string("parsing candidate window failed: " + line);
    }
}

CandidateWindow::CandidateWindow(const CandidateWindow &cand) {
    info = cand.info;
    target = cand.target;
    closeVariants = cand.closeVariants;
}

CandidateWindow &CandidateWindow::operator=(const CandidateWindow &cand) {
    info = cand.info;
    target = cand.target;
    closeVariants = cand.closeVariants;
    return *this;
}
