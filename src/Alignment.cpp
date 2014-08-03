/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <algorithm>
#include "Alignment.h"
#include "log.h"

Alignment::Alignment(const Haplotype &hap,
                     std::string &_seq,
                     std::vector<double> &_baseQualities,
                     int _leftMostPos,
                     double _mapQuality,
                     bool _hasIndel,
                     bool _onReverse,
                     std::string _seqName) : hap(hap) {
    readSeq = _seq;
    baseQualities = _baseQualities;
    leftMostPos = _leftMostPos;
    mapQuality = _mapQuality;
    seqName = _seqName;
    hasIndel = _hasIndel;
    onReverse = _onReverse;
}

void Alignment::print() {
    LOG(logDEBUG) << "lik: "  << likelihood;
    LOGP(logDEBUG) << ", Hap: " << hap.getStartPosInGenome() << " " << hap.getEndPosInGenome() << std::endl;
    LOG(logDEBUG) << "Read(" << seqName << ") " << leftMostPos << " mapQual " << mapQuality << " containIndel " << hasIndel <<  std::endl;
    LOG(logDEBUG) << "onReverse: "  << onReverse << std::endl;
    for (int s = 0; s < history.size(); s++) {
        int x = history[s].first;
        LOGP(logDEBUG) << x << " ";
    }
    LOGP(logDEBUG) << std::endl;
    for (int s = 0; s < history.size(); s++) {
        int i = history[s].second;
        LOGP(logDEBUG) << i;
    }
    LOGP(logDEBUG) << std::endl;
    std::vector<std::string> alignmentStrings = getAlignmentString();
    for (int s = 0; s < alignmentStrings.size(); s++) {
        LOGP(logDEBUG) << alignmentStrings[s];
    }
    LOGP(logDEBUG) << std::endl;
    if (hap.variants.empty()) {
        std::vector<Variant> variants = getAllVariants();
        for (int i = 0; i < variants.size(); i++) {
            LOGP(logDEBUG) << variants[i];
        }
        LOGP(logDEBUG) << std::endl;
    }
}

// convert state X=x in haplotype to genome scale position
int Alignment::getPosXInGenome(int x) {
    if (onReverse) {
        return ((int)getHaplotypeSeq().size() - x) + hap.getStartPosInGenome();
    } else {
        return x + hap.getStartPosInGenome();
    }
}

Variant Alignment::makeVariant(int leftX, int rightX, std::string ref, std::string obs) {
    int leftInGenome = std::min(getPosXInGenome(leftX), getPosXInGenome(rightX));
    if (onReverse) {
        std::reverse(ref.begin(), ref.end());
        std::reverse(obs.begin(), obs.end());
    }
    return Variant(makeVariantString(ref, obs), leftInGenome);
}

std::string Alignment::makeVariantString(std::string ref, std::string obs) {
    std::stringstream ss;
    if (ref == "-") {
        ss << "+" << obs;
    } else if (obs == "-") {
        ss << "-" << ref;
    } else {
        ss << ref << "=>" << obs;
    }
    return ss.str();
}

std::vector<Variant> Alignment::getAllVariants() {
    if (!hap.variants.empty()) {
        LOG(logERROR) << "getAllVariants only work for haplotypes, which do not have variants." << std::endl;
        return std::vector<Variant>();
    }
    std::vector<Variant> variants;
    int s = 0;
    while (s < getReadSeq().size()) {
        std::pair<int, int> &xi = history[s];
        if (xi.second == 1) {
            int i = s;
            // find continuous insertion
            while (i < getReadSeq().size()) {
                if (history[i].second != 1) break;
                i++;
            }
            variants.push_back(makeVariant(xi.first + 1, xi.first + 1, "-", getReadSeq().substr(s, i - s)));
            s = i + 1;
        } else {
            if (s != 0 && xi.first - history[s-1].first != 1) {
                // deletion
                int leftX = history[s-1].first + 1, rightX = xi.first;
                int leftG = getPosXInGenome(leftX), rightG = getPosXInGenome(rightX);
                LOG(logDEBUG) << leftX << " " << rightX << " " << leftG << " " << rightG << std::endl;
                variants.push_back(makeVariant(history[s-1].first + 1, xi.first, getHaplotypeSeq().substr(history[s-1].first + 1, xi.first - history[s-1].first - 1), "-"));
            }
            // REF or SNV
            std::string rbs = getReadSeq().substr(s, 1);
            // convert char to std::string
            char hb = getHapBase(xi.first);
            std::string hbs;
            std::stringstream ss;
            ss << hb;
            hbs = ss.str();
            // if haplotype base and read base did not match, add SNV
            if (hbs != rbs) {
                variants.push_back(makeVariant(xi.first, xi.first + 1, hbs, rbs));
            }
            s++;
        }
    }
    return variants;
}


inline std::string Alignment::getReadSeq() {
    return readSeq;
}

inline char Alignment::getHapBase(int x) {
    if (x == getHaplotypeSeq().size()) {
        return 'R';
    } else {
        return getHaplotypeSeq()[x];
    }
}

inline std::string Alignment::getHaplotypeSeq() {
    if (hapSeq != "") return hapSeq;
    if (!onReverse) {
        hapSeq = hap.seq;
        return hap.seq;
    } else {
        std::string tmp = hap.seq;
        std::reverse(tmp.begin(), tmp.end());
        hapSeq = tmp;
        return hapSeq;
    }
}

std::string Alignment::getHistoryBases(int b, int x1, int x0, int i) {
    char hb = getHapBase(x1), rb = getReadSeq()[b];
    if (i == 1) {
        return "I";
    } else {
        if (b != 0 && x1 != x0 + 1) {
            std::stringstream ss;
            ss << x1 - x0 - 1 << "D";
            return ss.str();
        } else {
            if (hb == rb) {
                return "-";
            } else {
                return std::string(1, rb);
            }
        }
    }
}

bool Alignment::isAligned() {
    const double th1 = -70, th2 = -120;
    if (likelihood > th1) {
        return true;
    } else {
        if (likelihood < th2) {
            LOG(logWARNING) << likelihood << " < " << th2 << std::endl;
            print();
            return false;
        } else {
            LOG(logDEBUG) << likelihood << " < " << th1 << std::endl;
            print();
            return true;
        }
    }
}

std::vector<std::string> Alignment::getAlignmentString() {
    // store std::strings for maximum-likelihood alignment
    std::vector<std::string> history_bases;
    history_bases.reserve(getReadSeq().size());
    for (int s = 0; s < getReadSeq().size(); s++) {
        if (s == 0) {
            history_bases.push_back(getHistoryBases(s, history[0].first, -1, history[0].second));
        } else {
            history_bases.push_back(getHistoryBases(s, history[s].first, history[s - 1].first, history[s].second));
        }
    }
    return history_bases;
}
