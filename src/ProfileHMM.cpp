/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <float.h>
#include <cmath>
#include <algorithm>
#include "log.h"
#include "ReadIndelErrorModel.hpp"
#include "ProfileHMM.h"
#include "Haplotype.hpp"
#include "SamRead.h"
#include "Alignment.h"

// initialize class variables
int ProfileHMM::maxIndelLength = 25;
int ProfileHMM::maxHomopolymerLength = 30;
double ProfileHMM::priorDelByError = 5e-4;

// returns index for vectors (historyStoreX, historyStoreI, currentDelta, prevDelta)
inline int ProfileHMM::idx(int x, int i) {
    return x * 2 + i;
}

// retrieve pair of X and I from index
inline std::pair<int, int> ProfileHMM::from_idx(int index) {
    return std::make_pair(index / 2, index % 2);
}

ProfileHMM::ProfileHMM(const Haplotype &_hap,
                       const std::string &_seq,
                       const std::vector<double> &_baseQualities,
                       int _leftMostPos,
                       double _mapQuality,
                       bool _hasIndel,
                       bool _onReverse,
                       std::string _seqName) : haplotype(_hap) {
    LOG(logDEBUG2) << "initialize ProfileHMM:" << std::endl;
    mapQuality = _mapQuality;
    seqName = _seqName;
    leftMostPos = _leftMostPos;
    hasIndel = _hasIndel;
    // initialize haplotype and read
    if (_onReverse) {
        // reverse both haplotype and read because
        // we'd like to consider the direction that the sequencer generated this read
        onReverseStrand = true;
        hapSeq = _hap.seq;
        reverse(hapSeq.begin(), hapSeq.end());
        readSeq = _seq;
        reverse(readSeq.begin(), readSeq.end());
        baseQualities = _baseQualities;
        reverse(baseQualities.begin(), baseQualities.end());
        anchor = haplotype.getStartPosInGenome() + (int)hapSeq.length() - (leftMostPos + (int)readSeq.length());
    } else {
        onReverseStrand = false;
        hapSeq = _hap.seq;
        readSeq = _seq;
        baseQualities = _baseQualities;
        anchor = leftMostPos - haplotype.getStartPosInGenome();
    }

    // initialize memories for X and I
    size_X = (int)hapSeq.size() + 1; // Haplotype positions + 'R' state
    historyStoreX.reserve(readSeq.size());
    historyStoreI.reserve(readSeq.size());
    for (int i = 0; i < readSeq.size(); i++) {
        historyStoreX.push_back(std::vector<int>(size_X * 2));
        historyStoreI.push_back(std::vector<int>(size_X * 2));
        fill(historyStoreX[i].begin(), historyStoreX[i].end(), -1);
        fill(historyStoreI[i].begin(), historyStoreI[i].end(), -1);
    }
    currentDelta = std::vector<double>(size_X * 2); // TODO: it can be optimized
    prevDelta = std::vector<double>(size_X * 2);
    
    // initialize probabilities for deletion
    logDelProbStore = std::vector<double>(maxIndelLength + 2);
    logDelProbStore[1] = log(1 - priorDelByError);
    double norm = 0.0;
    for (int x = 0; x < maxIndelLength + 2; x++) if (x != 1) {
        double tmp = -fabs(x - 1);
        logDelProbStore[x] = tmp;
        norm += exp(tmp);
    }
    norm = log(norm / priorDelByError);
    for (int x = 0; x < logDelProbStore.size(); x++) if (x != 1) {
        logDelProbStore[x] -= norm;
    }

    // setting search window size
    // we will check states in the range [anchor-window, anchor+window]
    // using anchor = the position written in bam file
    // if there indels exist, the window is rather bigger and insertion state (I=1) is considered
    if (haplotype.hasIndel || hasIndel) {
        startSearchWindowSize = maxIndelLength;
        searchWindowSize = maxIndelLength;
        insT = 2;
    } else {
        startSearchWindowSize = 4;
        searchWindowSize = 4;
        insT = 1;
    }
    
    initInsStartProb();
}

// get the transition history of read position b, X=x, insertion state I=i
inline std::pair<int, int> ProfileHMM::getHistory(int b, int x, int i) {
    return std::pair<int, int> (historyStoreX[b][idx(x, i)], historyStoreI[b][idx(x, i)]);
}

inline void ProfileHMM::setHistory(int b, int x1, int i1, int x0, int i0) {
    historyStoreX[b][idx(x1, i1)] = x0;
    historyStoreI[b][idx(x1, i1)] = i0;
}

// returns haplotype base of state X=x.
// Note: we don't have the L state.
inline char ProfileHMM::getHapBase(int x) {
    if (x == hapSeq.size()) {
        return 'R';
    } else {
        return hapSeq[x];
    }
}

// returns the length of homopolymer in haplotype at state X=x
inline int ProfileHMM::getHomopolymerLength(int x) {
    int c = 1;
    const char hb = getHapBase(x);
    if (hb == 'R' || hb == 'N') {
        return c;
    }
    for (int i = x - 1; i >= 0; i--) {
        if (hb == getHapBase(i)) {
            c++;
        } else {
            break;
        }
    }
    if (c >= logInsStartProb.size()) {
        LOG(logDEBUG2) << "long homopolymer (x, homopolymer_length) " << x << " " << c << std::endl;
        c = (int)logInsStartProb.size() - 1;
    }
    return c;
}

// returns the log probability of emission 
inline double ProfileHMM::logEmit(char hapBase, char readBase, double bq, int x , int i) {
    if (i == 1) {
        return log(mapQuality * bq);
    }
    if (x == 'R') {
        return log(mapQuality * bq);
    }
    if (hapBase == readBase) {
        return log(mapQuality * bq + (1.0 - mapQuality * bq) / 4.0);
    } else {
        return log((1.0 - mapQuality * bq) / 4.0);
    }
}

void ProfileHMM::initInsStartProb() {
    logInsStartProb.reserve(maxHomopolymerLength + 1);
    one_logInsStartProb.reserve(maxHomopolymerLength + 1);
    for (int i = 0; i < maxHomopolymerLength + 1; i++) {
        double d = readIndelErrorModel.getViterbiHPError(i);
        logInsStartProb.push_back(log(d));
        one_logInsStartProb.push_back(log(1.0 - d));
    }
}

inline double ProfileHMM::logTransI(int i1, int i0, int homopolymer_length) {
    if (i0 == 0) {
        if (i1 == 1) {
            return logInsStartProb[homopolymer_length];
        } else {
            return one_logInsStartProb[homopolymer_length];
        }
    } else {
        if (i1 == 1) {
            return -1.0;
        } else {
            return log(1 - exp(-1.0));
        }
    }
}

inline double ProfileHMM::logDelProb(int d) {
    if (d > logDelProbStore.size() - 1) {
        return *(logDelProbStore.end() - 1);
    } else {
        return logDelProbStore[d];
    }
}

inline double ProfileHMM::logTransX(int x1, int x0, int i1, int i0) {
    if (x1 < x0) {
        return -DBL_MAX;
    }
    if (i0 == 0) {
        if (i1 == 0) {
            return logDelProb(abs(x1 - x0));
        } else {
            return (x1 == x0) ? 0.0 : -DBL_MAX;
        }
    } else {
        if (i1 == 0) {
            return (x0 + 1 == x1) ? 0.0 : -DBL_MAX;
        } else {
            return (x0 == x1) ? 0.0 : -DBL_MAX;
        }
    }
}

inline std::pair<int, int> ProfileHMM::getRangeX1(int b) { // get search window for state X
    return std::pair<int, int> (anchor - startSearchWindowSize + b, anchor + startSearchWindowSize + b + 1);
}

// for checking the result
// ex. ---II--G---3D---
std::string ProfileHMM::getHistoryBases(int b, int x1, int x0, int i) {
    char hb = getHapBase(x1), rb = readSeq[b];
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

// perform viterbi algorithm to get maximum alignment and its likelihood
Alignment ProfileHMM::viterbi() {
    LOG(logDEBUG2) << "initialize viterbi:" << std::endl;
    // start from read position b0 = 0 (left most)
    int b0 = 0;
    std::pair<int, int> currentRangeX, prevRangeX;
    
    // search window of possible states at b0 = 0
    currentRangeX = getRangeX1(b0);
    
    // calculate probabilities for b0 = 0 and store them in vector prevDelta
    fill(prevDelta.begin(), prevDelta.end(), -DBL_MAX);
    for (int x = currentRangeX.first; x < currentRangeX.second; x++) {
        for (int i = 0; i < insT; i++) {
            prevDelta[idx(x, i)] = logEmit(getHapBase(x), readSeq[b0], baseQualities[b0], x, i) + logTransI(i, 0, getHomopolymerLength(x - 1));
        }
    }

    //
    // for every read-base, search best states (X, I) and the transisions.
    //
    for (int b1 = b0 + 1; b1 < readSeq.size(); b1++) {
        // calculate maximum-liklihood states in position b1
        prevRangeX = currentRangeX;
        currentRangeX = getRangeX1(b1); // search window for b1
        fill(currentDelta.begin(), currentDelta.end(), -DBL_MAX);
        for (int x1 = currentRangeX.first; x1 < currentRangeX.second; x1++) {
            int homopolymerLength = getHomopolymerLength(x1 - 1);
            for (int i1 = 0; i1 < insT; i1++) {
                double e = logEmit(getHapBase(x1), readSeq[b1], baseQualities[b1], x1, i1);
                // calculate maximum-liklihood for X=X1, I=i1
                double maxLogLik = -DBL_MAX;
                int max_x0, max_i0;
                for (int x0 = prevRangeX.first; x0 < prevRangeX.second; x0++) {
                    for (int i0 = 0; i0 < insT; i0++) {
                        double ti = logTransI(i1, i0, homopolymerLength);
                        double tx = logTransX(x1, x0, i1, i0);
                        double prev = prevDelta[idx(x0, i0)];
                        double logLik = e + ti + tx + prev;
                        if (logLik > 0) {
                            throw std::string("bug in profileHMM");
                        }
                        if (logLik > maxLogLik) {
                            max_x0 = x0;
                            max_i0 = i0;
                            maxLogLik = logLik;
                        }
                    }
                }
                // store the best transition for X=X1, I=i1, and its likelihood
                currentDelta[idx(x1, i1)] = maxLogLik;
                setHistory(b1, x1, i1, max_x0, max_i0);
            }
        }
        swap(prevDelta, currentDelta);
    }
    
    //
    // get maximum alignment from transition histories
    // first, get maximum likelihood state for the last base
    // 
    std::vector<double>::iterator max_it = max_element(prevDelta.begin(), prevDelta.end());
    int lastIndex = (int)(max_it - prevDelta.begin());

    // store maximum-likelihood states X, I in vector history
    // NOTE: backward
    std::vector<std::pair<int, int> > history;
    history.reserve(readSeq.size());
    history.push_back(from_idx(lastIndex));
    for (int s = 0; s < readSeq.size() - 1; s++) {
        history.push_back(getHistory((int)readSeq.size() - s - 1, history[s].first, history[s].second));
    }
    reverse(history.begin(), history.end());
    
    Alignment result(haplotype, readSeq, baseQualities, leftMostPos, mapQuality, hasIndel, onReverseStrand, seqName);
    result.likelihood = *max_it;
    result.history = history;
    result.onReverse = onReverseStrand;
    return result;
}
