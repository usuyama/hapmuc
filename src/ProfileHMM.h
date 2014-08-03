/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef __ProfileHMM__
#define __ProfileHMM__

#include <utility>
#include <vector>
#include <string>
#include "ReadIndelErrorModel.hpp"
#include "SamRead.h"
#include "Haplotype.hpp"
#include "Alignment.h"

class ProfileHMM {
public:
    ProfileHMM(const Haplotype &_hap,
               const std::string &_seq,
               const std::vector<double> &_baseQualities,
               int _leftMostPos,
               double _mapQuality,
               bool _hasIndel,
               bool _onReverse,
               std::string _seqName);
    Alignment viterbi();
    static int maxIndelLength;
    static int maxHomopolymerLength;
    static double priorDelByError;
private:
    int size_X;
    int leftMostPos;
    double mapQuality;
    bool hasIndel;
    std::string seqName;
    std::string hapSeq, readSeq;
    std::vector<double> baseQualities;
    std::vector<std::vector<int> > historyStoreX; // x1, i1 -> x0
    std::vector<std::vector<int> > historyStoreI; // x1, i1 -> i0
    std::vector<double> currentDelta, prevDelta;
    ReadIndelErrorModel readIndelErrorModel;
    int startSearchWindowSize, searchWindowSize; // for initial guess of state X in b == 0
    int anchor; // left-most index of read in haplotype
    const Haplotype &haplotype;
    bool onReverseStrand;
    int insT; // if 1, do not consider insertion state
    std::vector<double>logInsStartProb; // log(insStartProb(homopolymer_length))
    std::vector<double>one_logInsStartProb; // log(1-insStartProb(homopolymer_length))
    inline int idx(int x, int i);
    inline std::pair<int, int> from_idx(int index);
    inline std::pair<int, int> getHistory(int b, int x, int i); // get previous states (X0, I0) when X1=x, I0=i in position b
    inline void setHistory(int b, int x1, int i1, int x0, int i0); // set history
    inline double logEmit(char h, char r, double bq, int x, int i);
    inline double logTransI(int i1, int i0, int homopolymer_length);
    inline double logTransX(int x1, int x0, int i1, int i0);
    inline char getHapBase(int x);
    inline int getHomopolymerLength(int x);
    void initInsStartProb();
    std::vector<double> logDelProbStore;
    inline double logDelProb(int d); // d == 1 means not having a deletion
    inline std::pair<int, int> getRangeX1(int b); // get search window for state X
    inline std::pair<int, int> getRangeX0(int b0);// get max and min possible state X
    std::string getHistoryBases(int b, int x1, int x0, int i);
    inline int idxDelta(int x, int i, std::vector<double> delta);
    inline std::pair<int, int>from_idxDelta(int idx, int base, std::vector<double> delta);
    std::vector<Variant> getVariants();
};

#endif /* defined(__ProfileHMM__) */
