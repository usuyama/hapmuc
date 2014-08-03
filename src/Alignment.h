/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef __Alignment__
#define __Alignment__

#include <vector>
#include <utility>
#include "Haplotype.hpp"
#include "SamRead.h"

class Alignment {
    std::string getReadSeq();
    std::string getHaplotypeSeq();
    std::string getHistoryBases(int b, int x1, int x0, int i);
    std::string hapSeq, readSeq;
    inline char getHapBase(int x);
    std::string makeVariantString(std::string ref, std::string obs);
    int getPosXInGenome(int x);
    Variant makeVariant(int leftX, int rightX, std::string ref, std::string obs);
public:
    Alignment(const Haplotype &hap,
              std::string &_seq,
              std::vector<double> &_baseQualities,
              int _leftMostPos,
              double _mapQuality,
              bool _hasIndel,
              bool _onReverse,
              std::string _seqName);
    std::string seqName;
    std::vector<double> baseQualities;
    int leftMostPos;
    bool hasIndel;
    double mapQuality;
    const Haplotype &hap;
    double likelihood;
    bool onReverse;
    bool isAligned();
    std::vector<std::pair<int, int> > history; // pair of X, I
    void print();
    std::vector<std::string> getAlignmentString();
    std::vector<Variant> getAllVariants();
};

#endif /* defined(__Alignment__) */
