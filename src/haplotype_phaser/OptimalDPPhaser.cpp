#include "OptimalDPPhaser.h"

#include "Read.h"
#include "History.h"

#include <iostream>
#include <sstream>

namespace haplotype_phaser {

std::vector<Read *> OptimalDPPhaser::findReadsStartAt(int i) {
    std::vector<Read *> set;
    for (int s = 0;s < readMatrix.size(); s++) {
        if (readMatrix[s].getStartIndex() == i) {
            set.push_back(&readMatrix[s]);
        }
    }
    return set;
}

std::vector<std::string> OptimalDPPhaser::generatePartialHaplotypes(int size) {
    std::vector<std::string> haplotypes;
    haplotypes.push_back("");
    for (int i = 0;i < size;i++) {
        std::vector<std::string> haps1 = haplotypes;
        for (int j = 0;j < haps1.size();j++) {
            haplotypes[j] += '0';
            haps1[j] += '1';
        }
        haplotypes.insert(haplotypes.end(), haps1.begin(), haps1.end());
    }
    return haplotypes;
}

OptimalDPPhaser::OptimalDPPhaser(std::vector<Read> &_readMatrix): readMatrix(_readMatrix) {
}

int OptimalDPPhaser::calcDistance(std::vector<Read *> &reads, std::string &seq, int pos) {
    int distance = 0;
    std::string compSeq = seq;
    for (int i = 0;i< compSeq.length();i++) {
        if (compSeq[i] == '1') {
            compSeq[i] = '0';
        } else if (compSeq[i] == '0') {
            compSeq[i] = '1';
        }
    }
    for (int i = 0;i < reads.size();i++) {
        int score0 = reads[i]->distance(seq, pos);
        int score1 = reads[i]->distance(compSeq, pos);
        distance += std::min(score0, score1);
    }
    return distance;
}
std::string OptimalDPPhaser::toString(std::string &haplotype, int MEC) {
    std::stringstream ss;
    ss << haplotype << std::endl;
    for (int i = 0;i < readMatrix.size();i++) {
        ss << readMatrix[i].seq << " " <<  readMatrix[i].minDistance(haplotype, 0);
        ss << " " << readMatrix[i].name << std::endl;
    }
    ss << "MEC = " << MEC << std::endl;
    return ss.str();
}

void OptimalDPPhaser::printResult(std::string &haplotype, int MEC) {
    std::cout << toString(haplotype, MEC);
}

void OptimalDPPhaser::inferHaplotype(std::string *haplotype, int *MEC) {
    int n = (int)readMatrix[0].seq.length();

    std::vector<History> histories;

    for (int i = 0;i < n;i++) {
        std::vector<Read *> reads = findReadsStartAt(i);
        int endIndex = Read::getEndIndex(reads);
        if (i != 0) {
            int prevEndIndex = histories[i-1].haplotypeSize + i - 2;
            if (prevEndIndex > endIndex) {
                endIndex = prevEndIndex;
            }
        }
        int hapLength = endIndex - i + 1;
        History currentState = History(hapLength);
        std::vector<std::string> partialHaplotypes = generatePartialHaplotypes(hapLength);
        for (int j = 0;j < partialHaplotypes.size();j++) {
            std::string &ph = partialHaplotypes[j];
            int currentScore = calcDistance(reads, ph, i);
            if (i == 0) {
                currentState.setS(ph, currentScore, '-');
            } else {
                std::string h0 = '0' + ph;
                int score0 = histories[i-1].getS(h0);
                std::string h1 = '1' + ph;
                int score1 = histories[i-1].getS(h1);
                if (score0 < score1) {
                    currentState.setS(ph, currentScore + score0, '0');
                } else {
                    currentState.setS(ph, currentScore + score1, '1');
                }
            }
        }
        histories.push_back(currentState);
    }

    // get best state
    std::string bestPartialHaplotype;
    std::pair<std::string, int> bestState = (histories.end()-1)->getMinimumState();
    bestPartialHaplotype = bestState.first;

    // trace the results
    for (int i = (int)histories.size()-1;i > 0;i--) {
        char base = histories[i].getHistory(bestPartialHaplotype);
        bestPartialHaplotype = base + bestPartialHaplotype;
    }
    *haplotype = bestPartialHaplotype;
    *MEC = bestState.second;
}
}

