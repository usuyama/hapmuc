#ifndef __haplotype_hapser__OptimalDPPhaser__
#define __haplotype_hapser__OptimalDPPhaser__

#include "Read.h"

#include <vector>
#include <string>

namespace haplotype_phaser {

//
// implements a dynamic-programming assembler for haplotype-phasing
// http://bioinformatics.oxfordjournals.org/content/26/12/i183.full
//
class OptimalDPPhaser {
public:
    std::vector<Read> &readMatrix;
    OptimalDPPhaser(std::vector<Read> &_readMatrix);
    std::vector<Read *> findReadsStartAt(int i);
    std::vector<std::string> generatePartialHaplotypes(int size);
    void inferHaplotype(std::string *haplotype, int *MEC);
    int calcDistance(std::vector<Read *> &reads, std::string &seq, int pos);
    void printResult(std::string &haplotype, int MEC);
    std::string toString(std::string &haplotype, int MEC);
};
}

#endif /* defined(__haplotype_hapser__OptimalDPPhaser__) */

