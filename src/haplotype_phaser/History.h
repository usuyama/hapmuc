#ifndef __haplotype_hapser__History__
#define __haplotype_hapser__History__

#include <string>
#include <map>

namespace haplotype_phaser {
class History {
public:
    int haplotypeSize;
    std::map<std::string, int> scores;
    std::map<std::string, char> history;
    History(int _haplotypeSize) {
        haplotypeSize = _haplotypeSize;
    }

    inline int getS(std::string &partialHaplotype) {
        std::string hap = partialHaplotype.substr(0, haplotypeSize);
        return scores[hap];
    }

    inline void setS(std::string &partialHaplotype, int score, char prevBase) {
        scores[partialHaplotype] = score;
        history[partialHaplotype] = prevBase;
    }

    inline char getHistory(std::string &partialHaplotype) {
        std::string hap = partialHaplotype.substr(0, haplotypeSize);
        return history[hap];
    }

    std::pair<std::string, int> getMinimumState();
};
}


#endif /* defined(__haplotype_hapser__History__) */
