#ifndef __haplotype_hapser__Read__
#define __haplotype_hapser__Read__

#include <string>
#include <vector>

namespace haplotype_phaser {
class Read {
    int start, end;
    bool _isHaplotypeInformative;
public:
    std::string seq, name;
    bool isHaplotypeInformative() const;
    int getStartIndex() const;
    int getEndIndex() const;
    int distance(std::string &partialHaplotype, int pos) const;
    int minDistance(std::string &partialHaplotype, int pos) const;
    static Read fromString(std::string line, std::string name);
    static int getEndIndex(std::vector<Read *> &reads);
};
}

#endif /* defined(__haplotype_hapser__Read__) */

