#include "Read.h"

namespace haplotype_phaser {

Read Read::fromString(std::string line, std::string name) {
    Read r;
    r.seq = line;
    r.name = name;
    bool first = true;
    int cover_count = 0; // number of positions that is != '-'
    for (int i = 0;i < r.seq.size(); i++) {
        if (first && r.seq[i] != '-') {
            r.start = i;
            first = false;
        }
        if (r.seq[i] != '-') {
            cover_count++;
            r.end = i;
        }
    }
    if (cover_count >= 2) {
        r._isHaplotypeInformative = true;
    } else {
        r._isHaplotypeInformative = false;
    }
    return r;
}
    
bool Read::isHaplotypeInformative() const {
    return _isHaplotypeInformative;
}

int Read::getEndIndex() const {
    return end;
}

int Read::getStartIndex() const {
    return start;
}

int Read::distance(std::string &partialHaplotype, int pos) const {
    int dist = 0;
    for (int i = 0;i < partialHaplotype.size(); i++) {
        if (seq[pos + i] != '-' && partialHaplotype[i] != seq[pos + i]) {
            dist ++;
        }
    }
    return dist;
}

int Read::minDistance(std::string &partialHaplotype, int pos) const {
    std::string compHap = partialHaplotype;
    for (int i = 0;i< compHap.length();i++) {
        if (compHap[i] == '1') {
            compHap[i] = '0';
        } else if (compHap[i] == '0') {
            compHap[i] = '1';
        }
    }
    return std::min(distance(partialHaplotype, pos), distance(compHap, pos));
}

int Read::getEndIndex(std::vector<Read *> &reads) {
    int endIndex = 0;
    for (int i = 0;i < reads.size(); i++) {
        if (endIndex < reads[i]->getEndIndex()) {
            endIndex = reads[i]->getEndIndex();
        }
    }
    return endIndex;
}
}

