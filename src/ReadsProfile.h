#ifndef __haplotype_builder__ReadsProfile__
#define __haplotype_builder__ReadsProfile__

#include <vector>

#include "IRead.h"
#include "Variant.h"

class ReadsProfile {
public:
    class Profile {
    public:
        int numReads;
        double averageMapQuality;
        double softClipPosessionFreq;
        double indelPosesssionFreq;
        friend std::ostream& operator<<(std::ostream &os, const Profile &r) {
            os << "# reads = " << r.numReads;
            os << ", avg. map quality: " << r.averageMapQuality;
            os << ", soft-clipping posession frequency: " << r.softClipPosessionFreq;
            os << ", indel posession frequency: " << r.indelPosesssionFreq;
            return os;
        }
    };
    Profile withVariant, withoutVariant, all;
    Profile makeData(const std::vector<IRead *> &reads);
    ReadsProfile(const std::vector<IRead *> &reads, const Variant &variant);
    friend std::ostream& operator<<(std::ostream &os, const ReadsProfile &r) {
        os << "all:  " << r.all << std::endl;
        os << "variant-supporting: " << r.withVariant << std::endl;
        os << "not variant-supporting: " << r.withoutVariant;
        return os;
    }
};

#endif /* defined(__haplotype_builder__ReadsProfile__) */
