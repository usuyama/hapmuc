#include "History.h"

#include <algorithm>

namespace haplotype_phaser {

struct CompareSecond
{
    bool operator()(const std::pair<std::string, int>& left, const std::pair<std::string, int>& right) const
    {
        return left.second < right.second;
    }
};

std::pair<std::string, int> History::getMinimumState() {
    std::pair<std::string, int> min = *std::min_element(scores.begin(), scores.end(), CompareSecond());
    return min;
}

}
