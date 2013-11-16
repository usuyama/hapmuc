#ifndef HAPS_HPP_
#define HAPS_HPP_

#include <vector>
#include "HapMuC.hpp"
#include "Read.hpp"
#include "Haplotype.hpp"
#include "MLAlignment.hpp"

namespace Haps {
    
    void computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap, Parameters params);
    
}
#endif
