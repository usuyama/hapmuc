/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
/*
 * Dindel
 * http://www.sanger.ac.uk/resources/software/dindel/
 *
 * Copyright 2010, Kees Albers
 * Released under the GPL Version 3.
 */
#ifndef HAPS_CPP_
#define HAPS_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include "bam.h"
#include "Haplotype.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
#include "Utils.hpp"
#include "faidx.h"
#include "ObservationModelSeqAn.hpp"
#include "VariantFile.hpp"
#include <ext/hash_map>
#include <exception>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "MutationCall.hpp"
#include "Haps.hpp"
#include "EMBasic.hpp"
#include "log.h"

using namespace seqan;

namespace Haps {
    
    void computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap, Parameters params)
    {
        LOG(logDEBUG) << "### Computing likelihoods for all reads and haplotypes.\n";
        onHap = vector<int>(reads.size(),0); // records whether a read was aligned onto at least one haplotype
        
        typedef map<size_t, vector<size_t> >::const_iterator hapsCIt;
        
        liks=vector<vector<MLAlignment> >(haps.size(),vector<MLAlignment>(reads.size()));
        for (size_t r=0;r<reads.size();r++) {
            for (size_t hidx=0;hidx<haps.size();hidx++) {
                const Haplotype & hap=haps[hidx];
                ObservationModelFBMaxErr oms(hap, reads[r], leftPos, params.obsParams);
                liks[hidx][r]=oms.calcLikelihood();
                if (!liks[hidx][r].offHapHMQ) onHap[r]=1;
                /*
                 LOG(logDEBUG) << "---" << endl;
                 LOG(logDEBUG) <<  "read: " << bam1_qname(reads[r].getBam()) << ", hidx: " << hidx << " mpos: " << reads[r].matePos << endl;
                 LOG(logDEBUG) << "isUnmapped: " << reads[r].isUnmapped() << endl;
                 LOG(logDEBUG) << string(50,' ') << haps[hidx].seq << endl;
                 oms.printAlignment(50);*/
                if (liks[hidx][r].ll>0.1) {
                    LOG(logDEBUG) << "warning" << endl;
                    ObservationModelFBMaxErr om(hap, reads[r], leftPos, params.obsParams);
                    liks[hidx][r]=om.calcLikelihood();
                    LOG(logDEBUG) << string(25,' ') << hap.seq << endl;
                    om.printAlignment(25);
                    LOG(logDEBUG) << "hidx: " << hidx << " r: " << r << endl;
                    LOG(logDEBUG) << bam1_qname(reads[r].getBam()) << endl;
                    cerr << "Likelihood>0" << endl;
                    exit(1);
                }
                if (isnan(liks[hidx][r]) || isinf(liks[hidx][r])) {
                    LOG(logDEBUG) << "NAN/Inf error" << endl;
                    throw string("Nan detected");
                }
            }
        }
        LOG(logDEBUG) << "computeLikelihoods done" << endl;
    }

}
#endif