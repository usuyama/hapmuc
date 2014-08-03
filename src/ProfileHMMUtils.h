#ifndef __haplotype_builder__ProfileHMMUtils__
#define __haplotype_builder__ProfileHMMUtils__

#include "float.h"

#include <typeinfo>
#include <vector>

#include "Alignment.h"
#include "PairedRead.h"
#include "SingleRead.h"
#include "Haplotype.hpp"
#include "ProfileHMM.h"
#include "Parameters.h"

class ProfileHMMUtils {
public:
    static inline double calc(const std::vector<Variant> &allVariants,
                              const Haplotype &hap, const PairedRead &read) {
        double output = 0.0;
        bool covered1 = false, covered2 = false;
        for (int i = 0; i < allVariants.size(); i++) {
            if (read.first->cover(allVariants[i].startInGenome)) {
                covered1 = true;
            }
            if (read.second != NULL && read.second->cover(allVariants[i].startInGenome)) {
                covered2 = true;
            }
            if (covered1 && covered2) {
                break;
            }
        }
        if (covered1) {
            output += calc(hap, read.first);
        }
        if (covered2) {
            try {
                output += calc(hap, read.second);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "mate-pair mapping fail; continue;" << std::endl;
            }
        }
        return output;
    }

    static inline double calc(const std::vector<Variant> &allVariants,
                              const Haplotype &hap, const SingleRead &read) {
        bool covered = false;
        for (int i = 0; i < allVariants.size(); i++) {
            if (read.read->cover(allVariants[i].startInGenome)) {
                covered = true;
                break;
            }
        }
        if (covered) {
            return calc(hap, read.read);
        } else {
            return 0.0;
        }
    }
    
    static inline double calc(const std::vector<Variant> &allVariants,
                              const Haplotype &hap, const IRead *read) {
        if (typeid(*read) == typeid(PairedRead)) {
            return calc(allVariants, hap, *dynamic_cast<const PairedRead *>(read));
        } else if (typeid(*read) == typeid(SingleRead)) {
            return calc(allVariants, hap, *dynamic_cast<const SingleRead *>(read));
        } else {
            throw std::string("ProfileHMMUtils::calc not_implemented");
        }
    }
    
    static inline double calc(const Haplotype &hap, const SamRead *sam) {
        ProfileHMM pHMM0 = genProfileHMM(hap, sam);
        Alignment a0 = pHMM0.viterbi();
        return a0.likelihood;
        /* TODO
        int th1 = -100, th2 = -200;
        if (a0.likelihood > th1) {
           // a0.print();
            return a0.likelihood;
        } else {
            if (sam->hasSoftClip) {
                std::string newSeq;
                std::vector<double> newBQ;
                int newPos;
                sam->getReadWithoutSoftClipping(&newSeq, &newBQ, &newPos);
                ProfileHMM pHMM1 = ProfileHMM(hap, newSeq, newBQ, newPos,
                                              sam->mapQual,
                                              true,
                                              sam->onReverseStrand, sam->seq_name);
                Alignment a1 = pHMM1.viterbi();
                if (a1.likelihood > th2) {
                    if (a1.likelihood < th1) {
                        a1.print();
                    }
                    LOG(logDEBUG) << "aligned without soft-clip (" << sam->seq_name << "): " << a1.likelihood << std::endl;
                    return a1.likelihood;
                } else {
                    a0.print();
                    a1.print();
                    throw std::string("did_not_map_profileHMM2");
                }
            } else if (a0.likelihood > th2) {
                return a0.likelihood;
            } else {
                a0.print();
                throw std::string("did_not_map_profileHMM1");
            }
        }
         */
    }
    
    static std::vector<std::vector<double> > calcAllLikelihoods(const std::vector<Haplotype> &haps,
                                                                const std::vector<IRead *> &reads) {
        std::vector<std::vector<double> > liks;
        liks.reserve(reads.size());
        std::vector<Variant> allVariants = Haplotype::getAllVariants(haps);
    
        for (int i = 0; i < reads.size(); i++) {
            try {
                std::vector<double> tmpLiks(haps.size());
                for (int j = 0; j < haps.size(); j++) {
                    tmpLiks[j] = calc(allVariants, haps[j], reads[i]);
                }
                liks.push_back(tmpLiks);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }
    
    static inline ProfileHMM genProfileHMM(const Haplotype &hap, const SamRead *sam) {
        return ProfileHMM(hap, sam->seq, sam->qual, sam->leftMostPos, sam->mapQual, sam->hasIndel || sam->hasSoftClip || sam->hasHardClip, sam->onReverseStrand, sam->seq_name);
    }
    
    static std::vector<std::vector<double> > calcLikelihoods(const std::vector<Haplotype> &haps,
                                                             const std::vector<IRead *> &reads) {
        std::vector<std::vector<double> > liks = calcAllLikelihoods(haps, reads);
        return liks;
    }

    static std::vector<std::vector<double> > calcBasicLikelihoods(const std::vector<Haplotype> &haps,
                                                                  const std::vector<IRead *> &reads) {
        std::vector<std::vector<double> > liks;
        liks.reserve(reads.size());

        for (int i = 0; i < reads.size(); i++) {
            try {
                std::vector<double> tmpLiks(haps.size());
                for (int j = 0; j < haps.size(); j++) {
                    tmpLiks[j] = calc(haps[j], reads[i]->getPrimarySamRead());
                }
                liks.push_back(tmpLiks);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }

};



#endif /* defined(__haplotype_builder__ProfileHMMUtils__) */
