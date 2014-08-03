#ifndef ReadTest_Parameters_h
#define ReadTest_Parameters_h

#include <cstdlib>
#include <string>
#include <vector>
#include "log.h"
#include "Utils.h"

class Parameters {
public:
    class BayesEMParameters { // for Bayesian EM algorithm
    public:
        std::vector<double> mut_a0, mut_b0, mut_c0;
        std::vector<double> err_a0, err_b0, err_c0;
        BayesEMParameters() {};
        ~BayesEMParameters() {};
        void clear_all() {
            mut_a0.clear();
            mut_b0.clear();
            mut_c0.clear();
            err_a0.clear();
            err_b0.clear();
            err_c0.clear();
        }

        template <class X> void print_veci(std::vector<X> vec) {
            for (int i = 0; i < vec.size(); i++) {
                LOGP(logINFO) << vec[i] << " ";
            }
        }

        void print() {
            LOG(logINFO) << "\tMutationModel:" << std::endl;
            LOG(logINFO) << "\t\ta0: ";
            print_veci(mut_a0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tb0: ";
            print_veci(mut_b0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tc0: ";
            print_veci(mut_c0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\tErrorModel:" << std::endl;
            LOG(logINFO) << "\t\ta0: ";
            print_veci(err_a0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tb0: ";
            print_veci(err_b0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tc0: ";
            print_veci(err_c0);
            LOGP(logINFO) << std::endl;
        }

        BayesEMParameters &operator=(const BayesEMParameters &t) {
            mut_a0 = t.mut_a0;
            mut_b0 = t.mut_b0;
            mut_c0 = t.mut_c0;
            err_a0 = t.err_a0;
            err_b0 = t.err_b0;
            err_c0 = t.err_c0;
            return *this;
        }
    };

    Parameters();
    Parameters(int argc, const char *argv[]);
    void setDefaultValues();
    void getFromCommandLineArguments(int argc, const char *argv[]);
    static void parseHyperParameters(std::vector<double> &vec, std::string s);
    int maxReads, maxReadLength, minReads;
    int maxInsertSize, maxIndelLength;
    int minDistanceSNP, minDistanceGermlineIndel;
    double priorDelByError; //baseQualThreshold, 
    double averageMapQualThreshold, mapQualThreshold, priorIndel, priorSNP, EMtol;
    double pError, pMut;
    std::string outFilePrefix, refFileName, tumorBam, normalBam, windowFile;
//    std::string filterReadAux;
    bool showHapAlignments, quiet;
    bool singleReads;
    BayesEMParameters hap3_params, hap2_params;
    double indelPosessionFreqThreshold, softClipPosessionFreqThreshold;
    bool withoutBayesFactor;
};
#endif
