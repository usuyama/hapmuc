/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef MM_HPP_
#define MM_HPP_

#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include "Variant.h"
#include "Parameters.h"
#include "Haplotype.hpp"

class MutationModel {
    static inline double sumVec(std::vector<double> ar);
    static inline std::vector<double> calDirExp(std::vector<double> ar);
    static inline std::vector<double> normVec(std::vector<double> ar);
    static inline std::vector<double> expVec(std::vector<double> ar);
    static inline std::vector<double> sumEachCol(std::vector<double> ar);
    static inline double sumKthCol(std::vector<double> ar, int k);
    static void printVarInHap(const Haplotype &hap);
    static inline double calELogDir(std::vector<double> logTheta, std::vector<double> param);
    static inline double computeLowerBound(const std::vector<double> &rlt, const std::vector<double> &rln, const std::vector<double> &ahat, const std::vector<double> &bhat, const std::vector<double> &chat, const std::vector<double> &zt, const std::vector<double> &zn, Parameters::BayesEMParameters params);
    // void output_params(std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, std::vector<double> & zt, std::vector<double> & zn, std::string fname_prefix);
    // void print_params(const std::vector<ReadPair> & tumorReadPairs, const std::vector<ReadPair> & normalReadPairs, const std::vector<double> & rlt, const std::vector<double> & rln, std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, std::vector<double> & zt, std::vector<double> & zn, Parameters::BayesEMParameters params, std::ofstream &log, bool print_z);
    void printPosteriors(const std::vector<double> &rlt, const std::vector<double> &rln, std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, std::vector<double> &zt, std::vector<double> &zn);
    inline std::vector<double> makeNewrho(const std::vector<double> &rl, std::vector<double> &l_e_pi, std::vector<double> &l_e_ep, int b);
    inline void updateHat(const std::vector<double> &nkt, const std::vector<double> &nkn, std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, Parameters::BayesEMParameters params);
public:
    const std::vector<Haplotype> &haps;
    const std::vector<std::vector<double> > &tumorLiks, &normalLiks;
    std::vector<double> rlt, rln; // 4 * (# of reads); index(i, h) = i * 4 + h
    //i番目のリード, k番目のhap;i行k列 -> i*4+k
    const Parameters::BayesEMParameters &params;
    std::string estimationType;
    MutationModel(const std::vector<Haplotype> &haps, const std::vector<std::vector<double> > &tumorLiks, const std::vector<std::vector<double> > &normalLiks, const Parameters::BayesEMParameters &params, std::string estimationType);  //est_type=mutation or non-mutation
    typedef std::pair<int, Variant> PAV;
    typedef std::map<int, Variant>::const_iterator It;
    typedef std::map<int, std::set<PAV> >::const_iterator PIt;
    double estimate();
};

#endif
