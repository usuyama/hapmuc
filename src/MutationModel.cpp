/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include "math.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/foreach.hpp>
#include "MutationModel.hpp"
#include "log.h"
#include "Haplotype.hpp"

MutationModel::MutationModel(const std::vector<Haplotype> &haps_,
                             const std::vector<std::vector<double> > &tumorLiks_,
                             const std::vector<std::vector<double> > &normalLiks_,
                             const Parameters::BayesEMParameters &params_,
                             const std::string estimationType_) : haps(haps_), tumorLiks(tumorLiks_), normalLiks(normalLiks_), params(params_), estimationType(estimationType_) {
    if (tumorLiks.size() == 0) {
        throw std::string("tumorLiks.size() == 0");
    }
    if (normalLiks.size() == 0) {
        throw std::string("normalLiks.size() == 0");
    }
    if (tumorLiks[0].size() != 4) {
        throw std::string("tumorLiks[0].size() != 4");
    }
    if (normalLiks[0].size() != 4) {
        throw std::string("normalLiks[0].size() != 4");
    }
    rlt.reserve(tumorLiks.size() * 4);
    for (int i = 0; i < tumorLiks.size(); i++) {
        rlt.insert(rlt.end(), tumorLiks[i].begin(), tumorLiks[i].end());
    }
    rln.reserve(normalLiks.size() * 4);
    for (int i = 0; i < normalLiks.size(); i++) {
        rln.insert(rln.end(), normalLiks[i].begin(), normalLiks[i].end());
    }
}

inline double MutationModel::sumKthCol(std::vector<double> ar, int k) {
    double sum = 0.0;
    int nr = (int)ar.size() / 4;
    for (int i = 0; i < nr; i++) {
        sum += ar[i * 4 + k];
    }
    return sum;
}

inline std::vector<double> MutationModel::sumEachCol(std::vector<double> ar) {
    std::vector<double> ans(4, 0.0);
    for (int h = 0; h < 4; h++) {
        ans[h] = sumKthCol(ar, h);
    }
    return ans;
}

inline double MutationModel::sumVec(std::vector<double> ar) {
    double sum = 0.0;
    for (int i = 0; i < ar.size(); i++) {
        sum += ar[i];
    }
    return sum;
}

//l_e_pit <- digamma(new_al) - digammma(sum(new_al))のための関数
inline std::vector<double> MutationModel::calDirExp(std::vector<double> ar) {
    double sum = sumVec(ar);
    std::vector<double> ans(ar.size(), 0.0);
    try {
        double dig_sum = boost::math::digamma(sum);
        for (int i = 0; i < ar.size(); i++) {
            ans[i] = boost::math::digamma(ar[i]) - dig_sum;
        }
    } catch (std::exception &e) {
        std::string message = std::string("error_exception_").append(e.what());
        LOG(logDEBUG) << message << std::endl;
        LOG(logDEBUG) << "calDirExp" << std::endl;
        BOOST_FOREACH(double x, ar) {
            LOG(logDEBUG) << x << " ";
        }
        LOG(logDEBUG) << std::endl;
        LOG(logDEBUG) << "sum: " << sum << std::endl;
        LOG(logDEBUG).flush();
        throw;
    }
    return ans;
}

inline std::vector<double> MutationModel::normVec(std::vector<double> ar) {
    double sum = sumVec(ar);
    std::vector<double> ans(ar.size(), 0.0);
    for (int i = 0; i < ar.size(); i++) {
        ans[i] = ar[i] / sum;
    }
    return ans;
}

inline std::vector<double> MutationModel::expVec(std::vector<double> ar) {
    std::vector<double> ans(ar.size(), 0.0);
    for (int i = 0; i < ar.size(); i++) {
        ans[i] = exp(ar[i]);
    }
    return ans;
}


inline void MutationModel::printVarInHap(const Haplotype &hap) {
    for (std::vector<Variant>::const_iterator it = hap.variants.begin(); it != hap.variants.end(); it++) {
        LOGP(logDEBUG) << it->originalString << std::endl;
    }
    LOG(logDEBUG) << std::endl;
}

template <class X> void printVec(std::vector<X> vec) {
    BOOST_FOREACH(X x, vec) {
        LOGP(logDEBUG) << x << " ";
    }
}

inline double MutationModel::calELogDir(std::vector<double> log_theta, std::vector<double> param) {
    //E[log p(theta|param)] : Dir
    //- sum(lgamma(param)) + lgamma(sum(param)) + sum((param-1) * log_theta)
    double sum_lgam_param = 0.0;
    for (int i = 0; i < param.size(); i++) {
        sum_lgam_param -= lgamma(param[i]);
    }
    double lgam_sum_param = lgamma(sumVec(param));
    double sum_par_lt = 0.0;
    for (int i = 0; i < param.size(); i++) {
        sum_par_lt += (param[i] - 1) * log_theta[i];
    }
    if (log_theta.size() != param.size()) {
        LOG(logDEBUG) << "--- calELogDir ---" << std::endl;
        printVec(log_theta);
        LOG(logDEBUG) << std::endl;
        printVec(param);
        LOG(logDEBUG) << std::endl;
        LOG(logDEBUG) << sum_lgam_param << " " << lgam_sum_param << " " << sum_par_lt << std::endl;
        throw std::string("calELogDir");
    }
    return sum_lgam_param + lgam_sum_param + sum_par_lt;
}

inline double MutationModel::computeLowerBound(const std::vector<double> &rlt, const std::vector<double> &rln, const std::vector<double> &ahat, const std::vector<double> &bhat, const std::vector<double> &chat, const std::vector<double> &zt, const std::vector<double> &zn, Parameters::BayesEMParameters params) {
    int nh = 4; //number of haps
    int nrt = (int)rlt.size() / nh; //number of reads of tumor
    int nrn = (int)rln.size() / nh; //number of reads of normal
    std::vector<double> l_e_pit = calDirExp(ahat);
    std::vector<double> l_e_ep = calDirExp(bhat);
    std::vector<double> l_e_pin = calDirExp(chat);
    double e_log_p_Xt = 0.0;
    for (int i = 0; i < rlt.size(); i++) {
        e_log_p_Xt += rlt[i] * zt[i];
    }
    double e_log_p_Xn = 0.0;
    for (int i = 0; i < rln.size(); i++) {
        e_log_p_Xn += rln[i] * zn[i];
    }
    double e_log_p_Zt = 0.0, e_log_p_Zn = 0.0;
    std::vector<double> log_rho_t(4, 0), log_rho_n(4, 0);
    if (ahat.size() == 3) {
        log_rho_t[0] = l_e_pit[0];
        log_rho_t[1] = l_e_pit[1] + l_e_ep[1];
        log_rho_t[2] = l_e_pit[2];
        log_rho_t[3] = l_e_pit[1] + l_e_ep[0];
    } else {
        log_rho_t[0] = l_e_pit[0] + l_e_ep[1];
        log_rho_t[1] = l_e_pit[1] + l_e_ep[1];
        log_rho_t[2] = l_e_pit[0] + l_e_ep[0];
        log_rho_t[3] = l_e_pit[1] + l_e_ep[0];
    }
    for (int i = 0; i < nrt; i++) {
        for (int k = 0; k < 4; k++) {
            e_log_p_Zt += zt[i * 4 + k] * log_rho_t[k];
        }
    }
    log_rho_n[0] = l_e_pin[0] + l_e_ep[1];
    log_rho_n[1] = l_e_pin[1] + l_e_ep[1];
    log_rho_n[2] = l_e_pin[0] + l_e_ep[0];
    log_rho_n[3] = l_e_pin[1] + l_e_ep[0];
    for (int i = 0; i < nrn; i++) {
        for (int k = 0; k < 4; k++) {
            e_log_p_Zn += zn[i * 4 + k] * log_rho_n[k];
        }
    }
    double e_log_p_pit, e_log_p_pin, e_log_p_ep;
    if (ahat.size() == 2) {
        e_log_p_pit = calELogDir(l_e_pit, params.err_a0);
        e_log_p_pin = calELogDir(l_e_pin, params.err_c0);
        e_log_p_ep = calELogDir(l_e_ep, params.err_b0);
    } else {
        e_log_p_pit = calELogDir(l_e_pit, params.mut_a0);
        e_log_p_pin = calELogDir(l_e_pin, params.mut_c0);
        e_log_p_ep = calELogDir(l_e_ep, params.mut_b0);
    }
    double e_log_q_Zt = 0.0, e_log_q_Zn = 0.0;
    for (int i = 0; i < zt.size(); i++) {
        e_log_q_Zt += zt[i] * log(zt[i]);
    }
    for (int i = 0; i < zn.size(); i++) {
        e_log_q_Zn += zn[i] * log(zn[i]);
    }
    double e_log_q_pit = calELogDir(l_e_pit, ahat);
    double e_log_q_pin = calELogDir(l_e_pin, chat);
    double e_log_q_ep = calELogDir(l_e_ep, bhat);
    double lower_bound = e_log_p_Xt + e_log_p_Xn + e_log_p_Zt + e_log_p_Zn + e_log_p_pin + e_log_p_pit + e_log_p_ep - e_log_q_Zt - e_log_q_Zn - e_log_q_pit - e_log_q_pin - e_log_q_ep;
    //LOG(logDEBUG) << "lb: " << lower_bound << std::endl;
    //LOG(logDEBUG) << e_log_p_Xt << " " << e_log_p_Xn << " " << e_log_p_Zt << " " << e_log_p_Zn << " " << e_log_p_pit << " " << e_log_p_pin << " " << e_log_p_ep << std::endl;
    //LOG(logDEBUG) << e_log_q_Zt << " " << e_log_q_Zn << " " << e_log_q_pit << " " << e_log_q_pin << " " << e_log_q_ep << std::endl;
    return lower_bound;
}

/*
 void MutationModel::output_params(std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, std::vector<double> & zt, std::vector<double> & zn, std::string fname_prefix) {
 #ifdef LOGDEBUG
 std::ofstream zt_ofs((fname_prefix+".zt").c_str(), std::ios::out|std::ios::trunc);
 std::ofstream zn_ofs((fname_prefix+".zn").c_str(), std::ios::out|std::ios::trunc);
 int nrt = zt.size()/4;
 int nrn = zn.size()/4;
 for(int i=0;i<nrt;i++) {
 zt_ofs << zt[i*4];
 for(int j=1;j<4;j++) {
 zt_ofs << "\t" << zt[i*4+j];
 }
 zt_ofs << "\n";
 }
 for(int i=0;i<nrn;i++) {
 zn_ofs << zn[i*4];
 for(int j=1;j<4;j++) {
 zn_ofs << "\t" << zn[i*4+j];
 }
 zn_ofs << "\n";
 }
 #endif
 }

 void MutationModel::print_params(const std::vector<ReadPair> & tumorReadPairs, const std::vector<ReadPair> & normalReadPairs, const std::vector<double> & rlt, const std::vector<double> & rln, std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, std::vector<double> & zt, std::vector<double> & zn, Parameters::BayesEMParameters params, std::ofstream &log, bool print_z=false) {
 #ifdef LOGMM
 int nh = 4; //number of haps
 int nrt = rlt.size() / nh; //number of reads of tumor
 int nrn = rln.size() / nh; //number of reads of normal
 log << "##print params" << std::endl;
 if(print_z) {
 log << "###zt###" << std::endl;
 for (int r=0;r<nrt;r++) {
 #ifndef MMTEST
 log << tumorReadPairs[r].first->seq_name;
 #endif
 log << " zt[" << r << "]:";
 for (int h=0;h<nh;h++) log << " " << zt[r*nh+h];
 log << std::endl;
 }
 log << "###zn###" << std::endl;
 for (int r=0;r<nrn;r++) {
 #ifndef MMTEST
 log << normalReadPairs[r].first->seq_name;
 #endif
 log << " zn[" << r << "]:";
 for (int h=0;h<nh;h++) log << " " << zn[r*nh+h];
 log << std::endl;
 }
 }
 log << "ahat: ";printVec(ahat, log);
 log << endl << "bhat: ";printVec(bhat, log);
 log << endl << "chat: ";printVec(chat, log);
 log << std::endl;
 #endif
 }
*/

inline void MutationModel::printPosteriors(const std::vector<double> &rlt, const std::vector<double> &rln, std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, std::vector<double> &zt, std::vector<double> &zn) {
    //LOG(logDEBUG) << "## print posteriors ##" << std::endl;
    std::vector<double> l_e_pit = calDirExp(ahat);
    std::vector<double> l_e_ep = calDirExp(bhat);
    std::vector<double> l_e_pin = calDirExp(chat);
    LOG(logDEBUG) << "pit: ";
    printVec(expVec(l_e_pit));
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "ep: ";
    printVec(expVec(l_e_ep));
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG)  << "pin: ";
    printVec(expVec(l_e_pin));
    LOGP(logDEBUG) << std::endl;
    std::vector<double> nkt = sumEachCol(zt);
    std::vector<double> nkn = sumEachCol(zn);
    LOG(logDEBUG) << "nkt: ";
    printVec(nkt);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "nkn: ";
    printVec(nkn);
    LOGP(logDEBUG) << std::endl;
}


inline std::vector<double> MutationModel::makeNewrho(const std::vector<double> &rl, std::vector<double> &l_e_pi, std::vector<double> &l_e_ep, int b) {
    std::vector<double> new_rho(4, 0.0);
    if (l_e_pi.size() == 3) {
        new_rho[0] = rl[b + 0] + l_e_pi[0];
        new_rho[1] = rl[b + 1] + l_e_pi[1] + l_e_ep[1];
        new_rho[2] = rl[b + 2] + l_e_pi[2];
        new_rho[3] = rl[b + 3] + l_e_pi[1] + l_e_ep[0];
    } else {
        new_rho[0] = rl[b + 0] + l_e_pi[0] + l_e_ep[1];
        new_rho[1] = rl[b + 1] + l_e_pi[1] + l_e_ep[1];
        new_rho[2] = rl[b + 2] + l_e_pi[0] + l_e_ep[0];
        new_rho[3] = rl[b + 3] + l_e_pi[1] + l_e_ep[0];
    }
    return new_rho;
}

inline void MutationModel::updateHat(const std::vector<double> &nkt, const std::vector<double> &nkn, std::vector<double> &ahat, std::vector<double> &bhat, std::vector<double> &chat, Parameters::BayesEMParameters params) {
    if (ahat.size() == 3) {
        ahat[0] = nkt[0] + params.mut_a0[0];
        ahat[1] = nkt[1] + nkt[3] + params.mut_a0[1];
        ahat[2] = nkt[2] + params.mut_a0[2];

        bhat[0] = nkt[3] + nkn[2] + nkn[3] + params.mut_b0[0];
        bhat[1] = nkt[1] + nkn[0] + nkn[1] + params.mut_b0[1];

        chat[0] = nkn[0] + nkn[2] + params.mut_c0[0];
        chat[1] = nkn[1] + nkn[3] + params.mut_c0[1];
    } else {
        ahat[0] = nkt[0] + nkt[2] + params.err_a0[0];
        ahat[1] = nkt[1] + nkt[3] + params.err_a0[1];

        bhat[0] = nkt[2] + nkt[3] + nkn[2] + nkn[3] + params.err_b0[0];
        bhat[1] = nkt[0] + nkt[1] + nkn[0] + nkn[1] + params.err_b0[1];

        chat[0] = nkn[0] + nkn[2] + params.err_c0[0];
        chat[1] = nkn[1] + nkn[3] + params.err_c0[1];
    }
}

double MutationModel::estimate() {
    LOG(logINFO) << "==== MutationModel.estimate() ====" << std::endl;
    LOG(logINFO) << "ModelType: " << estimationType << std::endl;
    LOG(logDEBUG) << "mut a0: ";
    printVec(params.mut_a0);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "mut b0: ";
    printVec(params.mut_b0);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "mut c0: ";
    printVec(params.mut_c0);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "err a0: ";
    printVec(params.err_a0);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "err b0: ";
    printVec(params.err_b0);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "err c0: ";
    printVec(params.err_c0);
    LOGP(logDEBUG) << std::endl;
#ifndef MMTEST
    if (haps.size() != 4) {
        throw std::string("error: only for 4 haps");
    }
#endif
    size_t nh = 4;
    size_t nrt = rlt.size() / 4;
    size_t nrn = rln.size() / 4;

    std::vector<double> zt(nh * nrt, 0.0); // expectations of read-haplotype indicator variables
    std::vector<double> zn(nh * nrn, 0.0); // expectations of read-haplotype indicator variables
    std::vector<double> nkt(nh, 0.0), nkn(nh, 0.0);

#ifndef MMTEST
    LOG(logDEBUG) << "#haplotype list" << std::endl;
    for (size_t th = 0; th < nh; th++) {
        const Haplotype &hap = haps[th];
        LOG(logDEBUG) << "hap[" << th << "] ";
        for (std::vector<Variant>::const_iterator it = hap.variants.begin(); it != hap.variants.end(); it++) {
            LOGP(logDEBUG) << "[" <<  it->originalString << " " << it->startInGenome << "]";
        }
        LOGP(logDEBUG) << std::endl;
    }
#endif
    // check haplotypes
    // run EM for this set of active variants
    bool converged = false;
    double tol = 1e-5;
    // initialize frequencies
    std::vector<double> ahat, bhat, chat;
    if (estimationType == "mutation") {
        ahat = params.mut_a0;
        bhat = params.mut_b0;
        chat = params.mut_c0;
    } else {
        ahat = params.err_a0;
        bhat = params.err_b0;
        chat = params.err_c0;
    }
    LOG(logDEBUG) << "init:" << std::endl;
    LOG(logDEBUG) << "ahat:";
    printVec(ahat);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "bhat:";
    printVec(bhat);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "chat:";
    printVec(chat);
    LOGP(logDEBUG) << std::endl;
    double llNew, llOld = -HUGE_VAL;
    int iter = 0;
    while (!converged) {
        std::vector<double> l_e_pit = calDirExp(ahat);
        std::vector<double> l_e_ep = calDirExp(bhat);
        for (int i = 0; i < nrt; i++) {
            std::vector<double> new_rho = makeNewrho(rlt, l_e_pit, l_e_ep, i * 4);
            std::vector<double> zt_i = normVec(expVec(new_rho));
            for (int k = 0; k < 4; k++) {
                zt[i * 4 + k] = zt_i[k];
            }
        }
        std::vector<double> l_e_pin = calDirExp(chat);
        for (int i = 0; i < nrn; i++) {
            std::vector<double> new_rho = makeNewrho(rln, l_e_pin, l_e_ep, i * 4);
            std::vector<double> zn_i = normVec(expVec(new_rho));
            for (int k = 0; k < 4; k++) {
                zn[i * 4 + k] = zn_i[k];
            }
        }
        nkt = sumEachCol(zt);
        nkn = sumEachCol(zn);
        updateHat(nkt, nkn, ahat, bhat, chat, params);
        llNew = computeLowerBound(rlt, rln, ahat, bhat, chat, zt, zn, params);
        converged = (fabs(llOld - llNew)) < tol || iter > 500;
        printPosteriors(rlt, rln, ahat, bhat, chat, zt, zn);
        LOG(logDEBUG) << "==== iter: " << iter << " llOld: " << llOld << " llNew: " << llNew << " ====" << std::endl;
        //print_params(tumorReadPairs, normalReadPairs, rlt, rln, ahat, bhat, chat, zt, zn, params, true);
        //LOG(logDEBUG) << std::endl;
        llOld = llNew;
        iter++;
    }
    double lower_bound = llNew;
    LOG(logINFO) << "---- finished[" << iter << "] ----" << std::endl;
    //print_params(tumorReadPairs, normalReadPairs, rlt, rln, ahat, bhat, chat, zt, zn, params, log, true);
    //output_params(ahat, bhat, chat, zt, zn, log_prefix);

    LOG(logDEBUG) << "estimated parameters" << std::endl;
    LOG(logDEBUG) << "ahat:";
    printVec(ahat);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "bhat:";
    printVec(bhat);
    LOGP(logDEBUG) << std::endl;
    LOG(logDEBUG) << "chat:";
    printVec(chat);
    
    std::vector<double> l_e_pin = calDirExp(chat);
    if (!(l_e_pin[0] > log(0.15) && l_e_pin[0] < log(0.85))) {
        throw std::string("umbalanced_normal_haplotypes");
    }

    LOG(logINFO) << "lower bound = " << lower_bound << std::endl;
    return lower_bound;
}

