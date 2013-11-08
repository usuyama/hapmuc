#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include "foreach.hpp"
#include "bam.h"
#include "HapMuC.hpp"
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
#include <boost/algorithm/string.hpp>
#include "MutationModel.hpp"
#include "log.h"

namespace MutationModel {
    typedef pair<int, AlignedVariant> PAV;
    typedef map<int, AlignedVariant>::const_iterator It;
    typedef map<int, std::set<PAV> >::const_iterator PIt;

    double sum_kth_col(vector<double> ar, int k) {
        //i番目のリード, k番目のhap;i行k列 -> i*4+k
        double sum=0.0;
        int nr=ar.size()/4;
        for(int i=0;i<nr;i++) sum+=ar[i*4+k];
        return sum;
    }

    vector<double> sum_each_col(vector<double> ar) {
        vector<double> ans(4, 0.0);
        for(int h=0;h<4;h++) ans[h] = sum_kth_col(ar, h);
        return ans;
    }

    double sum_vec(vector<double> ar) {
        double sum = 0.0;
        for(int i = 0;i < ar.size();i++) sum += ar[i];
        return sum;
    }

    //l_e_pit <- digamma(new_al) - digammma(sum(new_al))のための関数
    vector<double> cal_dir_exp(vector<double> ar) {
        double sum = sum_vec(ar);
        vector<double> ans(ar.size(), 0.0);
        try {
            double dig_sum = boost::math::digamma(sum);
            for(int i=0;i< ar.size();i++) {
                ans[i]=boost::math::digamma(ar[i]) - dig_sum;
            }
        } catch (std::exception& e) {
            string message = string("error_exception_").append(e.what());
            LOG(logDEBUG) << message << endl;
            LOG(logDEBUG) << "cal_dir_exp" << endl;
            BOOST_FOREACH(double x, ar) { LOG(logDEBUG) << x << " "; }
            LOG(logDEBUG) << endl;
            LOG(logDEBUG) << "sum: " << sum << endl;
            LOG(logDEBUG).flush();
            throw;
        }
        return ans;
    }

    vector<double> norm_vec(vector<double> ar) {
        double sum = sum_vec(ar);
        vector<double> ans(ar.size(), 0.0);
        for(int i=0;i<ar.size();i++) {
            ans[i] = ar[i] / sum;
        }
        return ans;
    }

    vector<double> exp_vec(vector<double> ar) {
        vector<double> ans(ar.size(), 0.0);
        for(int i = 0;i<ar.size();i++) {
            ans[i] = exp(ar[i]);
        }
        return ans;
    }

    void print_var_in_hap(Haplotype hap, std::ofstream &log) {
#ifndef LOGDEBUG
        return;
#endif
        for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
            if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                log << "[" << it->second.getString() << " " << (it->first) << "]";
            }
        }
        log << endl;
    }

    template <class X> void print_vec(vector<X> vec, std::ofstream &log) {
#ifdef LOGDEBUG
        BOOST_FOREACH(X x, vec) { log << x << " "; }
#endif
    }

    double cal_e_log_dir(vector<double> log_theta, vector<double> param, std::ofstream &log) {
        //E[log p(theta|param)] : Dir
        //- sum(lgamma(param)) + lgamma(sum(param)) + sum((param-1) * log_theta)
        double sum_lgam_param=0.0;
        for(int i=0;i<param.size();i++) sum_lgam_param-=lgamma(param[i]);
        double lgam_sum_param=lgamma(sum_vec(param));
        double sum_par_lt=0.0;
        for(int i=0;i<param.size();i++) sum_par_lt+=(param[i]-1)*log_theta[i];
        if(log_theta.size() != param.size()) {
            log << "--- cal_e_log_dir ---" << endl;
            print_vec(log_theta, log);log << endl;
            print_vec(param, log);log << endl;
            log << sum_lgam_param << " " << lgam_sum_param << " " << sum_par_lt << endl;
            throw string("cal_e_log_dir");
        }
        return sum_lgam_param + lgam_sum_param + sum_par_lt;
    }

    void compute_lower_bound(const vector<double> & rlt, const vector<double> & rln, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, Hap3Param params, std::ofstream &logofs) {
        int nh = 4; //number of haps
        int nrt = rlt.size() / nh; //number of reads of tumor
        int nrn = rln.size() / nh; //number of reads of normal
        vector<double> l_e_pit = cal_dir_exp(ahat);
        vector<double> l_e_ep = cal_dir_exp(bhat);
        vector<double> l_e_pin = cal_dir_exp(chat);
        double e_log_p_Xt = 0.0;
        for(int i=0;i<rlt.size();i++) e_log_p_Xt += rlt[i] * zt[i];
        double e_log_p_Xn = 0.0;
        for(int i=0;i<rln.size();i++) e_log_p_Xn += rln[i] * zn[i];
        double e_log_p_Zt = 0.0, e_log_p_Zn = 0.0;
        vector<double> log_rho_t(4,0), log_rho_n(4,0);
        if (ahat.size()==3) {
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
        for(int i=0;i<nrt;i++) {
            for(int k=0;k<4;k++) e_log_p_Zt += zt[i*4+k] * log_rho_t[k];
        }
        log_rho_n[0] = l_e_pin[0] + l_e_ep[1];
        log_rho_n[1] = l_e_pin[1] + l_e_ep[1];
        log_rho_n[2] = l_e_pin[0] + l_e_ep[0];
        log_rho_n[3] = l_e_pin[1] + l_e_ep[0];
        for(int i=0;i<nrn;i++) {
            for(int k=0;k<4;k++) e_log_p_Zn += zn[i*4+k] * log_rho_n[k];
        }
        double e_log_p_pit, e_log_p_pin, e_log_p_ep;
        if(ahat.size()==2) {
            e_log_p_pit = cal_e_log_dir(l_e_pit, params.err_a0, logofs);
            e_log_p_pin = cal_e_log_dir(l_e_pin, params.err_c0, logofs);
            e_log_p_ep=cal_e_log_dir(l_e_ep, params.err_b0, logofs);
        } else {
            e_log_p_pit = cal_e_log_dir(l_e_pit, params.mut_a0, logofs);
            e_log_p_pin = cal_e_log_dir(l_e_pin, params.mut_c0, logofs);
            e_log_p_ep=cal_e_log_dir(l_e_ep, params.mut_b0, logofs);
        }
        double e_log_q_Zt=0.0, e_log_q_Zn=0.0;
        for(int i=0;i<zt.size();i++) e_log_q_Zt+=zt[i]*log(zt[i]);
        for(int i=0;i<zn.size();i++) e_log_q_Zn+=zn[i]*log(zn[i]);
        double e_log_q_pit=cal_e_log_dir(l_e_pit, ahat, logofs);
        double e_log_q_pin=cal_e_log_dir(l_e_pin, chat, logofs);
        double e_log_q_ep=cal_e_log_dir(l_e_ep, bhat, logofs);
        lb.lower_bound = e_log_p_Xt + e_log_p_Xn + e_log_p_Zt + e_log_p_Zn + e_log_p_pin + e_log_p_pit + e_log_p_ep - e_log_q_Zt - e_log_q_Zn - e_log_q_pit - e_log_q_pin - e_log_q_ep;
#ifdef LOGDEBUG
        logofs << "lb: " << lb.lower_bound << endl;
        logofs << e_log_p_Xt<< " " <<e_log_p_Xn<< " " <<e_log_p_Zt<< " " <<e_log_p_Zn<< " " <<e_log_p_pit<< " " <<e_log_p_pin << " " << e_log_p_ep << endl;
        logofs  <<e_log_q_Zt<< " " <<e_log_q_Zn<< " " <<e_log_q_pit<< " " <<e_log_q_pin<< " " <<e_log_q_ep << endl;
#endif

    }

    void output_params(vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, string fname_prefix) {
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


    void print_params(const vector<Read> & tumor_reads, const vector<Read> & normal_reads, const vector<double> & rlt, const vector<double> & rln, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, Hap3Param params, std::ofstream &log, bool print_z=false) {
#ifdef LOGDEBUG
        int nh = 4; //number of haps
        int nrt = rlt.size() / nh; //number of reads of tumor
        int nrn = rln.size() / nh; //number of reads of normal
        log << "##print params" << endl;
        if(print_z) {
            log << "###zt###" << endl;
            for (int r=0;r<nrt;r++) {
#ifndef MMTEST
                log << tumor_reads[r].seq_name;
#endif
                log << " zt[" << r << "]:";
                for (int h=0;h<nh;h++) log << " " << zt[r*nh+h];
                log << endl;
            }
            log << "###zn###" << endl;
            for (int r=0;r<nrn;r++) {
#ifndef MMTEST
                log << normal_reads[r].seq_name;
#endif
                log << " zn[" << r << "]:";
                for (int h=0;h<nh;h++) log << " " << zn[r*nh+h];
                log << endl;
            }
        }
        log << "ahat: ";print_vec(ahat, log);
        log << endl << "bhat: ";print_vec(bhat, log);
        log << endl << "chat: ";print_vec(chat, log);
        log << endl;
#endif
    }

    void print_posteriors(const vector<double> & rlt, const vector<double> & rln, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, Hap3Param params, std::ofstream &log) {
#ifdef LOGDEBUG
        log << "##print posteriors" << endl;
        int nh = 4; //number of haps
        int nrt = rlt.size() / nh; //number of reads of tumor
        int nrn = rln.size() / nh; //number of reads of normal
        vector<double> l_e_pit = cal_dir_exp(ahat);
        vector<double> l_e_ep = cal_dir_exp(bhat);
        vector<double> l_e_pin = cal_dir_exp(chat);
        log << "pit: ";print_vec(exp_vec(l_e_pit), log);
        log << endl << "ep: ";print_vec(exp_vec(l_e_ep), log);
        log << endl << "pin: ";print_vec(exp_vec(l_e_pin), log);
        vector<double> nkt = sum_each_col(zt);
        vector<double> nkn = sum_each_col(zn);
        log << endl << "nkt: ";print_vec(nkt, log);
        log << endl << "nkn: ";print_vec(nkn, log);
        log << endl;
#endif
    }

    vector<double> copyVecFast(const vector<double>& original)
    {
        vector<double> newVec;
        newVec.reserve(original.size());
        copy(original.begin(),original.end(),back_inserter(newVec));
        return newVec;
    }

     vector<double> make_newrho(const vector<double> &rl, vector<double> &l_e_pi, vector<double> &l_e_ep, int b) {
         vector<double> new_rho(4, 0.0);
         if (l_e_pi.size()==3) {
             new_rho[0] = rl[b+0] + l_e_pi[0];
             new_rho[1] = rl[b+1] + l_e_pi[1] + l_e_ep[1];
             new_rho[2] = rl[b+2] + l_e_pi[2];
             new_rho[3] = rl[b+3] + l_e_pi[1] + l_e_ep[0];
         } else {
             new_rho[0] = rl[b+0] + l_e_pi[0] + l_e_ep[1];
             new_rho[1] = rl[b+1] + l_e_pi[1] + l_e_ep[1];
             new_rho[2] = rl[b+2] + l_e_pi[0] + l_e_ep[0];
             new_rho[3] = rl[b+3] + l_e_pi[1] + l_e_ep[0];
         }
         return new_rho;
    }

    void update_hat(vector<double> nkt, vector<double>nkn, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, Hap3Param params) {
        if (ahat.size()==3) {
            ahat[0]=nkt[0]+params.mut_a0[0];
            ahat[1]=nkt[1]+nkt[3]+params.mut_a0[1];
            ahat[2]=nkt[2]+params.mut_a0[2];

            bhat[0]=nkt[3]+nkn[2]+nkn[3]+params.mut_b0[0];
            bhat[1]=nkt[1]+nkn[0]+nkn[1]+params.mut_b0[1];

            chat[0]=nkn[0]+nkn[2]+params.mut_c0[0];
            chat[1]=nkn[1]+nkn[3]+params.mut_c0[1];
        } else {
            ahat[0]=nkt[0]+nkt[2]+params.err_a0[0];
            ahat[1]=nkt[1]+nkt[3]+params.err_a0[1];

            bhat[0]=nkt[2]+nkt[3]+nkn[2]+nkn[3]+params.err_b0[0];
            bhat[1]=nkt[0]+nkt[1]+nkn[0]+nkn[1]+params.err_b0[1];

            chat[0]=nkn[0]+nkn[2]+params.err_c0[0];
            chat[1]=nkn[1]+nkn[3]+params.err_c0[1];
        }
    }

    void estimate(const vector<Haplotype> & haps, const vector<Read> & tumor_reads, const vector<Read> & normal_reads, const vector<double> & rlt, const vector<double> & rln, vector<double> & tumorHapFreqs, vector<double> & normalHapFreqs, vector <HapEstResult > & tumorPosteriors, vector <HapEstResult > & normalPosteriors, uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, Hap3Param params, string est_type, string log_prefix)
    {
        std::ofstream log;
#ifdef LOGDEBUG
        log.open((log_prefix+".log").c_str(), std::ios::out | std::ios::app);
        //est_type=mutation or non-mutation
        log << "MutationModel: " << est_type << endl;
        log << "mut a0:";print_vec(params.mut_a0, log);
        log << endl << "mut b0:";print_vec(params.mut_b0, log);
        log << endl << "mut c0:";print_vec(params.mut_c0, log);
        log << endl << "err a0:";print_vec(params.err_a0, log);
        log << endl << "err b0:";print_vec(params.err_b0, log);
        log << endl << "err c0:";print_vec(params.err_c0, log);
        log << endl;
#endif
        normalHapFreqs.clear();
        tumorHapFreqs.clear();
#ifndef MMTEST
        if(haps.size() != 4) { LOG(logDEBUG) << "error: only for 4 haps" << endl; return; }
#endif
        size_t nh=4;
        size_t nrt=rlt.size()/4;
        size_t nrn=rln.size()/4;

        vector<double> zt(nh*nrt,0.0); // expectations of read-haplotype indicator variables
        vector<double> zn(nh*nrn,0.0); // expectations of read-haplotype indicator variables
        vector<double> nkt(nh, 0.0), nkn(nh, 0.0);

        normalHapFreqs=nkn;tumorHapFreqs=nkt;

#ifndef MMTEST
        //filter reads
        std::set< PAV > allVariants;
        map<int, std::set<PAV> > allVariantsByPos;

        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    allVariants.insert(PAV(it->first,it->second));
                    allVariantsByPos[it->first].insert(PAV(it->first,it->second));
                }
            }
        }
#ifdef LOGDEBUG
        log << "#haplotype list" << endl;
        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            log << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    log << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            log << endl;
        }
#endif
        // create matrix of which variant is active in which set
        int idx=0;
        int nv=int(allVariants.size());
        vector<int> hapHasVar(nh*nv,0);

        BOOST_FOREACH(PAV pav, allVariants) {
            for (size_t h=0;h<nh;h++) {
                It it=haps[h].indels.find(pav.first);
                if (it!=haps[h].indels.end() && it->second.getString()==pav.second.getString()) hapHasVar[h*nv+idx]=1;
            }
            idx++;
        }
#endif
        // check haplotypes
        // run EM for this set of active variants
        bool converged=false;
        double tol=0.000001;double eOld=-HUGE_VAL, eNew;
        // initialize frequencies
        vector<double> ahat, bhat, chat;
        if (est_type=="mutation") {
            BOOST_FOREACH(double x, params.mut_a0) { ahat.push_back(x); }
            BOOST_FOREACH(double x, params.mut_b0) { bhat.push_back(x); }
            BOOST_FOREACH(double x, params.mut_c0) { chat.push_back(x); }
        } else {
            BOOST_FOREACH(double x, params.err_a0) { ahat.push_back(x); }
            BOOST_FOREACH(double x, params.err_b0) { bhat.push_back(x); }
            BOOST_FOREACH(double x, params.err_c0) { chat.push_back(x); }
        }

#ifdef LOGDEBUG
        log << "init:" << endl;
        log << "ahat:";print_vec(ahat, log);
        log << endl << "bhat:";print_vec(bhat, log);
        log << endl << "chat:";print_vec(chat, log);
#endif
        double llNew, llOld=-HUGE_VAL;
        int iter=0;
        lower_bound_t lb;
        while (!converged) {
            vector<double> l_e_pit = cal_dir_exp(ahat);
            vector<double> l_e_ep = cal_dir_exp(bhat);
            for(int i=0;i<nrt;i++) {
                vector<double> new_rho = make_newrho(rlt, l_e_pit, l_e_ep, i*4);
                vector<double> zt_i = norm_vec(exp_vec(new_rho));
                for(int k=0;k<4;k++) zt[i*4+k] = zt_i[k];
            }
            vector<double> l_e_pin = cal_dir_exp(chat);
            for(int i=0;i<nrn;i++) {
                vector<double> new_rho = make_newrho(rln, l_e_pin, l_e_ep, i*4);
                vector<double> zn_i = norm_vec(exp_vec(new_rho));
                for(int k=0;k<4;k++) zn[i*4+k] = zn_i[k];
            }
            nkt = sum_each_col(zt);nkn = sum_each_col(zn);
            update_hat(nkt, nkn, ahat, bhat, chat, params);
            compute_lower_bound(rlt, rln, ahat, bhat, chat, zt, zn, lb, params, log);
            llNew=lb.lower_bound;
            converged=(fabs(llOld-llNew))<tol || iter > 500;
#ifdef LOGDEBUG
            log << "### iter: " << iter << " llOld: " << llOld << " llNew: " << llNew << endl;

            print_posteriors(rlt, rln, ahat, bhat, chat, zt, zn, lb, params, log);
            print_params(tumor_reads, normal_reads, rlt, rln, ahat, bhat, chat, zt, zn, lb, params, log);
            log << endl;
#endif
            llOld=llNew;
            iter++;
        }
#ifdef LOGDEBUG
        log << "----------------finished[" << iter << "]---------------" << endl;
#endif
        print_posteriors(rlt, rln, ahat, bhat, chat, zt, zn, lb, params, log);
        print_params(tumor_reads, normal_reads, rlt, rln, ahat, bhat, chat, zt, zn, lb, params, log);

        vector<int> tmp_vec; //nonsense
        lb.compatible = tmp_vec;
        lb.ln_prior = 1.0;
        // check sum
        double zc=0.0;
        for (size_t h=0;h<nh;h++) {
            tumorHapFreqs[h]=nkt[h]/(double)nrt;
            normalHapFreqs[h]=nkn[h]/(double)nrn;
        }
        output_params(ahat, bhat, chat, zt, zn, lb, log_prefix);
#ifdef LOGDEBUG
        log << "==================results====================" << endl;
        log << "lb= " << lb.lower_bound << endl;
#endif
#ifndef MMTEST
#ifdef LOGDEBUG
        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            log << "hap[" << th << "] " << hap.seq << endl;
            log << tumorHapFreqs[th] << " " << normalHapFreqs[th];
            print_var_in_hap(hap,log);
        }

        log << endl;
#endif

        // compute marginal posteriors for the individual variants
        idx=-1;
        //log << "push posteriors" << endl;
        BOOST_FOREACH(PAV pav, allVariants) {
            idx++;
            double tumor_freq=0.0, normal_freq=0.0;
            for (size_t h=0;h<nh;h++) {
                if (hapHasVar[h*nv+idx]) {
                    tumor_freq+=tumorHapFreqs[h];
                    normal_freq+=normalHapFreqs[h];
                }
            }
            const AlignedVariant & avar = pav.second;
            const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
            int totnf=0, totnr=0;
            //LOG(logDEBUG) << pav.first << " " << pav.second.getString() << endl;
            if(av!=NULL) {
                tumorPosteriors.push_back(HapEstResult(pav.second, pav.first, -1.0, tumor_freq, totnf, totnr, av->info, nkt[0], nkt[1], nkt[2], nkt[3]));
                normalPosteriors.push_back(HapEstResult(pav.second, pav.first, -1.0, normal_freq, totnf, totnr, av->info, nkn[0], nkn[1], nkn[2], nkn[3]));
            }
        }
#endif
        best_lower_bound = lb;
    }
}
