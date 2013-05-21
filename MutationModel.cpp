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
#include "DInDel.hpp"
#include "Haplotype.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
#include "Utils.hpp"
#include "faidx.h"
#include "GetCandidates.hpp"
#include "ObservationModelSeqAn.hpp"
#include "VariantFile.hpp"
#include "Faster.hpp"
#include <ext/hash_map>
#include <exception>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/algorithm/string.hpp>
#include "MutationModel.hpp"

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
        double dig_sum = boost::math::digamma(sum);
        for(int i=0;i< ar.size();i++) {
            ans[i]=boost::math::digamma(ar[i]) - dig_sum;
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
    
    void print_var_in_hap(Haplotype hap) {
        for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
            if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                cout << "[" << it->second.getString() << " " << (it->first) << "]";
            }
        }
        cout << endl;
    }
    
    template <class X> void print_vec(vector<X> vec) {
        BOOST_FOREACH(X x, vec) { cout << x << " "; }
    }
    
    double cal_e_log_dir(vector<double> log_theta, vector<double> param) {
        //E[log p(theta|param)] : Dir
        //- sum(lgamma(param)) + lgamma(sum(param)) + sum((param-1) * log_theta)
        double sum_lgam_param=0.0;
        for(int i=0;i<param.size();i++) sum_lgam_param-=lgamma(param[i]);
        double lgam_sum_param=lgamma(sum_vec(param));
        double sum_par_lt=0.0;
        for(int i=0;i<param.size();i++) sum_par_lt+=(param[i]-1)*log_theta[i];
        if(log_theta.size() != param.size()) {
            cout << "--- cal_e_log_dir ---" << endl;
            print_vec(log_theta);cout << endl;
            print_vec(param);cout << endl;
            cout << sum_lgam_param << " " << lgam_sum_param << " " << sum_par_lt << endl;
            throw string("cal_e_log_dir");
        }    
        return sum_lgam_param + lgam_sum_param + sum_par_lt;
    }
    
    

    
    void compute_lower_bound(const vector<double> & rlt, const vector<double> & rln, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, Parameters params) {
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
        double e_log_p_pit = (ahat.size()==2) ? cal_e_log_dir(l_e_pit, params.c0) : cal_e_log_dir(l_e_pit, params.a0);
        double e_log_p_pin=cal_e_log_dir(l_e_pin, params.c0);
        double e_log_p_ep=cal_e_log_dir(l_e_ep, params.b0);
        
        double e_log_q_Zt=0.0, e_log_q_Zn=0.0;
        for(int i=0;i<zt.size();i++) e_log_q_Zt+=zt[i]*log(zt[i]);
        for(int i=0;i<zn.size();i++) e_log_q_Zn+=zn[i]*log(zn[i]);
        double e_log_q_pit=cal_e_log_dir(l_e_pit, ahat);
        double e_log_q_pin=cal_e_log_dir(l_e_pin, chat);
        double e_log_q_ep=cal_e_log_dir(l_e_ep, bhat);
        cout << "error";
        print_vec(l_e_pit);cout << endl;
        print_vec(ahat);cout << endl;
        cout << "lb: ";
        cout << e_log_p_Xt<< " " <<e_log_p_Xn<< " " <<e_log_p_Zt<< " " <<e_log_p_Zn<< " " <<e_log_p_pit<< " " <<e_log_p_pin << " " << e_log_p_ep << endl;
        cout  <<e_log_q_Zt<< " " <<e_log_q_Zn<< " " <<e_log_q_pit<< " " <<e_log_q_pin<< " " <<e_log_q_ep << endl;
        lb.lower_bound = e_log_p_Xt + e_log_p_Xn + e_log_p_Zt + e_log_p_Zn + e_log_p_pin +e_log_p_pit + e_log_p_ep - e_log_q_Zt - e_log_q_Zn - e_log_q_pit - e_log_q_pin - e_log_q_ep;
    }
   
    
    void print_params(const vector<Read> & tumor_reads, const vector<Read> & normal_reads, const vector<double> & rlt, const vector<double> & rln, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, Parameters params, bool print_z=false) {
        int nh = 4; //number of haps 
        int nrt = rlt.size() / nh; //number of reads of tumor
        int nrn = rln.size() / nh; //number of reads of normal
        cout << "##print params" << endl;
        if(print_z) {
            cout << "###zt###" << endl;
            for (int r=0;r<nrt;r++) {
#ifndef MMTEST
                cout << tumor_reads[r].seq_name;
#endif
                cout << " zt[" << r << "]:";
                for (int h=0;h<nh;h++) cout << " " << zt[r*nh+h];
                cout << endl;
            }
            cout << "###zn###" << endl;
            for (int r=0;r<nrn;r++) {
#ifndef MMTEST
                cout << normal_reads[r].seq_name;
#endif
                cout << " zn[" << r << "]:";
                for (int h=0;h<nh;h++) cout << " " << zn[r*nh+h];
                cout << endl;
            }
        }
        cout << "ahat: ";print_vec(ahat);
        cout << endl << "bhat: ";print_vec(bhat);
        cout << endl << "chat: ";print_vec(chat);
        cout << endl;
    }
    
    void print_posteriors(const vector<double> & rlt, const vector<double> & rln, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, vector<double> & zt, vector<double> & zn, lower_bound_t & lb, Parameters params) {
        cout << "##print posteriors" << endl;
        int nh = 4; //number of haps 
        int nrt = rlt.size() / nh; //number of reads of tumor
        int nrn = rln.size() / nh; //number of reads of normal
        vector<double> l_e_pit = cal_dir_exp(ahat);
        vector<double> l_e_ep = cal_dir_exp(bhat);
        vector<double> l_e_pin = cal_dir_exp(chat);
        cout << "pit: ";print_vec(exp_vec(l_e_pit));
        cout << endl << "ep: ";print_vec(exp_vec(l_e_ep));
        cout << endl << "pin: ";print_vec(exp_vec(l_e_pin));
        vector<double> nkt = sum_each_col(zt);
        vector<double> nkn = sum_each_col(zn);
        cout << endl << "nkt: ";print_vec(nkt);
        cout << endl << "nkn: ";print_vec(nkn);
        cout << endl;
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
    
    void update_hat(vector<double> nkt, vector<double>nkn, vector<double> &ahat, vector<double> &bhat, vector<double> &chat, Parameters params) {
        if (ahat.size()==3) {
            ahat[0]=nkt[0]+params.a0[0];
            ahat[1]=nkt[1]+nkt[3]+params.a0[1];
            ahat[2]=nkt[2]+params.a0[2];
            
            bhat[0]=nkt[3]+nkn[2]+nkn[3]+params.b0[0];
            bhat[1]=nkt[1]+nkn[0]+nkn[1]+params.b0[1];
            
            chat[0]=nkn[0]+nkn[2]+params.c0[0];
            chat[1]=nkn[1]+nkn[3]+params.c0[1];
        } else {
            ahat[0]=nkt[0]+nkt[2]+params.c0[0];
            ahat[1]=nkt[1]+nkt[3]+params.c0[1];
            
            bhat[0]=nkt[2]+nkt[3]+nkn[2]+nkn[3]+params.b0[0];
            bhat[1]=nkt[0]+nkt[1]+nkn[0]+nkn[1]+params.b0[1];
            
            chat[0]=nkn[0]+nkn[2]+params.c0[0];
            chat[1]=nkn[1]+nkn[3]+params.c0[1];
        }
    }
    
    void estimate(const vector<Haplotype> & haps, const vector<Read> & tumor_reads, const vector<Read> & normal_reads, const vector<double> & rlt, const vector<double> & rln, vector<double> & tumorHapFreqs, vector<double> & normalHapFreqs, vector <HapEstResult > & tumorPosteriors, vector <HapEstResult > & normalPosteriors, uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, int index, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & tumorVariantPosteriors, 
                  map<AlignedVariant, double> & normalVariantPosteriors, Parameters params, string est_type)
    {
        //est_type=mutation or non-mutation
        cout << "MutationModel: " << est_type << endl;
        cout << "a0:";print_vec(params.a0);
        cout << endl << "b0:";print_vec(params.b0);
        cout << endl << "c0:";print_vec(params.c0);
        cout << endl;
        normalHapFreqs.clear();
        tumorHapFreqs.clear();
#ifndef MMTEST
        if(haps.size() != 4) { cout << "error: only for 4 haps" << endl; return; }
#endif
        size_t nh=4;
        size_t nrt=rlt.size()/4;
        size_t nrn=rln.size()/4;
        
        vector<double> zt(nh*nrt,0.0); // expectations of read-haplotype indicator variables
        vector<double> zn(nh*nrn,0.0); // expectations of read-haplotype indicator variables
        vector<double> nkt(nh, 0.0), nkn(nh, 0.0);
        
        normalHapFreqs=nkn;tumorHapFreqs=nkt; 
       
#ifndef MMTEST
        // filter reads
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
        cout << "#haplotype list" << endl;
        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    cout << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            cout << endl;
        }

        cout << "#log haplotype and read likelihood" << endl;
        cout << "##tumor" << endl;
        for (size_t r=0;r<nrt;r++) {
            cout << tumor_reads[r].seq_name << " rl[" << r << "]:";
            for (size_t h=0;h<nh;h++) cout << " " << rlt[r*nh+h];
            cout << endl;
        }
        cout << "##normal" << endl;
        for (size_t r=0;r<nrn;r++) {
            cout << normal_reads[r].seq_name << " rl[" << r << "]:";
            for (size_t h=0;h<nh;h++) cout << " " << rln[r*nh+h];
            cout << endl;
        }

        // create matrix of which variant is active in which set
        int idx=0;
        int nv=int(allVariants.size());
        vector<int> hapHasVar(nh*nv,0);
        
        BOOST_FOREACH(PAV pav, allVariants) {
            cout << pav.first << " " << pav.second.getString() << " ";
            for (size_t h=0;h<nh;h++) {
                It it=haps[h].indels.find(pav.first);
                if (it!=haps[h].indels.end() && it->second.getString()==pav.second.getString()) hapHasVar[h*nv+idx]=1;
            }
            idx++;
        }
        cout << "allVariants: ";
        BOOST_FOREACH(PAV pav, allVariants) {
            cout << " [" << pav.first << " " << pav.second.getString() << "]";
        }
        cout << endl;
#endif
        // check haplotypes
        // run EM for this set of active variants
        bool converged=false;
        double tol=0.000001;double eOld=-HUGE_VAL, eNew;
        // initialize frequencies
        vector<double> ahat;
        if (est_type=="mutation") {
            BOOST_FOREACH(double x, params.a0) { ahat.push_back(x); }
        } else {
            BOOST_FOREACH(double x, params.c0) { ahat.push_back(x); }
        }
        vector<double> bhat(params.b0);
        vector<double> chat(params.c0);
        cout << "init:" << endl;
        cout << "ahat:";print_vec(ahat);
        cout << endl << "bhat:";print_vec(bhat);
        cout << endl << "chat:";print_vec(chat);
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
            compute_lower_bound(rlt, rln, ahat, bhat, chat, zt, zn, lb, params);
            llNew=lb.lower_bound;
            converged=(fabs(llOld-llNew))<tol || iter > 500;
            cout << "### iter: " << iter << " llOld: " << llOld << " llNew: " << llNew << endl;            

            print_posteriors(rlt, rln, ahat, bhat, chat, zt, zn, lb, params);
            print_params(tumor_reads, normal_reads, rlt, rln, ahat, bhat, chat, zt, zn, lb, params);
            cout << endl;
            llOld=llNew;
            iter++;
        }
        cout << "----------------finished[" << iter << "]---------------" << endl;
        print_posteriors(rlt, rln, ahat, bhat, chat, zt, zn, lb, params);
        print_params(tumor_reads, normal_reads, rlt, rln, ahat, bhat, chat, zt, zn, lb, params, true);
       
        vector<int> tmp_vec; //nonsense
        lb.compatible = tmp_vec;
        lb.ln_prior = 1.0;
        // check sum
        double zc=0.0;
        for (size_t h=0;h<nh;h++) {
            tumorHapFreqs[h]=nkt[h]/(double)nrt;
            normalHapFreqs[h]=nkn[h]/(double)nrn;
        }
        
        cout << "==================results====================" << endl;
        cout << "lb= " << lb.lower_bound << endl;
#ifndef MMTEST    
        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            cout << tumorHapFreqs[th] << " " << normalHapFreqs[th];
            print_var_in_hap(hap);
        }
        
        cout << endl;

        // compute marginal posteriors for the individual variants
        idx=-1;
        cout << "push posteriors" << endl;
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
            cout << pav.first << " " << pav.second.getString() << endl;
            if(av!=NULL) {
                tumorPosteriors.push_back(HapEstResult(pav.second, pav.first, -1.0, tumor_freq, totnf, totnr, av->info, nkt[0], nkt[1], nkt[2], nkt[3]));
                normalPosteriors.push_back(HapEstResult(pav.second, pav.first, -1.0, normal_freq, totnf, totnr, av->info, nkn[0], nkn[1], nkn[2], nkn[3]));
            }
        }
#endif
        best_lower_bound = lb;
    }
}