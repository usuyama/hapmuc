//
//  EMBasic.cpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/12/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//

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
#include "EMBasic.hpp"
#include "log.h"

namespace EMBasic {
    void computeLowerBound(const vector<double> & resp, const double a0, const vector<double> & ak, const vector<double> & ln_p_x_given_h, const vector<double> & lpi, const vector<int> & compatible, lower_bound_t & lb) {
        LOG(logDEBUG) << endl << "###computeLowerBound###" << endl;
        int nh = ak.size(); //number of haps
        int nr = resp.size() / nh; //number of reads
        LOG(logDEBUG) << "nh, nr = " << nh << ", " << nr << endl;
        LOG(logDEBUG) << "###compatible" << endl;
        for(int i=0;i < nh;i++) LOGP(logDEBUG) << " " << compatible[i];
        LOGP(logDEBUG) << endl;
        LOG(logDEBUG) << "lpi = ";
        for(int i=0;i < nh;i++) LOGP(logDEBUG) << " " << lpi[i];
        LOGP(logDEBUG) << endl << "a0=" << a0 << ", a = ";
        for(int i=0;i < nh;i++) LOGP(logDEBUG) << " " << ak[i];
        LOGP(logDEBUG) << endl;
        //ln_p_x_given_z
        double tmp = 0.0;
        for(int i=0;i < nr;i++) {
            for(int j=0;j < nh;j++) {
                if(!compatible[j]) continue;
                tmp += resp[i*nh+j] * ln_p_x_given_h[i*nh+j];
            }
        }
        lb.ln_p_x_given_z = tmp;
        //ln_p_z_given_pi
        tmp = 0.0;
        for(int i=0;i < nr;i++) {
            for(int j=0;j < nh;j++) {
                if(!compatible[j]) continue;
                tmp += resp[i*nh+j] * lpi[j];
            }
        }
        lb.ln_p_z_given_pi = tmp;
        //ln_p_pi
        tmp = 0.0;
        //lgamma(sum(a0)) - sum(lgamma(a0))
        double ahat = 0.0;
        for(size_t h=0;h<nh;h++) if (compatible[h]) {
            ahat += a0;
            tmp += lgamma(a0);
        }
        double ln_C = lgamma(ahat) - tmp;
        tmp = 0.0;
        // (a0 -1)Σln(pi)
        for(size_t h=0;h<nh;h++) if (compatible[h]) {
            tmp += (a0 - 1.0) * lpi[h];
        }
        LOG(logDEBUG) << "ln_C(a0) + tmp = " << ln_C << " + " << tmp << endl;
        lb.ln_p_pi = ln_C + tmp;
        //ln_q_z
        tmp = 0.0;
        for(int i=0;i < nr;i++) {
            for(int j=0;j < nh;j++) {
                if(!compatible[j]) continue;
                double t = resp[i*nh+j] * log(resp[i*nh+j]);
                if(!isnan(t)) tmp += t; //check NaN
            }
        }
        lb.ln_q_z = tmp;
        //ln_q_pi
        tmp = 0.0;
        //lgamma(sum(a)) - sum(lgamma(a))
        ahat = 0.0;
        for(size_t h=0;h<nh;h++) if (compatible[h]) {
            ahat += ak[h];
            tmp += lgamma(ak[h]);
        }
        ln_C = lgamma(ahat) - tmp;
        tmp = 0.0;
        // (a -1)Σln(pi)
        for(size_t h=0;h<nh;h++) if (compatible[h]) {
            tmp += (ak[h] - 1.0) * lpi[h];
        }
        LOG(logDEBUG) << "ln_C(a) + tmp = " << ln_C << " + " << tmp << endl;
        lb.ln_q_pi = ln_C + tmp;
        //
        lb.lower_bound = lb.ln_p_x_given_z + lb.ln_p_z_given_pi + lb.ln_p_pi - lb.ln_q_z - lb.ln_q_pi;
        printf("lower_bound=%lf, ln_p_x_given_z=%lf, ln_p_z_given_pi=%lf, ln_p_pi=%lf, ln_q_z=%lf, ln_q_pi=%lf\n", lb.lower_bound, lb.ln_p_x_given_z, lb.ln_p_z_given_pi, lb.ln_p_pi, lb.ln_q_z, lb.ln_q_pi);
        //
        LOG(logDEBUG) << "########" << endl << endl;
    }

    void estimate(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors, uint32_t candPos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & variantPosteriors, double a0, string program, Parameters params)
    {
        LOG(logDEBUG) << "EMBasic::estimate" << endl;
        LOG(logDEBUG).flush();
        // estimate haplotype frequencies using EM
        hapFreqs.clear();

        vector<lower_bound_t> lower_bounds;
        size_t nh=haps.size();
        size_t nr=reads.size();
        LOG(logDEBUG) << "nh, nr:" << nh << " " << nr << endl;
        if (nh == 0) {
            throw(string("num hap = 0 in EMBasic::estimate"));
        }
        vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods
        vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
        vector<double> pi(nh); // haplotype frequencies
        vector<double> lpi(nh); // expectation of log frequencies
        vector<double> nk(nh,0.0), ak(nh,0.0); // counts for each haplotype

        hapFreqs=nk;
        int numUnmappedRealigned=0;
        int idx=0;
        int numReadOffAllHaps=0;
        for (size_t r=0;r<nr;r++) {
            int offallhap=1;
            for (size_t h=0;h<nh;h++) {
                // initialize read-haplotype likelihoods
                rl[idx]=liks[h][r].ll;
                if (!liks[h][r].offHap) offallhap=0;
                idx++;
            }
            if (offallhap) {
                numReadOffAllHaps++;
            }else {
                if (reads[r].isUnmapped()) numUnmappedRealigned++;
            }
        }

        // filter reads
        vector<int> filtered(nh, 0);

        typedef pair<int, AlignedVariant> PAV;

        std::set< PAV > allVariants;
        map<int, std::set<PAV> > allVariantsByPos; //

        typedef map<int, AlignedVariant>::const_iterator It;
        typedef map<int, std::set<PAV> >::const_iterator PIt;

        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    allVariants.insert(PAV(it->first,it->second));
                    allVariantsByPos[it->first].insert(PAV(it->first,it->second));
                }
            }
        }

/*        LOG(logDEBUG) << "#read list" << endl;
        for (size_t r=0;r<nr;r++) {
            LOG(logDEBUG) << "r[" << r << "] " << reads[r].seq.seq << endl;
        }
        */

        LOG(logDEBUG) << "#haplotype list" << endl;
        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            LOG(logDEBUG) << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    LOGP(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            LOGP(logDEBUG) << endl;
        }

#ifdef LOGDEBUG
        LOG(logDEBUG) << "#log haplotype and read likelihood" << endl;
        for (size_t r=0;r<nr;r++) {
            LOGP(logDEBUG) << reads[r].seq_name << " rl[" << r << "]:";
            for (size_t h=0;h<nh;h++) {
                LOGP(logDEBUG) << " " << rl[r*nh+h];
            }
            LOGP(logDEBUG) << endl;
        }
#endif

        // set active variants, and divide into snps and indels
        vector< std::set< PAV > >  activeVariants, activeSNPs, activeIndels;

        int nav=0;
        int PRID=-1;
        if (program=="all") {
            std::set<PAV> snps, indels;
            BOOST_FOREACH(PAV pav, allVariants) {
                if (pav.second.isSNP()) {
                    snps.insert(pav);
                } else if (pav.second.isIndel()) {
                    indels.insert(pav);
                }
            }

            // both (double prior)
            activeVariants.push_back(allVariants);
            activeIndels.push_back(indels);
            activeSNPs.push_back(snps);
            nav++;
            PRID=1;
        } else if (program=="singlevariant") {
            std::set < std::set<PAV> > ssPAV;
            for (size_t h=0;h<haps.size();h++) if (filtered[h]==0) {
                const Haplotype & hap=haps[h];


                //LOG(logDEBUG) << "hap[" << h << "].seq: " << hap.seq << endl;

                //LOG(logDEBUG) << "vars:";
                std::set<PAV> act;
                for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                    if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                        act.insert(PAV(it->first, it->second));
                        //		LOG(logDEBUG) << "[" << it->first << "," << it->second.getString() << "]";
                    }
                }
                //LOG(logDEBUG) << endl;
                ssPAV.insert(act);
            }
            //make fullset
            std::set<PAV> act;
            BOOST_FOREACH(PAV pav, allVariants) {
                act.insert(pav);
            }
            ssPAV.insert(act);

            nav=0;
            BOOST_FOREACH(std::set<PAV> s, ssPAV) {
                std::set<PAV> snps, indels;
                BOOST_FOREACH(PAV pav, s) {
                    if (pav.second.isSNP()) {
                        snps.insert(pav);
                    } else if (pav.second.isIndel()) {
                        indels.insert(pav);
                    }
                }

                // both (double prior)
                activeVariants.push_back(s);
                activeIndels.push_back(indels);
                activeSNPs.push_back(snps);
                nav++;
            }

            PRID=2;
        } else if (program == "priorpersite") {
            nav = 0;

            // add reference haplotype
            activeVariants.push_back(std::set<PAV>());
            activeIndels.push_back(std::set<PAV>());
            activeSNPs.push_back(std::set<PAV>());



            for (map<int, std::set<PAV> >::const_iterator site_it = allVariantsByPos.begin();site_it!=allVariantsByPos.end();site_it++) {
                std::set<PAV> snps, indels;
                BOOST_FOREACH(PAV pav, site_it->second) {
                    if (pav.second.isSNP()) snps.insert(pav);
                    else if (pav.second.isIndel()) indels.insert(pav);
                }

                int maxStateSnp = (snps.size()==0)?1:2;
                int maxStateIndel = (indels.size()==0)?1:2;

                int prevNumActive = activeVariants.size();
                int sSnp = 1, sIndel = 1;
                //for (int sSnp = 0;sSnp<maxStateSnp;sSnp++) {
                //	for (int sIndel = 0; sIndel < maxStateIndel; sIndel++) {
                //
                //		if (sSnp == 1 || sIndel == 1) {
                // extend previous activeVariants

                for (int pna = 0;pna<prevNumActive;pna++) {
                    std::set<PAV> av = activeVariants[pna];
                    std::set<PAV> aIndels = activeIndels[pna];
                    std::set<PAV> aSNPs = activeSNPs[pna];
                    if (sSnp == 1) {
                        av.insert(snps.begin(),snps.end());
                        aSNPs.insert(snps.begin(),snps.end());
                    }
                    if (sIndel == 1) {
                        av.insert(indels.begin(),indels.end());
                        aIndels.insert(indels.begin(),indels.end());
                    }

                    activeVariants.push_back(av);
                    activeIndels.push_back(aIndels);
                    activeSNPs.push_back(aSNPs);
                }

                //		}
                //		}
                //	}

            }
            nav = activeVariants.size();

            PRID = 3;

        } else {
            cerr << "Unknown EM option" << endl;
            exit(1);
        }

        vector<int> compatible(nh,0);
        vector<double> logliks(nav,0.0);
        vector<double> logpriors(nav, 0.0);
        vector<double> post(nav,0.0);
        vector<double> freqs(nav*nh,0.0);

        // create matrix of which variant is active in which set
        idx=0;
        int nv=int(allVariants.size());
        vector<int> active(nav*nv,0), hapHasVar(nh*nv,0);

        //LOG(logDEBUG) << "active: " << endl;
        BOOST_FOREACH(PAV pav, allVariants) {
            //  LOG(logDEBUG) << pav.first << " " << pav.second.getString() << " ";
            for (int a=0;a<nav;a++) {
                if (activeVariants[a].find(pav)!=activeVariants[a].end()) active[a*nv+idx]=1;
                //    LOG(logDEBUG) << " " << active[a*nv+idx];
            }
           // LOG(logDEBUG) << endl;
            for (size_t h=0;h<nh;h++) {
                It it=haps[h].indels.find(pav.first);
                if (it!=haps[h].indels.end() && it->second.getString()==pav.second.getString()) hapHasVar[h*nv+idx]=1;
            }
            for(int i=0;i<nv;i++) {
                for(int j=0;j<nh;j++) {
                    //      LOG(logDEBUG) << "[" << active[i*nv+idx] << " " << hapHasVar[j*nv+idx] << "]";
                }
                //LOG(logDEBUG) << endl;
            }
            idx++;
        }


        LOG(logDEBUG) << "allVariants: ";
        BOOST_FOREACH(PAV pav, allVariants) {
            LOGP(logDEBUG) << " [" << pav.first << " " << pav.second.getString() << "]";
        }
        LOGP(logDEBUG) << endl;

        vector<vector<int> > compatibles;
        vector<vector<double> > aks;
        LOG(logDEBUG) << "#subsets" << endl;
        for (int th=0;th<nav;th++) {
            LOGP(logDEBUG) << "[" << th << "] ";
            for (size_t h=0;h<nh;h++) {
                compatible[h]=1;
                if (filtered[h]!=0) {
                    compatible[h]=0;
                } else {
                    for (It it=haps[h].indels.begin();it!=haps[h].indels.end();it++) {
                        if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D') && activeVariants[th].find(PAV(it->first,it->second))==activeVariants[th].end()) {
                            // haplotype h has a non-ref variant that is not one of the active variants
                            compatible[h]=0;
                            break;
                        }
                    }
                }
                LOGP(logDEBUG) << " " << compatible[h];
            }
            compatibles.push_back(compatible);
            LOGP(logDEBUG) << endl;
        }



        double logz=-HUGE_VAL;

        for (int th=0;th<nav;th++) {

            // set active variants

            double logprior=0.0;

            map <int, int> sites;
            BOOST_FOREACH(PAV pav, activeSNPs[th]) {
                sites[pav.first]=1;

                const AlignedVariant & avar = pav.second;
                const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
                int type = pav.second.getType();

                double lnf = 0.0;

                if (type==Variant::SNP) lnf = log(params.priorSNP); else if (type==Variant::DEL || type==Variant::INS) lnf = log(params.priorIndel);
                if (av==NULL) {
                    logprior += lnf;
                } else {
                    double prior = av->getFreq();
                    if (prior<0.0) logprior += lnf; else logprior+=log(prior);
                }

            }
            BOOST_FOREACH(PAV pav, activeIndels[th]) {

                const AlignedVariant & avar = pav.second;
                const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
                int type = pav.second.getType();

                double lnf = 0.0;

                if (type==Variant::SNP) lnf = log(params.priorSNP); else if (type==Variant::DEL || type==Variant::INS) lnf = log(params.priorIndel);
                if (av==NULL) {
                    logprior += lnf;
                } else {
                    double prior = av->getFreq();
                    if (prior<0.0) logprior += lnf; else logprior+=log(prior);
                }

                sites[pav.first]=2;
            }


            /*
             for (map<int,int>::const_iterator it=sites.begin();it!=sites.end();it++) {
             if (it->second==2) logprior+=log(params.priorIndel); else if (it->second==1) logprior+=log(params.priorSNP);
             }
             */

            logpriors[th]=logprior;

            //		LOG(logDEBUG) << "Number of indels: " << ni << " number of SNPs: " << ns << endl;

            // check haplotypes

            int numah=0; // number of haplotypes for which frequencies will be estimated
            for (size_t h=0;h<nh;h++) {
                compatible[h]=1;
                if (filtered[h]!=0) {
                    compatible[h]=0;
                } else {
                    for (It it=haps[h].indels.begin();it!=haps[h].indels.end();it++) {
                        if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D') && activeVariants[th].find(PAV(it->first,it->second))==activeVariants[th].end()) {
                            // haplotype h has a non-ref variant that is not one of the active variants
                            compatible[h]=0;
                            break;
                        }
                    }
                }
                if (compatible[h]) numah++;
            }


            // run EM for this set of active variants
            bool converged=false;
            double tol=params.EMtol;

            double eOld=-HUGE_VAL, eNew;

            // initialize frequencies
            for (size_t h=0;h<nh;h++) if (compatible[h]) lpi[h]=log(1.0/double(numah)); else lpi[h]=-100;


            double loglik, llNew, llOld=-HUGE_VAL;
            int iter=0;
            while (!converged) {
                //LOG(logDEBUG) << endl << "EM[" << iter << "]:" << endl;
                // compute expectation of indicator variables
                for (size_t h=0;h<nh;h++) nk[h]=0.0;

                loglik=0.0;

                int idx=0;
                for (size_t r=0;r<nr;r++) {
                    double lognorm=-HUGE_VAL;
                    // compute responsibilities
                    for (size_t h=0;h<nh;h++) {
                        z[idx]=lpi[h]+(rl[idx]);
                        lognorm=addLogs(lognorm, z[idx]);
                        idx++;
                    }
                    idx-=nh;
                    // normalize and compute counts
                    for (size_t h=0;h<nh;h++) {
                        z[idx]-=lognorm;
                        z[idx]=exp(z[idx]);
                        nk[h]+=z[idx];
                        idx++;
                    }
                    loglik+=lognorm;

                }
                // compute frequencies
                //LOG(logDEBUG) << "pi: ";

                double ahat=0.0;
                for (size_t h=0;h<nh;h++) if (compatible[h]) {
                    ak[h]=nk[h]+a0; // a0 is Dirichlet prior parameter
                    ahat+=ak[h];
                }
                try {
                double dahat=boost::math::digamma(ahat);

                for (size_t h=0;h<nh;h++) {

                    // do variational bayes
                    if (compatible[h]) {
                        lpi[h]=boost::math::digamma(ak[h])-dahat;
                        pi[h]=log((a0+nk[h])/(double(numah)*a0+double(nr)));
                    } else {
                        lpi[h]=-100;
                        pi[h]=-100;
                    }
                    //	LOG(logDEBUG) << " " << pi[h];
                    //	zp+=exp(pi[h]);
                }

                } catch (std::exception& e) {
                    string message = string("error_exception_").append(e.what());
                    LOG(logDEBUG) << message << endl;
                    LOGP(logDEBUG) << "EMBasic digamma" << endl << "ak: ";
                    BOOST_FOREACH(double x, ak) { LOGP(logDEBUG) << x << " "; }
                    LOGP(logDEBUG) << endl;
                    LOGP(logDEBUG) << "ahat: " << ahat << endl;
                    LOGP(logDEBUG).flush();
                    throw;
                }

                //LOG(logDEBUG) << " zp: " << zp << endl;


                idx=0;
                eNew=0.0;
                for (size_t r=0;r<nr;r++) {
                    for (size_t h=0;h<nh;h++) {
                        // compute responsibilities
                        eNew+=z[idx]*( pi[h]+rl[idx]);
                        idx++;
                    }
                }
                //LOG(logDEBUG) << "eOld: " << eOld << " " << " eNew: " << eNew; for (size_t h=0;h<nh;h++) { LOG(logDEBUG) << " " << pi[h]; } LOG(logDEBUG) << endl;
                /*
                 for (size_t r=0;r<nr;r++) {
                 LOG(logDEBUG) << "z[" << r << "]:";
                 for (size_t h=0;h<nh;h++) {
                 LOG(logDEBUG) << " " << z[r*nh+h];
                 }
                 LOG(logDEBUG) << endl;
                 }
                 */
                llNew=loglik;
                //LOG(logDEBUG) << "loglik: " << loglik << endl;
                if (0 && llOld>llNew+1e-5)  {
                    cerr << "ERROR: nr: " << nr << " eOld: " << eOld << " eNew: " << eNew << " diff: " << eOld-eNew << endl;
                    LOG(logDEBUG) << "ERROR: nr: " << nr << " eOld: " << eOld << " eNew: " << eNew << " diff: " << eOld-eNew << endl;
                    cerr << "ERROR: nr: " << nr << " llOld: " << llOld << " eNew: " << llNew << " diff: " << llOld-llNew << endl;
                    LOG(logDEBUG) << "ERROR: nr: " << nr << " llOld: " << llOld << " eNew: " << llNew << " diff: " << llOld-llNew << endl;

                    //throw string("EM Error in estimateHapFreqs");
                    //iter=100;

                }
                converged=(fabs(eOld-eNew))<tol || iter > 500;
                //LOG(logDEBUG) << "iter: " << iter << " eOld: " << eOld << " eNew: " << eNew << endl;

                eOld=eNew;
                llOld=llNew;
                iter++;


            }
            LOG(logDEBUG) << "----------------finished[" << iter << "]---------------" << endl;
            LOG(logDEBUG) << "###compatible" << endl;
            for(int i=0;i < nh;i++) LOGP(logDEBUG) << " " << compatible[i];
            LOGP(logDEBUG) << endl;
            LOG(logDEBUG) << "lpi = ";
            for(int i=0;i < nh;i++) LOGP(logDEBUG) << " " << lpi[i];
            LOGP(logDEBUG) << endl << "a0=" << a0 << ", a = ";
            for(int i=0;i < nh;i++) LOGP(logDEBUG) << " " << ak[i];
            aks.push_back(ak);

            LOG(logDEBUG) << endl;
            lower_bound_t lb;
            computeLowerBound(z, a0, ak, rl, lpi, compatible, lb);
            //  lower_bound_t another_lb;
            //  computeAnotherLowerBound(z, a0, ak, rl, lpi, compatible, another_lb);
            /*  if(abs(lb.lower_bound - another_lb.lower_bound) > 0.01) {
             LOG(logDEBUG) << "error: " << "another lower bound"  << another_lb.lower_bound << endl;
             LOG(logDEBUG) << "normal lower bound" << lb.lower_bound << endl;
             throw 20;
             }*/
            lb.compatible = compatible;
            lb.ln_prior = logprior;
            lower_bounds.push_back(lb);

#ifdef LOGDEBUG
            for (int r=0;r<nr;r++) {
                LOG(logDEBUG) << reads[r].seq_name << " z[" << r << "]:";
                for (int h=0;h<nh;h++) {
                    LOGP(logDEBUG) << " " << z[r*nh+h];
                }
                LOGP(logDEBUG) << endl;
            }
#endif

            // check sum

            double zc=0.0;
            for (size_t y=0;y<nh;y++) {
                zc+=exp(pi[y]);
            }

            if (0) {
                LOG(logDEBUG) << "th: " << th << endl;
                for (size_t y=0;y<nh;y++) {
                    LOGP(logDEBUG) << "[" << y << "," << exp(pi[y]) << "]";
                }
                LOGP(logDEBUG) << endl << endl;
            }

            logliks[th]=loglik;
            logz=addLogs(logz, lower_bounds[th].lower_bound+logprior);
            for (size_t h=0;h<nh;h++) { freqs[th*nh+h]=exp(pi[h])/zc; }
            for (size_t h=0;h<nh;h++) { LOGP(logDEBUG) << " " << freqs[th*nh+h]; } LOGP(logDEBUG) << endl;

            LOG(logDEBUG) << "lower_bound: " << lower_bounds[th].lower_bound << " loglik: " << loglik << " " << logliks[th] << " logprior: " << logprior << endl << endl;

        }


        for (int a=0;a<nav;a++) {
            post[a]=exp(lower_bounds[a].lower_bound+logpriors[a]-logz);
            LOG(logDEBUG) << "post[" << a << "]: " << post[a] << endl;
        }

        for (size_t h=0;h<nh;h++) {
            hapFreqs[h]=0.0;
        }

        for (int th=0;th<nav;th++) for (size_t h=0;h<nh;h++) {
            hapFreqs[h]+=exp(lower_bounds[th].lower_bound+logpriors[th]-logz)*freqs[th*nh+h];
        }

        LOG(logDEBUG) << "==================results====================" << endl;

        for (int a=0;a<nav;a++) {
            LOG(logDEBUG) << "[" << a << "] " << lower_bounds[a].lower_bound << " " << logpriors[a] << ";";
            for(int j=0;j<nh;j++) {
                LOGP(logDEBUG) << " " << compatibles[a][j];
            }
            LOGP(logDEBUG) << "; (freqs)";
            for(int j=0;j<nh;j++) {
                LOGP(logDEBUG) << " " << freqs[a*nh+j];
            }
            LOGP(logDEBUG) << endl;
        }

        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            LOG(logDEBUG) << "hap[" << th << "] " << hap.seq << endl;
            LOGP(logDEBUG) << hapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    LOGP(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            LOGP(logDEBUG) << endl;
        }

        LOG(logDEBUG) << endl;

        // compute marginal posteriors for the individual variants
        vector< std::set<int> > readidx(2); //TODO
        for (int r=0;r<int(nr);r++) readidx[reads[r].poolID].insert(r);

        vector<double> prior(nh*nh,0.0);

        idx=-1;
        BOOST_FOREACH(PAV pav, allVariants) {
            idx++;

            double logp=-HUGE_VAL;
            double freq=0.0;
            for (int th=0;th<nav;th++) {
                if (active[th*nv+idx]) {
                    logp=addLogs(logp, lower_bounds[th].lower_bound+logpriors[th]);
                }
            }

            for (size_t h=0;h<nh;h++) {
                if (hapHasVar[h*nv+idx]) {
                    freq+=hapFreqs[h];
                }
            }

            logp-=logz;
            bool doGLF=false; //(candPos==leftPos+pav.first)?true:false;
            const AlignedVariant & avar = pav.second;
            const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
            int totnf=0, totnr=0;
            if(av!=NULL) {
                if(nh == 1) {
                    posteriors.push_back(HapEstResult(pav.second, pav.first,exp(logp),freq, totnf, totnr, av->info, nk[0], 0.0, 0.0, 0.0));
                } else if(nh == 2) {
                    posteriors.push_back(HapEstResult(pav.second, pav.first,exp(logp),freq, totnf, totnr, av->info, nk[0], nk[1], 0.0, 0.0));
                } else {
                    posteriors.push_back(HapEstResult(pav.second, pav.first,exp(logp),freq, totnf, totnr, av->info, nk[0], nk[1], nk[2], 0.0));
                }
            }
        }

        //select best model
        double best = -HUGE_VAL;
        lower_bound_t lb;
        BOOST_FOREACH(lower_bound_t x, lower_bounds) {
            if(best < x.lower_bound) {
                best = x.lower_bound;
                lb = x;
            }
        }
        best_lower_bound = lb;


        stringstream os; os << params.fileName << "." << params.tid << "." << candPos;
        /* debug HMM */
        string oprefix = os.str();

        // string fname = oprefix;
        // fname.append(".hapvars");
        // ofstream of(fname.c_str());
        /* if (!of.is_open()) {
         throw string("Cannot open file ").append(fname).append(" for writing .hapvars file");
         }
         // output all likelihoods

         for (size_t r=0;r<nr;r++) {
         of << r << " " << bam1_qname(reads[r].getBam()) << " " << log(1.0-reads[r].mapQual) << " " << reads[r].poolID;
         for (size_t h=0;h<nh;h++) {
         of << " " << liks[h][r].ll;
         }
         for (size_t h=0;h<nh;h++) {
         of << " " << liks[h][r].offHap;
         }
         of << endl;
         }
         of.close();
         */
        // output alignments
        /*
         fname = oprefix;
         fname.append(".alignments");
         LOG(logDEBUG) << "fname: " << fname << endl;*/
        // of.open(fname.c_str());
        /*if (!of.is_open()) {
         throw string("Cannot open file ").append(fname).append(" for writing .liks file");
         }
         */

        /*
         for (size_t r=0;r<nr;r++) {
         LOG(logDEBUG) << "###" << endl;
         LOG(logDEBUG) << "read: " << bam1_qname(reads[r].getBam()) << " mpos: " << reads[r].matePos << endl;
         LOG(logDEBUG) << "isUnmapped: " << reads[r].isUnmapped() << endl;
         // compute maximum alignment likelihood
         double lq = 0.0;
         for (size_t b=0;b<reads[r].size();b++) lq += log(reads[r].qual[b]);

         LOG(logDEBUG) << "Max alignment loglik: " << lq << endl;

         double maxll = -HUGE_VAL;
         std::set <int> mlhaps;
         for (int h=nh-1;h>=0;h--) if (liks[h][r]>maxll) { maxll = liks[h][r]; }
         for (int h=nh-1;h>=0;h--) mlhaps.insert(h); //if (fabs(liks[h][r]-maxll)<0.01) mlhaps.insert(h);
         BOOST_FOREACH(int hidx, mlhaps) {
         LOG(logDEBUG) << "r: " << r << " hidx: " << hidx << " maxll:" << maxll << endl;
         ObservationModelFBMaxErr obs(haps[hidx], reads[r], leftPos, params.obsParams);
         LOG(logDEBUG) << string(50,' ') << haps[hidx].seq << endl;
         obs.printAlignment(50);
         }
         }
         */

        //of.close();
    }

}
