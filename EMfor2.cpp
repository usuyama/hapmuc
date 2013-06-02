//
//  EMfor2.cpp
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
#include "EMfor2.hpp"

namespace EMfor2 {
    void computeLowerBound(const vector<double> & resp, const double a0, const double b0, const double ahat, const double bhat, const vector<double> & ln_p_x_given_h, const vector<double> & lpi, lower_bound_t & lb) {
        cout << endl << "###computeAnotherLowerBound###" << endl;
        int nh = 2; //number of haps 
        int nr = resp.size() / nh; //number of reads
        cout << "lpi = "; 
        for(int i=0;i < nh;i++) cout << " " << lpi[i];
        cout << endl << "a0=" << a0 << ", b0 = " << b0 << ", ahat = " << ahat << ", bhat = " << bhat << endl;
        cout << endl;
        //ln_p_x_given_z
        double tmp = 0.0;
        for(int i=0;i < nr;i++) {
            for(int j=0;j < nh;j++) {
                tmp += resp[i*nh+j] * ln_p_x_given_h[i*nh+j];
            }
        }
        lb.ln_p_z_given_pi = tmp;
        //entropy
        double entropy = 0.0;
        for(int i=0;i < nr;i++) {
            for(int j=0;j < nh;j++) {
                double t = resp[i*nh+j] * log(resp[i*nh+j]); 
                if(!isnan(t)) entropy -= t; //check NaN
            }
        }
        //ln_c
        double ln_c_Bhat = lbeta(ahat, bhat);
        double ln_c_B0 = lbeta(a0, b0);
        lb.lower_bound = lb.ln_p_z_given_pi + entropy + ln_c_Bhat - ln_c_B0;
        cout << "another lb=" << lb.lower_bound << ", lb.ln_p_z_given_pi=" << lb.ln_p_z_given_pi << ", entropy=" << entropy << ", ln_c_Bhat(" << ln_c_Bhat << ") - ln_c_B0(" << ln_c_B0 << ")=" << ln_c_Bhat - ln_c_B0 << endl;
    }
    

    void estimate_basic(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, lower_bound_t & best_lower_bound, map<AlignedVariant, double> & variantPosteriors, Parameters params, double a0, double b0)
    {
        cout << "EMfor2" << endl;
        hapFreqs.clear();
        if(haps.size() != 2) { cout <<"error: only for 2 haps" << endl; return; }
        vector<lower_bound_t> lower_bounds;
        size_t nh=2;
        size_t nr=reads.size();
        
        vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods
        vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
        vector<double> pi(nh); // haplotype frequencies
        vector<double> lpi(nh); // expectation of log frequencies
        vector<double> nk(nh,0.0); // counts for each haplotype
        double ahat, bhat;
        
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
        /*
        cout << "#read list" << endl;
        for (size_t r=0;r<nr;r++) {
            cout << "r[" << r << "] " << reads[r].seq.seq << endl;
        }*/
        
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
        for (size_t r=0;r<nr;r++) {
            cout << reads[r].seq_name << " rl[" << r << "]:";
            for (size_t h=0;h<nh;h++) {
                cout << " " << rl[r*nh+h];
            }
            cout << endl;
        }
        
        
        
        double logprior = 1.0;
        double post;
        vector<double> freqs(nh,0.0);
        
        // create matrix of which variant is active in which set
        idx=0;
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
        double logz=-HUGE_VAL;
        //		cout << "Number of indels: " << ni << " number of SNPs: " << ns << endl;
        
        // check haplotypes
        // run EM for this set of active variants
        bool converged=false;
        double tol=params.EMtol;
        
        double eOld=-HUGE_VAL, eNew;
        
        // initialize frequencies
        double ave_ep = a0 / (a0 + b0);
        lpi[1] = log(ave_ep);
        lpi[0] = log(1-ave_ep);
        
        
        double loglik, llNew, llOld=-HUGE_VAL;
        int iter=0;
        while (!converged) {
            cout << endl << "EM[" << iter << "]:" << endl;
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
            cout << "pi: ";
            ahat = nk[1] + a0;
            bhat = nk[0] + b0;
            ave_ep = ahat / (ahat + bhat);
            pi[1] = ave_ep;
            pi[0] = 1 - ave_ep;
            for(int l=0;l < nh;l++) { lpi[l] = log(pi[l]); }
            idx=0;
            eNew=0.0;
            for (size_t r=0;r<nr;r++) {
                for (size_t h=0;h<nh;h++) {
                    eNew+=z[idx]*(lpi[h]+rl[idx]);
                    idx++;
                }
            }
            cout << "eOld: " << eOld << " " << " eNew: " << eNew; for (size_t h=0;h<nh;h++) { cout << " " << pi[h]; } cout << endl;
            cout << "a0:" << a0 << " b0:" << b0 << " ahat:" << ahat << " bhat:" << bhat << endl;
            /*
             for (size_t r=0;r<nr;r++) {
             cout << "z[" << r << "]:";
             for (size_t h=0;h<nh;h++) {
             cout << " " << z[r*nh+h];
             }
             cout << endl;
             }
            */
            llNew=loglik;
            cout << "loglik: " << loglik << endl;
            if (0 && llOld>llNew+1e-10)  {
                cerr << "ERROR: nr: " << nr << " eOld: " << eOld << " eNew: " << eNew << " diff: " << eOld-eNew << endl;
                cout << "ERROR: nr: " << nr << " eOld: " << eOld << " eNew: " << eNew << " diff: " << eOld-eNew << endl;
                cerr << "ERROR: nr: " << nr << " llOld: " << llOld << " eNew: " << llNew << " diff: " << llOld-llNew << endl;
                cout << "ERROR: nr: " << nr << " llOld: " << llOld << " eNew: " << llNew << " diff: " << llOld-llNew << endl;
                
                //throw string("EM Error in estimateHapFreqs");
                //iter=100;
                
            }
            converged=(fabs(eOld-eNew))<tol || iter > 500;
            cout << "iter: " << iter << " eOld: " << eOld << " eNew: " << eNew << endl;
            
            eOld=eNew;
            llOld=llNew;
            iter++;
            
            
        }
        cout << "----------------finished[" << iter << "]---------------" << endl;
        cout << "lpi = "; 
        for(int i=0;i < nh;i++) cout << " " << lpi[i];
        
        cout << endl;
        lower_bound_t lb;
        computeLowerBound(z, a0, b0, ahat, bhat, rl, lpi, lb);
        //  lower_bound_t another_lb;
        //  computeAnotherLowerBound(z, a0, ak, rl, lpi, compatible, another_lb);
        /*  if(abs(lb.lower_bound - another_lb.lower_bound) > 0.01) {
         cout << "error: " << "another lower bound"  << another_lb.lower_bound << endl;
         cout << "normal lower bound" << lb.lower_bound << endl;
         throw 20;
         }*/
        vector<int> tmp_vec; //nonsense
        lb.compatible = tmp_vec;
        lb.ln_prior = logprior;
        
        for (int r=0;r<nr;r++) {
            cout << reads[r].seq_name << " z[" << r << "]:";
            for (int h=0;h<nh;h++) {
                cout << " " << z[r*nh+h];
            }
            cout << endl;
        }
        
        // check sum
        
        double zc=0.0;
        for (size_t h=0;h<nh;h++) { freqs[h]=pi[h]; }
        for (size_t h=0;h<nh;h++) { cout << " " << freqs[h]; } cout << endl;
        
        cout << "lower_bound: " << lb.lower_bound << " loglik: " << loglik << " " << loglik << " logprior: " << logprior << endl << endl;
        post=exp(lb.lower_bound+logprior-logz);
        cout << "post["  << "]: " << post << endl;
        
        for (size_t h=0;h<nh;h++) {
            hapFreqs[h]+=freqs[h];
        }
        
        cout << "==================results====================" << endl;
        
        cout << "[" << "] " << lb.lower_bound << " " << logprior << ";";
        
        cout << "; (freqs)";
        for(int j=0;j<nh;j++) {
            cout << " " << freqs[j];
        }
        cout << endl;
        
        
        for (size_t th=0;th<nh;th++) {
            const Haplotype & hap=haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            cout << hapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    cout << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            cout << endl;
        }
        
        cout << endl;
        
        // compute marginal posteriors for the individual variants
        vector< std::set<int> > readidx(2); //TODO
        for (int r=0;r<int(nr);r++) readidx[reads[r].poolID].insert(r);
        idx=-1;
        cout << "push posteriors" << endl;
        BOOST_FOREACH(PAV pav, allVariants) {
            idx++;
            double freq=0.0;
            for (size_t h=0;h<nh;h++) {
                cout << hapHasVar[h*nv+idx] << " " << freq << endl;
                if (hapHasVar[h*nv+idx]) {
                    freq+=hapFreqs[h];
                }
            }
            bool doGLF=false; //(candPos==leftPos+pav.first)?true:false;
            const AlignedVariant & avar = pav.second;
            const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
            int totnf=0, totnr=0;
            cout << pav.first << " " << pav.second.getString() << endl;
            if(av!=NULL) {
                posteriors.push_back(HapEstResult(pav.second, pav.first, -1.0, freq, totnf, totnr, av->info, nk[0], nk[1], 0.0, 0.0));
            }
        }

        //select best model
        best_lower_bound = lb;
        
    }
    
}