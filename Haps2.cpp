//
//  Created by 直人 臼山 on 10/10/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//

#ifndef HAPS2_CPP_
#define HAPS2_CPP_

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
#include "foreach.hpp"
#include "bam.h"
#include "DInDel.hpp"
#include "Haplotype.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
#include "Utils.hpp"
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
#include "MutationCall.hpp"
#include "Haps.hpp"
#include "Haps2.hpp"
#include "EMBasic.hpp"

using namespace seqan;
namespace Haps2 {
    struct HapFreqPair{
        int index;
        double freq;
    };
    
    class GreaterHapFreq {
    public:
        bool operator()(const HapFreqPair& riLeft, const HapFreqPair& riRight) const {
            return riLeft.freq > riRight.freq;        }
    };
    
    void filter_reads(vector<Read> & reads, vector<vector<MLAlignment> > & liks) {
        cout << "##filter reads" << endl;
        int nh = liks.size();
        int i = 0;
        while(i < reads.size()) {
            double max, min;
            max = liks[0][i].ll;
            min = max;
            for(int k = 1;k < nh;k++) {
                double l = liks[k][i].ll;
                if(l > max) max = l;
                if(l < min) min = l;
            }
            cout << i << ": max, min " << max << ", " << min << endl;
            if(max - min < 0.01) {
                // we couldnt distinguish sources of reads
                reads.erase(reads.begin()+i);
                for(int k = 0;k < nh;k++)
                    liks[k].erase(liks[k].begin()+i);
            } else {
                i++;
            }
        }
    }
    
    string get_hap(string refseq, int left, int right, vector<AlignedVariant> variants) {
        vector<int> ref_pos(right-left+1);
        for(int i = 0; i < ref_pos.size()-1; i++) ref_pos[i] = left + i;
        string seq = refseq;
        BOOST_FOREACH(AlignedVariant var, variants) {
            vector<int>::iterator it = find(ref_pos.begin(), ref_pos.end(), var.getStartHap());
            if(it!=ref_pos.end()) {
                int idx = distance(ref_pos.begin(), it);
                bool changed = false;
                if (var.getType()==Variant::DEL) {
                    seq.erase(idx, var.size());
                    ref_pos.erase(ref_pos.begin()+idx, ref_pos.begin()+idx+var.size());
                    changed=true;
                } else if (var.getType()==Variant::INS) {
                    seq.insert(idx, var.getSeq());
                    ref_pos.insert(ref_pos.begin()+idx, (size_t) var.size(), -1);
                    changed=true;
                } else if (var.getType()==Variant::SNP) {
                    char nuc = var.getSeq()[3];
                    if (seq[idx]!=nuc) {
                        seq[idx]=nuc;
                        changed=true;
                    }
                }
            }
        }
        return seq;
    }
    
    vector<vector<AlignedVariant > > get_comb_vars(vector<vector<AlignedVariant> > current, vector<AlignedVariant> variants) {
        cout << "get_comb_vars" << endl;
        vector<vector<AlignedVariant> > out = current;
        BOOST_FOREACH(AlignedVariant var, variants) {
            vector<vector<AlignedVariant> > tmp;
            BOOST_FOREACH(vector<AlignedVariant> x, out) {
                vector<AlignedVariant> y1(x), y2(x);
                y1.push_back(var);
                tmp.push_back(y1);
                tmp.push_back(y2);
            }
            out.swap(tmp);
        }
        return out;
    }
    
    class LongerString {
    public:
        bool operator()(const string& riLeft, const string& riRight) const {
            return riLeft.size() > riRight.size();        }
    };
    
    vector<string> rm_dup_seqs(vector<string> seqs) {
        vector<string> tmp_seqs;
        sort(seqs.begin(), seqs.end(), LongerString());
        for (int i=0;i < seqs.size();i++) {
            string &s = seqs[i];
            bool found = false;
            for (int j=0;j < tmp_seqs.size();j++) {
                if(std::string::npos != tmp_seqs[j].find(s)) {
                    found = true;
                    break;
                }
            }
            if(!found)
                tmp_seqs.push_back(s);
        }
        return tmp_seqs;
    }
    
    vector<string> getAdditionalCombSeq(string refSeq, int left, int right, vector<AlignedVariant> variants, vector<vector<AlignedVariant> > current, Parameters params) {
        if(current.size() == 0) {
            vector<AlignedVariant> emp;
            current.push_back(emp);
        } 
        vector<vector<AlignedVariant> > potential_var_comb = get_comb_vars(current, variants);
        vector<string> haps;
        BOOST_FOREACH(vector<AlignedVariant> va, potential_var_comb) {
            haps.push_back(get_hap(refSeq, left, right, va));
        }
        return rm_dup_seqs(haps);
    }
    
    
    vector<Haplotype> getAdditionalHaps(string refSeq, string refSeqForAlign, uint32_t & leftPos, uint32_t & rightPos, vector<AlignedVariant> variants, vector<vector<AlignedVariant> > current, Parameters params) {
        vector<string> seqs = getAdditionalCombSeq(refSeqForAlign, leftPos, rightPos, variants, current, params);
        sort( seqs.begin(), seqs.end() );
        seqs.erase( unique( seqs.begin(), seqs.end() ), seqs.end() );
        vector<Haplotype> haps;
        BOOST_FOREACH(string s, seqs) {
            Haplotype h;
            h.seq = s;
            h.freq=1.0;
            h.nfreq=1.0;
            h.type=Haplotype::Normal;
            haps.push_back(h);
        }
        map<int, std::set<AlignedVariant> > var_map;
        alignHaplotypes(haps, leftPos, rightPos, var_map, params, refSeqForAlign);
        // remove duplicate reference-haplotypes of different length
        vector<Haplotype> tmp_haps;
		bool foundRef = false;
		for (size_t th=0;th<haps.size();th++) {
			const Haplotype & hap=haps[th];
			int num_indels =  hap.countIndels();
			int num_snps = hap.countSNPs();
			if (num_indels == 0 && num_snps == 0) {
				if (!foundRef) {
					tmp_haps.push_back(Haplotype(haps[th]));
					foundRef = true;
				}
			} else {
				tmp_haps.push_back(Haplotype(haps[th]));
			}
		}
        haps.swap(tmp_haps);
        BOOST_FOREACH(Haplotype h, haps) {
            cout << h.seq << endl;
        }
        return haps;
    }
    
    vector<AlignedVariant> get_hap_vars(Haplotype h, int leftPos) {
        vector<AlignedVariant> tav;
        typedef pair<int, AlignedVariant> PAV;
        BOOST_FOREACH(PAV pav, h.indels) {
            AlignedVariant &x = pav.second;
            if (!x.isRef() && !(x.isSNP() && x.getString()[3]=='D')) {
                AlignedVariant v(x.getString(), leftPos+x.getStartHap(), -1, 0);
                tav.push_back(v);
                cout << v.getString() << " " << leftPos+x.getStartHap();
            }
        }
        cout << "]" << endl;
        return tav;
    }
    
    bool getHaplotypes(vector<Haplotype> &haps, const vector<Read> & normal_reads, const vector<Read> & tumor_reads, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, vector<AlignedVariant> & close_somatic, vector<AlignedVariant>& close_germline, vector<vector<MLAlignment> > & normal_liks, vector<vector<MLAlignment> > & tumor_liks, OutputData & glfData, int index)
    {
        typedef pair<int, AlignedVariant> PAV;
        cout << "== target variant ==" << endl;
        BOOST_FOREACH(AlignedVariant v, candidateVariants.variants) {
            cout << " " << v.getStartHap() << "," << v.getString();
        }
        cout << endl;
        cout << "== close somatic ==" << endl;
        BOOST_FOREACH(AlignedVariant v, close_somatic) {
            cout << " " << v.getStartHap() << "," << v.getString();
        }
        cout << endl;
        cout << "== close germline ==" << endl;
        BOOST_FOREACH(AlignedVariant v, close_germline) {
            cout << " " << v.getStartHap() << "," << v.getString();
        }
        cout << endl;
        cout << "check vars" << endl;
        for(int i=close_germline.size()-1;i>=0;i--) {
            AlignedVariant &av = close_germline[i];
            if(av.getType() == Variant::DEL && av.size() + av.getStartHap() >= rightPos) {
                cout << "remove germline var: " << av.getString() << endl;
                close_germline.erase(close_germline.begin() + i);
            }
        }
        for(int i=close_somatic.size()-1;i>=0;i--) {
            AlignedVariant &av = close_somatic[i];
            if(av.getType() == Variant::DEL && av.size() + av.getStartHap() >= rightPos) {
                cout << "remove somatic var: " << av.getString() << endl;
                close_somatic.erase(close_somatic.begin() + i);
            }
        }
        vector<vector<AlignedVariant> > tmp_avs;
        vector<Haplotype> normal_haps = getAdditionalHaps(refSeq, refSeqForAlign, leftPos, rightPos, close_germline, tmp_avs, params);
        vector<Haplotype> tmp_haps;
        vector<vector<MLAlignment> > tmp_liks;
        if(normal_haps.size() > params.skipMaxHap) {
            cerr << "tid: " << params.tid << " pos: " << pos << " too many haplotypes [(" << normal_haps.size() << ")]" << endl;
            throw string("too many normal candidate haplotypes");
        }
        typedef map<int, AlignedVariant>::const_iterator It;
        vector<int> onHap(normal_reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality
        vector<double> normal_hap_freqs;vector<HapEstResult> normal_her;map<AlignedVariant, double> normal_vpp;lower_bound_t normal_lb;
        vector<HapFreqPair> freq_pairs;
        if(normal_haps.size() == 1) {
            cout << "no close germline vars: only ref hap in normal" << endl;
        } else {
            Haps::computeLikelihoods(normal_haps, normal_reads, normal_liks, leftPos, rightPos, onHap, params);
            EMBasic::estimate(normal_haps, normal_reads, normal_liks, normal_hap_freqs, normal_her, pos, leftPos, rightPos, glfData, index, candidateVariants, normal_lb, normal_vpp, 0.001, "all", params);
            int i = 0;
            for(int k=0;k < normal_haps.size();k++) {
                if(normal_hap_freqs[k]>0.12) { 
                    HapFreqPair hfp;hfp.index = i;i++;
                    hfp.freq = normal_hap_freqs[k];
                    freq_pairs.push_back(hfp);
                    tmp_haps.push_back(normal_haps[k]); 
                    tmp_liks.push_back(normal_liks[k]);
                }
            }
            normal_haps.swap(tmp_haps);
            normal_liks.swap(tmp_liks);
            cout << "test3" << endl;
            if(freq_pairs.size() > 1) {
                cout << "two normal haps" << endl;
                tmp_haps.clear();tmp_liks.clear();
                sort(freq_pairs.begin(), freq_pairs.end(), GreaterHapFreq());
                for(int k=0;k < freq_pairs.size();k++) {
                    HapFreqPair &hfp = freq_pairs[k];
                    cout << hfp.index << " " << hfp.freq;
                    if(k > 1) {
                        cout << " filtered!" << endl;                  
                    } else {
                        tmp_haps.push_back(normal_haps[freq_pairs[k].index]); 
                        tmp_liks.push_back(normal_liks[freq_pairs[k].index]);
                    }
                    cout << endl;
                } 
                normal_haps.swap(tmp_haps);
                normal_liks.swap(tmp_liks);
            } else if(freq_pairs.size() == 0) {
                cout << "weird region: couldnt determin normal haplotype(s)" << endl;
                throw string("any haplotypes didint meet 12% freqs in normal");
            } else {
                cout << "only one normal hap" << endl;
                cout << freq_pairs[0].index <<  " " << freq_pairs[0].freq << endl;
            }
        }
        vector<vector<AlignedVariant> > normal_vars;
        BOOST_FOREACH(Haplotype h, normal_haps) {
            normal_vars.push_back(get_hap_vars(h, leftPos));
        }
        cout << "#normal haplotype list" << endl;
        for (size_t th=0;th<normal_haps.size();th++) {
            const Haplotype & hap=normal_haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    cout << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            cout << endl;
        }
        cout << "normal hap decided" << endl;
        vector<AlignedVariant> somatic_vars(close_somatic);
        somatic_vars.push_back(candidateVariants.variants[0]);
        vector<Haplotype> merged_haps = getAdditionalHaps(refSeq, refSeqForAlign, leftPos, rightPos, somatic_vars, normal_vars, params);
        if(merged_haps.size() > params.skipMaxHap) {
            cerr << "tid: " << params.tid << " pos: " << pos << " too many haplotypes [(" << merged_haps.size() << ")]" << endl;
            throw string("too many merged candidate haplotypes");
        }
        vector<int> on_hap_tumor(tumor_reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality
        Haps::computeLikelihoods(merged_haps, tumor_reads, tumor_liks, leftPos, rightPos, on_hap_tumor, params);
        vector<double> tumor_hap_freqs;vector<HapEstResult> tumor_her;map<AlignedVariant, double> tumor_vpp;lower_bound_t tumor_lb;
        EMBasic::estimate(merged_haps, tumor_reads, tumor_liks, tumor_hap_freqs, tumor_her, pos, leftPos, rightPos, glfData, index, candidateVariants, tumor_lb, tumor_vpp, 0.001, "all", params);
        freq_pairs.clear();
        //find top tumor hap
        tmp_haps.clear();
        tmp_liks.clear();
        double max_freq=0.0;int max_index=0;bool is_normal_hap;
        int normal_only_hap_index=0;
        int changed_hap_index=0;
        vector<int> used_hap(4);
        for(int k=0;k < merged_haps.size();k++) {
            cout << "k:" << k << endl;
            is_normal_hap=false;
            for(int j=0;j < normal_haps.size();j++) {
                string &normal_seq = normal_haps[j].seq;
                string &merged_seq = merged_haps[k].seq;
                if(std::string::npos != normal_seq.find(merged_seq) || std::string::npos != merged_seq.find(normal_seq)) {                  
                    tmp_haps.push_back(merged_haps[k]); 
                    tmp_liks.push_back(tumor_liks[k]);
                    is_normal_hap=true;
                    cout << "from normal: " << merged_haps[k] << endl;
                    used_hap[k] = 1;
                }
            }
            if(!is_normal_hap) {
                if(max_freq < tumor_hap_freqs[k]) {
                    max_freq = tumor_hap_freqs[k];
                    max_index = k;
                }
                used_hap[k] = 0;
            }
        }
        used_hap[max_index] = 1;
        
        if(tmp_haps.size()==0) { throw string("no normal hap while determinig tumor hap?! something strange!"); }
        tmp_haps.push_back(merged_haps[max_index]);
        tmp_liks.push_back(tumor_liks[max_index]);
        cout << "tumor hap with max freq: " << merged_haps[max_index] << endl;
        haps.swap(tmp_haps);
        tumor_liks.swap(tmp_liks);
        if (haps.size() == 3) {
            //check which is the origin of the tumor hap
            Haplotype &th = haps[2];
            vector<AlignedVariant> th_vars = get_hap_vars(th, leftPos);
            Haplotype &h0 = haps[0];
            vector<AlignedVariant> h0_vars = get_hap_vars(h0, leftPos);
            if(th_vars.size()-1 == h0_vars.size()) {
                BOOST_FOREACH(AlignedVariant v, h0_vars) {
                    bool contains = false;
                    BOOST_FOREACH(AlignedVariant tv, th_vars) {
                        if(v.getString() == tv.getString() && v.getStartHap() == tv.getStartHap()) {contains = true;}
                    }
                    if(contains==false) {
                        //h0はtumor hapのoriginではないので入れ替える
                        swap(haps[0], haps[1]);
                        swap(tumor_liks[0], tumor_liks[1]);
                        cout << "tumor hap order change!" << endl;
                        break;
                    }
                }
            } else {
                //h0はtumor hapのoriginではないので入れ替える
                swap(haps[0], haps[1]);
                swap(tumor_liks[0], tumor_liks[1]);
                cout << "tumor hap order change!" << endl;
            }
            //h4の追加
            int h4_i = 0;
            for(int i=0;i < 4;i++) {
                if(used_hap[i] == 0) {
                    h4_i = i;
                    break;
                }
            }
            cout << "h4:" << h4_i << endl;
            haps.push_back(merged_haps[h4_i]);
            //tumor_liks.push_back(tmp_liks[h4_i]);
        }
        cout << "tumor hap decided" << endl;
        Haps::computeLikelihoods(haps, normal_reads, normal_liks, leftPos, rightPos, onHap, params);
        Haps::computeLikelihoods(haps, tumor_reads, tumor_liks, leftPos, rightPos, onHap, params);
        cout << "end Haps2" << endl;
        cout << "#normal haplotype list" << endl;
        for (size_t th=0;th<normal_haps.size();th++) {
            const Haplotype & hap=normal_haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    cout << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            cout << endl;
        }
        cout << "#tumor haplotype list" << endl;
        for (size_t th=0;th<haps.size();th++) {
            const Haplotype & hap=haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    cout << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            cout << endl;
        }
        cout << "Haps2::getHaplotypes done" << endl;
        return true;
    }
    
    bool getHaplotypesBasic(vector<Haplotype> &haps, const vector<Read> & normal_reads, const vector<Read> & tumor_reads, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, vector<vector<MLAlignment> > & normal_liks, vector<vector<MLAlignment> > & tumor_liks, OutputData & glfData, int index) {
        typedef pair<int, AlignedVariant> PAV;
        vector<vector<AlignedVariant> > tmp_avs;
        vector<AlignedVariant> somatic_var;
        typedef map<int, AlignedVariant>::const_iterator It;
        somatic_var.push_back(candidateVariants.variants[0]);
        haps = getAdditionalHaps(refSeq, refSeqForAlign, leftPos, rightPos, somatic_var, tmp_avs, params);
        //後半にtumor hapがあるかチェック
        Haplotype &fhap = haps[0];
        bool check = true;
        for (It it=fhap.indels.begin();it!=fhap.indels.end();it++) {
            if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                check = false;
            }
        }
        if(!check) {
            std::reverse(haps.begin(), haps.end());
        }
        vector<int> on_hap_tumor(tumor_reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality
        Haps::computeLikelihoods(haps, tumor_reads, tumor_liks, leftPos, rightPos, on_hap_tumor, params);
        vector<int> onHap(normal_reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality
        Haps::computeLikelihoods(haps, normal_reads, normal_liks, leftPos, rightPos, onHap, params);
        cout << "#basic haplotypes for comparison" << endl;
        for (size_t th=0;th<haps.size();th++) {
            const Haplotype & hap=haps[th];
            cout << "hap[" << th << "] " << hap.seq << endl;
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    cout << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            cout << endl;
        }
        return true;
    }
    
    bool alignHaplotypes(vector<Haplotype> & haps,  uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants, Parameters params, string refSeq)
    {
        cout << "alignHaplotypes" << endl;
        cout << "num hap: " << haps.size() << endl;
        uint32_t start=leftPos;
        uint32_t end=rightPos+1;
        
        variants.clear();
        
        int print=0;
        
        seqan::Score<int> score(-1, -460, -100,-960);
        
        Read rh1;
        rh1.pos=0;
        rh1.posStat.first=0;
        rh1.mapQual=1.0-1e-32;
        ObservationModelParameters alignParams("probabilistic");
        alignParams.pError=0.0001;
        alignParams.pMut=0.01;
        alignParams.maxLengthDel=50;
        alignParams.forceReadOnHaplotype=true;
        alignParams.bMid=0;
        //alignParams.maxLengthIndel=12;
        //alignParams.numIndels=2;
        //alignParams.indelDist="uniform";
        
        vector<Haplotype> tmp_haps;
        for (size_t h=0;h<haps.size();h++) {
            rh1.seq.seq=haps[h].seq;
            cout << h << " " << haps[h].seq << endl;
            rh1.setAllQual(1.0-1e-16);
            
            Haplotype hRef;
            uint32_t start=leftPos;
            uint32_t end=rightPos;
            
            
            hRef.append(refSeq);
            /*
             char lc = (haps[h].seq[haps[h].seq.size()-1]);
             char lcl;
             if (lc == 'T') lcl = 'A'; else if (lc == 'A') lcl = 'T'; else if (lc=='G') lcl = 'C'; else if (lc=='C') lcl = 'G';
             
             hRef.seq+= lcl;
             */
            /*
             ObservationModelFBMax om(hRef, rh1, 0, alignParams);
             */
            ObservationModelSeqAn om(hRef, rh1, 0, alignParams, score);
            haps[h].indels.clear();
            haps[h].snps.clear();
            //om.reportVariants(haps[h].indels, haps[h].snps, haps[h].align);
            //om.calcLikelihood();
            om.align();
            const MLAlignment & ml=om.getMLAlignment();
            haps[h].indels=ml.indels;
            haps[h].snps=ml.snps;
            haps[h].align=ml.align;
            haps[h].ml=ml;
            bool hasStartEndIndel = false;
            if (ml.hpos[0] == MLAlignment::LO) hasStartEndIndel = true;
            int hs = ml.hpos.size()-1;
            if (hs>0 && ml.hpos[hs] == MLAlignment::RO) hasStartEndIndel = true;
            //if (params.showCandHap) {
            //			cout << "hap " << h << endl;om.printAlignment(20);
            //			cout << string(20,' ') << haps[h].align << endl;
            //	}
            
            for (map<int, AlignedVariant>::const_iterator it=haps[h].indels.begin(); it!=haps[h].indels.end();it++) variants[it->first].insert(it->second);
            for (map<int, AlignedVariant>::const_iterator it=haps[h].snps.begin(); it!=haps[h].snps.end();it++) variants[it->first].insert(it->second);
            if (!hasStartEndIndel) {
                tmp_haps.push_back(haps[h]);
            }
            
            
        }
        
        haps.swap(tmp_haps);
        typedef map<int, AlignedVariant>::const_iterator It;
        // add REF allele as variant to each haplotype in order to compute coverage statistics
        for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
            for (size_t h=0;h<haps.size();h++) haps[h].addRefVariant(it->first);
        }
        
        if (!params.quiet) {
            for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
                cout << "aligned_var@pos " << " " << leftPos+it->first;
                BOOST_FOREACH(AlignedVariant av, it->second) {
                    cout << " " << av.getString();
                }
                cout << endl;
            }
        }
        
        
        return true;
    }
    
    
}
#endif
