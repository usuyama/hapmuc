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
#include <boost/program_options.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include "foreach.hpp"
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
#include <boost/algorithm/string.hpp>
#include "MutationCall.hpp"
#include "Haps.hpp"
#include "EMBasic.hpp"
#include "log.h"

using namespace seqan;

namespace Haps {
    
    
    bool getHaplotypes(vector<Haplotype> &haps, const vector<Read> & reads,uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign)
    {
        
        LOG(logDEBUG) << "getHaplotypes" << endl;
        uint32_t rs=(int(leftPos)>params.minReadOverlap)?(leftPos-params.minReadOverlap):0;
        uint32_t re=rightPos+params.minReadOverlap;
        HaplotypeDistribution hd(pos, refSeq, rs);
        
        // infer empirical haplotype distribution from WS alignment
        for (size_t r=0;r<reads.size();r++) {
            hd.insertRead(reads[r].getBam());
        }
        hd.setFrequencies();
        
        haps.clear();
        vector<Variant> indelVariants;
        
        /*
         if (!(candidateVariants.size()>0 && params.checkAllCIGARs==0)) {
         indelVariants=hd.getIndelVariantsAtMidPos();
         }
         
         
         // add any prespecified variants
         //indelVariants.insert(indelVariants.begin(), variants.begin(), variants.end());
         */
        if (!params.quiet) {
            LOG(logDEBUG) << "candidate_var@pos: " << pos ;
            BOOST_FOREACH(AlignedVariant v, candidateVariants.variants) {
                LOG(logDEBUG) << " " << v.getStartHap() << "," << v.getString();
            }
            LOG(logDEBUG) << endl;
        }
        
        
        // get haplotypes from empirical distribution
        
        try {
            HDIterator2 hdi(hd, params.maxHap, pos, leftPos, rightPos, params.noIndelWindow);
            
            double logNumHaps=hdi.getLogNumHaps();
            if (logNumHaps>log(params.skipMaxHap)) {
                cerr << "tid: " << params.tid << " pos: " << pos << " too many haplotypes [>exp(" << logNumHaps << ")]" << endl;
                return true;
            }
            
            //hdi.generateHapsWithIndels(haps, indels);
            vector<Haplotype> tmp_haps;
            hdi.generateHapsWithAlignedVariants(haps, candidateVariants, 0, params.changeINStoN);
            
            
            
            if (haps.size()>params.skipMaxHap || haps.size()*reads.size()>params.maxHapReadProd) {
                cerr << "tid: " << params.tid << " pos: " << pos << " too many haplotypes [>" << haps.size() << "]" << " numreads: " << reads.size() << endl;
                return true;
            }
            
            //	if (params.showHapDist) {
            /*
             LOG(logDEBUG) << endl << "Empirical distribution: " << endl;
             LOG(logDEBUG) << hdi << endl; 
             */
            //	}
            
            leftPos=hdi.start();
            rightPos=hdi.end();
            typedef map<int, AlignedVariant>::const_iterator It;
            
            LOG(logDEBUG) << "#haplotype list[debug]" << endl;
            for (size_t th=0;th<haps.size();th++) {
                const Haplotype & hap=haps[th];
                LOG(logDEBUG) << "hap[" << th << "] " << hap.seq << endl;
                for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                    if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                        LOG(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
                    }
                }
                LOG(logDEBUG) << endl;
            }
            
            
            map<int, std::set<AlignedVariant> > variants;
            alignHaplotypes(haps,pos, leftPos, rightPos, variants, params, refSeqForAlign);
            
            LOG(logDEBUG) << "#haplotype list[debug2]" << endl;
            for (size_t th=0;th<haps.size();th++) {
                const Haplotype & hap=haps[th];
                LOG(logDEBUG) << "hap[" << th << "] " << hap.seq << endl;
                for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                    if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                        LOG(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
                    }
                }
                LOG(logDEBUG) << endl;
            }
            
            // remove duplicate reference-haplotypes of different length
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
            
            // if (params.showCandHap) {
            for (size_t i=0;i<haps.size();i++) {
                // LOG(logDEBUG) << "PRE FILTER hdi[" << i << "]:" << haps[i] << endl;
            }
            //}
            
            
            typedef map<int, AlignedVariant>::const_iterator It;
            haps.swap(tmp_haps);
            
            int nh=0;
            //if (params.showCandHap) {
			for (size_t i=0;i<haps.size();i++) {
                //		LOG(logDEBUG) << "POSTFILTER hdi[" << nh++ << "]:" << haps[i] << endl;
			}
            //}
        }
        catch (string s) {
            if (s=="Blocks are not consecutive.") {
                cerr << "tid: " << params.tid <<  "pos: " << pos << s << endl;
                //return true;
                throw string("hapblock");
            } else {
                throw string(s);
            }
        }
        return false;
    }
    
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
    
    void filterReads(vector<Read> & mergedReads, vector<vector<MLAlignment> > & mergedLiks, int & normal_count, int & tumor_count) {
        LOG(logDEBUG) << "filterReads" << endl;
        int numhap = (int)mergedLiks.size();
        if(numhap < 2) return;
        for(int i = 0;i < mergedReads.size();) {
            LOG(logDEBUG) << "read " << i << " " << mergedReads[i].seq_name;
            double max = -HUGE_VAL;
            double min = 0.0;
            int max_overlap = 0;
            for(int k = 0;k < numhap;k++) {
                double ll = mergedLiks[k][i].ll	;
                if(max < ll) max = ll;
                if(min > ll) min = ll;
                int overlap = mergedLiks[k][i].lastBase - mergedLiks[k][i].firstBase;
                if(max_overlap < overlap) max_overlap = overlap;
            }
            LOG(logDEBUG) << " " << max << " " << min << " " << max_overlap;
            if(max_overlap < 25 || (max - min) < 0.5) {
                LOG(logDEBUG) << " filter reads!";
                if(i < normal_count)
                    normal_count--;
                else                
                    tumor_count--;
                vector<Read>::iterator it = mergedReads.begin();
                advance(it, i);
                mergedReads.erase(it);
                for(int k = 0;k < numhap;k++) {
                    vector<MLAlignment>::iterator it = mergedLiks[k].begin();
                    advance(it, i);
                    mergedLiks[k].erase(it);
                }
            } else {
                i++;
            }
            LOG(logDEBUG) << endl;
        }
        LOG(logDEBUG) << mergedReads.size() << " " << normal_count << " " << tumor_count << endl;
    }


    void selectHaplotypesAndReads(vector<Haplotype> & haps, vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, int & normal_count, int & tumor_count, Parameters params) {
        LOG(logDEBUG) << "############ select haplotypes ##############" << endl;
        //filter duplicate haps
        int i = 0;
        while(i < haps.size() - 1) {
            string seq = haps[i].seq;
            for(int j = i + 1;j < haps.size();j++) {
                if(seq == haps[j].seq) {
                    //      LOG(logDEBUG) << i << " " << j << endl;
                    haps.erase(haps.begin()+j);
                    break;
                }
            }
            i++;
        }
        //filter haplotypes with variants on edge
        
        i = 0;
        typedef map<int, AlignedVariant>::const_iterator It;
        
        int margin =  params.width - 20;
        while(i < haps.size()) {
            Haplotype hap = haps[i];
            bool flag = false;
            //LOG(logDEBUG) << "i, margin, len " << i << ", " << margin << ", " << hap.seq.length();
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    //   LOG(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
                    if(it->first < margin || it->first > hap.seq.length() - margin) {
                        flag = true;
                    }
                }
            }
            // LOG(logDEBUG) << endl;
            if(flag) {
                //   LOG(logDEBUG) << "filter " << i << endl;
                haps.erase(haps.begin()+i);
            } else {
                i++;
            }
            
        }     
        
        vector<int> onHap(reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality
        vector<int> filtered(haps.size(), 0);
        map<pair<int, AlignedVariant>, VariantCoverage> varCoverage;
        computeLikelihoods(haps, reads, liks, leftPos, rightPos, onHap, params);
        filterHaplotypes(haps, reads, liks, filtered, varCoverage, params.filterHaplotypes, params);
        vector<Haplotype> tmp_haps;
        vector<vector<MLAlignment> > tmp_liks;
        for(int k=0;k < haps.size();k++) {
            if(!filtered[k]) { 
                tmp_haps.push_back(haps[k]); 
                tmp_liks.push_back(liks[k]);
            }
        }
        haps = tmp_haps;
        liks = tmp_liks;
        
        LOG(logDEBUG) << "############ vb for filter haplotypes ##############" << endl;
        vector<double> hapFreqs;
        map <int, vector<tuple<AlignedVariant, double, double> > > posteriors;
        lower_bound_t lb;
        vector<HapEstResult> her;
        map<AlignedVariant, double> vpp;
        vector<double> filter_th;
        filter_th.push_back(0.01);
        filter_th.push_back(0.02);
        filter_th.push_back(0.01);
        filter_th.push_back(0.01);
        filter_th.push_back(0.01);
        for(int i=0;i < 1;i=i++) {
            LOG(logDEBUG) << "#### filter haps[" << i  << " ] ###" << endl;
            EMBasic::estimate(haps, reads, liks, hapFreqs, her, pos, leftPos, rightPos, candidateVariants, lb, vpp, 0.1, "all", params);
            vector<Haplotype> tmp_haps;
            vector<vector<MLAlignment> > tmp_liks;
            int passed = 1;
            for(int k=0;k < haps.size();k++) {
                if(hapFreqs[k]>filter_th[i]) { 
                    tmp_haps.push_back(haps[k]); 
                    tmp_liks.push_back(liks[k]);
                } else {
                    passed = 0;
                }
            }
            haps = tmp_haps;
            liks = tmp_liks;
            if(haps.size() < 2) return;
            filterReads(reads, liks, normal_count, tumor_count);
            if(haps.size() < 3 && i > 2) return;
        }
        LOG(logDEBUG) << endl;
        LOG(logDEBUG) << "############ select haplotypes done ##############" << endl;
    }

    
    bool alignHaplotypes(vector<Haplotype> & haps,  uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants, Parameters params, string refSeq)
    {
        LOG(logDEBUG) << "alignHaplotypes" << endl;
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
            //			LOG(logDEBUG) << "hap " << h << endl;om.printAlignment(20);
            //			LOG(logDEBUG) << string(20,' ') << haps[h].align << endl;
            //	}
            
            for (map<int, AlignedVariant>::const_iterator it=haps[h].indels.begin(); it!=haps[h].indels.end();it++) variants[it->first].insert(it->second);
            for (map<int, AlignedVariant>::const_iterator it=haps[h].snps.begin(); it!=haps[h].snps.end();it++) variants[it->first].insert(it->second);
            if (!hasStartEndIndel) {
                tmp_haps.push_back(haps[h]);
            }
            
            
        }
        
        haps.swap(tmp_haps);
        typedef map<int, AlignedVariant>::const_iterator It;
        /*   LOG(logDEBUG) << "#haplotype list[debug4]" << endl;
         for (size_t th=0;th<haps.size();th++) {
         const Haplotype & hap=haps[th];
         LOG(logDEBUG) << "hap[" << th << "] " << hap.seq << endl;
         for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
         if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
         LOG(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
         }
         }
         LOG(logDEBUG) << endl;
         }*/
        
        
        // add REF allele as variant to each haplotype in order to compute coverage statistics
        for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
            for (size_t h=0;h<haps.size();h++) haps[h].addRefVariant(it->first);
        }
        
        /*  LOG(logDEBUG) << "#haplotype list[debug5]" << endl;
         for (size_t th=0;th<haps.size();th++) {
         const Haplotype & hap=haps[th];
         LOG(logDEBUG) << "hap[" << th << "] " << hap.seq << endl;
         for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
         if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
         LOG(logDEBUG) << "[" << it->second.getString() << " " << (it->first) << "]";
         }
         }
         LOG(logDEBUG) << endl;
         }*/
        
        if (!params.quiet) {
            for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
                LOG(logDEBUG) << "aligned_var@pos " << pos << " " << leftPos+it->first;
                BOOST_FOREACH(AlignedVariant av, it->second) {
                    LOG(logDEBUG) << " " << av.getString();
                }
                LOG(logDEBUG) << endl;
            }
        }
        
        
        return true;
    }


    
#define FILTERHAPS
#ifdef FILTERHAPS
    
    void filterHaplotypes(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks,  vector<int> & filtered, map<pair<int, AlignedVariant>, VariantCoverage> & varCoverage, bool doFilter, Parameters params)
    {
        
        const int debugfh = 0;
        int numFiltered = 0;
        int numHaps = int(haps.size());
        filtered = vector<int>(haps.size(),0);
        
        varCoverage.clear();
        
        
        typedef pair<int, AlignedVariant> PAV;
        map<PAV, vector< std::set<int> > > hVarCoverage; // coverage of each variant per haplotype, so that eventually only coverage of reads that were not filtered out are reported
        
        for (int h=0;h<int(haps.size());h++) {
            // check all reads aligned to this haplotype and select the ones that are not off-haplotype and aligned without indels and at most two high-quality mismatches
            std::set<int> selReads;
            for (size_t r=0;r<reads.size();r++) {
                int sel = 0;
                if (!liks[h][r].offHapHMQ && liks[h][r].numIndels == 0) { // && liks[h][r].numMismatch<=3) {
                    selReads.insert(int(r));
                    sel = 1;
                }
                if (debugfh) LOG(logDEBUG) << "sel: " << "h: " << h << " " << bam1_qname(reads[r].getBam()) << " mpos: " << reads[r].matePos << " selected: " << sel << endl;
                
            }
            
            // check each variant in the haplotype
            
            bool allCovered = true; // all variants in haplotype should be covered by at least one read.
            for (map<int, AlignedVariant>::const_iterator it = haps[h].indels.begin();it!= haps[h].indels.end();it++) {
                const AlignedVariant & av = it->second;
                
                PAV pav(it->first, av);
                
                map<PAV,vector<std::set<int> > >::iterator pit  = hVarCoverage.find(pav);
                if (pit == hVarCoverage.end()) {
                    hVarCoverage[pav] = vector< std::set<int> >(haps.size()*2);
                }
                
                if (av.getType() == Variant::INS || av.getType() == Variant::DEL) {
                    int left = av.getLeftFlankRead() - params.obsParams.padCover;     // readFlankLeft is the first left unique base flanking the indel
                    int right = av.getRightFlankRead() + params.obsParams.padCover;
                    int leftV = av.getLeftFlankRead();
                    int rightV = av.getRightFlankRead();
                    
                    
                    int len = right-left+1;
                    bool covered = false;
                    int numdelcovered = 0;
                    //LOG(logDEBUG) << "left: " << left << " right: " << right << " len: " << len << endl;
                    if (av.getType() == Variant::DEL) {
                        // see if there is at least one read spanning the interval with at most one mismatch
                        BOOST_FOREACH(int r, selReads) {
                            std::set <int> c;
                            int strand = 0;
                            if (reads[r].isUnmapped()) {
                                if (!reads[r].mateIsReverse()) strand = 1;
                            } else {
                                if (reads[r].isReverse()) strand = 1;
                            }
                            
                            int nmm = 0;
                            for (int b=0;b<=int(liks[h][r].hpos.size());b++) {
                                int hb = liks[h][r].hpos[b];
                                if (hb>=left && hb<=right) {
                                    c.insert(hb);
                                    if ( haps[h].seq[hb]!='N' && haps[h].seq[hb]!=reads[r].seq.seq[b]) nmm++;
                                }
                                
                            }
                            int cov = 0;
                            if (int(c.size())>= len && nmm<=params.obsParams.maxMismatch) {
                                cov = 1;
                                numdelcovered++;
                                hVarCoverage[pav][h+strand*numHaps].insert(r);
                            }
                            //LOG(logDEBUG) << "RC" << bam1_qname(reads[r].getBam()) << " cov: " << cov << endl;
                        }
                        
                        if (numdelcovered>=1) {
                            covered = true;
                        }
                    } else if (av.getType() == Variant::INS) {
                        // see if all bases in the haplotype from left to right are covered by at least one read that matches the haplotype without indels and at most 2 mismatches
                        vector<int> hapBaseCovered(len, 0);
                        vector<int> thisReadCovered(len, 0);
                        int lenins = av.getSeq().size();
                        
                        BOOST_FOREACH(int r, selReads) {
                            for (int x=0;x<len;x++) thisReadCovered[x]=0;
                            int nmm = 0;
                            std::set <int> c;
                            
                            // determine read strand
                            int strand = 0;
                            if (reads[r].isUnmapped()) {
                                if (!reads[r].mateIsReverse()) strand = 1;
                            } else {
                                if (reads[r].isReverse()) strand = 1;
                            }
                            
                            for (int b=0;b<=int(liks[h][r].hpos.size());b++) {
                                int hb = liks[h][r].hpos[b];
                                if (hb>=left && hb<=right) {
                                    // covered even if there is a mismatch
                                    thisReadCovered[hb-left]+=1;
                                    c.insert(hb);
                                    // count number of mismatches
                                    if (haps[h].seq[hb]!=reads[r].seq.seq[b]) nmm++;
                                }
                            }
                            
                            bool thisread_covered = false;
                            // for insertion <= 10 bp, whole insertion+padCover must be covered with at most one error by at least one read (ie not just covered by multiple reads together)
                            if ( (lenins>10 && nmm<=params.obsParams.maxMismatch) || (lenins<=10 && int(c.size())>=len && nmm<=params.obsParams.maxMismatch)) {
                                thisread_covered = true;
                                for (size_t x=0;x<thisReadCovered.size();x++) {
                                    hapBaseCovered[x] += thisReadCovered[x];
                                    if (thisReadCovered[x]==0) {
                                        thisread_covered = false;
                                    }
                                    if (debugfh) LOG(logDEBUG) << " " << hapBaseCovered[x];
                                }
                                if (thisread_covered) {
                                    hVarCoverage[pav][h+strand*numHaps].insert(r);
                                }
                            }
                            
                            
                            if (0) LOG(logDEBUG) << " hap " << h << " var: " << av.getString() << " len: " << len << " " << bam1_qname(reads[r].getBam()) << " nmm: " << nmm << " c.size(): " << c.size() << " mpos: " << reads[r].matePos << " covered: " << thisread_covered << endl;
                            if (thisread_covered) covered=true;
                            
                            
                        }
                    }
                    
                    if (!covered) {
                        allCovered = false;
                        break;
                    }
                    if (debugfh) LOG(logDEBUG) << "hap" << h << " var: " << av.getString() << " COVERED:" << covered << endl;
                    
                }
                
            }
            if (doFilter) {
                if (!allCovered) {
                    numFiltered++;
                    filtered[h]=1;
                }
            }
            //LOG(logDEBUG) << "Haplotype[" << h << "]: filtered " << filtered[h] << endl;
            
        }
        LOG(logDEBUG) << "Filtered " << numFiltered << " haplotypes." << endl;
        // determine coverage of each variant
        for (map<PAV, vector <std::set<int> > >::const_iterator it = hVarCoverage.begin();it != hVarCoverage.end(); it++) {
            const PAV & pav = it->first;
            std::set<int> rf, rr; // forward and reverse strand reads
            for (int h=0;h<numHaps;h++) if (filtered[h]!=1) {
                rf.insert(hVarCoverage[pav][h].begin(), hVarCoverage[pav][h].end());
                rr.insert(hVarCoverage[pav][h+numHaps].begin(), hVarCoverage[pav][h+numHaps].end());
            }
            varCoverage[pav]=VariantCoverage(int(rf.size()), int(rr.size()));
        }
        
        
    }
#endif

}
#endif