/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
//#include <stdlib.h>
#include <iostream>
//#include <iomanip>
//#include <string>
//#include <sstream>
//#include <fstream>
//#include <set>
//#include <algorithm>
//#include <boost/program_options.hpp>
//#include <boost/assign.hpp>
//#include <boost/math/special_functions/digamma.hpp>
//#include <seqan/align.h>
//#include <seqan/graph_align.h>
//#include "foreach.hpp"
//#include "bam.h"
#include "HapMuC.hpp"
//#include "Haplotype.hpp"
//#include "HaplotypeDistribution.hpp"
//#include "ObservationModelFB.hpp"
//#include "Utils.hpp"
//#include "faidx.h"
//#include "ObservationModelSeqAn.hpp"
//#include "VariantFile.hpp"
//#include <ext/hash_map>
//#include <exception>
//#include <math.h>
//#include <sys/stat.h>
//#include <sys/types.h>
//#include <boost/algorithm/string.hpp>
#include "Haps2.hpp"
#include "EMfor2.hpp"
#include "MutationCall.hpp"
#include "EMBasic.hpp"
#include "MutationModel.hpp"
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <float.h>
#include <log.h>

namespace MutationCall
{
    string itos(int num) {
        std::string str = boost::lexical_cast<std::string>(num);
        return(str);
    }

    void output_liktable(int index, string id, const vector<double> & liks, Parameters params) {
        //hap size=4
#ifndef LOGDEBUG
        return;
#endif
        std::ofstream ofs((params.fileName+"logs/"+itos(index)+"."+id+".liks").c_str(), std::ios::out);
        int nr = liks.size() / 4;
        for(int i=0;i<nr;i++) {
            ofs << liks[i*4];
            for(int j=1;j<4;j++)
                ofs << "\t" << liks[i*4+j];
            ofs << endl;
        }
    }

    void output_alignments(int index, string id, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, Parameters params) {
        if (!params.showHapAlignments) return;
        std::ofstream ofs((params.fileName+"logs/"+itos(index)+"."+id+".alignments").c_str(), std::ios::out);
        int nr = reads.size();
        int nh = liks.size();
        typedef map<int, AlignedVariant>::const_iterator It;
        for(size_t i=0;i<nr;i++) {
            ofs << reads[i].seq_name;
            for(size_t j=0;j<nh;j++)  ofs << "\t" << liks[j][i].ll;
            ofs << endl;
            for(size_t j=0;j<nh;j++)  ofs << "\t" << liks[j][i].align << endl;
            for(size_t j=0;j<nh;j++) {
                const MLAlignment &tmp = liks[j][i];
                ofs << "#" << j << "\t" << "relPos: " << liks[j][i].relPos << " offHap: " << (int)liks[j][i].offHap;
                ofs << " firstBase:" << tmp.firstBase << " lastBase:" << tmp.lastBase << " numMismatch:" << tmp.numMismatch << " numIndels:" << tmp.numIndels << " ";
                for (map<int, AlignedVariant>::const_iterator it=tmp.indels.begin(); it!=tmp.indels.end();it++) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
                ofs << " ";
                for (map<int, AlignedVariant>::const_iterator it=tmp.snps.begin(); it!=tmp.snps.end();it++) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
                ofs << endl;
            }
        }
    }

    void filter_reads(vector<Read> &reads, vector<double> &liks) {
        if (liks.size() % 4 != 0) throw(string("haps size should be 4!"));
        // TODO: can be more efficient.
        LOG(logDEBUG) << "filter_reads: " << endl;
        int count = 0;
        //h1とh3で差がでないreadを消す。
        for(int i=0;i<reads.size();i++) {
            double max=-DBL_MAX, min=DBL_MAX;
            for(int j=0;j<2;j++) {
                if(liks[i*4+j*2]>max) max = liks[i*4+j*2];
                if(liks[i*4+j*2]<min) min = liks[i*4+j*2];
            }
            if((max-min)<0.0001) {
                // erase
            } else {
                if (count != i) {
                    reads[count] = reads[i];
                    std::copy(liks.begin()+i*4, liks.begin()+i*4+4, liks.begin()+count*4);
                }
                count++;
            }
        }
        reads.erase(reads.begin()+count, reads.end());
        liks.erase(liks.begin()+count*4, liks.end());
        LOG(logDEBUG) << "filter_reads done" << endl;
    }

    void output_mm(vector<Haplotype> &haps, string fname, uint32_t leftPos, uint32_t rightPos, double bf, lower_bound_t& lb, vector<double> normalHapFreqs, vector<double> tumorHapFreqs) {
#ifndef LOGDEBUG
        return;
#endif
        typedef map<int, AlignedVariant>::const_iterator It;
        std::ofstream ofs(fname.c_str(), std::ios::out | std::ios::app);
        ofs << "-------- window[" << leftPos << "-" << rightPos <<"] ---------------" << endl;
        ofs << "#tumor " << lb.lower_bound << endl;
        for (size_t th=0;th<tumorHapFreqs.size();th++) {
            const Haplotype & hap=haps[th];
            ofs << "hap[" << th << "] ";
            ofs << tumorHapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs << endl;
        }
        ofs << "#normal " << endl;
        for (size_t th=0;th<normalHapFreqs.size();th++) {
            const Haplotype & hap=haps[th];
            ofs << "hap[" << th << "] ";
            ofs << normalHapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs << endl;
        }
    }


    void print_freqs(vector<Haplotype> &haps, string fname, uint32_t leftPos, uint32_t pos, uint32_t rightPos, const vector<Read> & normal_reads, const vector<Read> & tumor_reads, const AlignedCandidates & candidateVariants, const vector<vector<MLAlignment> > & normal_liks, const vector<vector<MLAlignment> > & tumor_liks, Parameters params) {
#ifndef LOGDEBUG
        return;
#endif
        vector<double> tumor_hap_freqs, normal_hap_freqs;
        vector<HapEstResult> tumor_her, normal_her;
        map<AlignedVariant, double> normal_vpp, tumor_vpp;
        lower_bound_t normal_lb, tumor_lb;
        EMBasic::estimate(haps, tumor_reads, tumor_liks, tumor_hap_freqs, tumor_her, pos, leftPos, rightPos, candidateVariants, tumor_lb, tumor_vpp, 1.0, "all", params);
        EMBasic::estimate(haps, normal_reads, normal_liks, normal_hap_freqs, normal_her, pos, leftPos, rightPos, candidateVariants, normal_lb, normal_vpp, 1.0, "all", params);
        std::ofstream ofs(fname.c_str(), std::ios::out | std::ios::app);
        ofs << pos << "\t";
        if(haps.size() == 2) {
            ofs << tumor_hap_freqs[0] << "\t" << tumor_hap_freqs[1] << "\t-\t";
            ofs << normal_hap_freqs[0] << "\t" << normal_hap_freqs[1] << "\t-" << endl;
        } else {
            ofs << tumor_hap_freqs[0] << "\t" << tumor_hap_freqs[1] << "\t" << tumor_hap_freqs[2] << "\t";
            ofs << normal_hap_freqs[0] << "\t" << normal_hap_freqs[1] << "\t" << normal_hap_freqs[2] << endl;
        }
        typedef map<int, AlignedVariant>::const_iterator It;
        std::ofstream ofs2((params.fileName+".general_free_haps.txt").c_str(), std::ios::out | std::ios::app);
        ofs2 << "-------- window[" << leftPos << "-" << rightPos <<"] ---------------" << endl;
        ofs2 << "#tumor " << tumor_lb.lower_bound << endl;
        for (size_t th=0;th<haps.size();th++) {
            const Haplotype & hap=haps[th];
            ofs2 << "hap[" << th << "] ";
            ofs2 << tumor_hap_freqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs2 << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs2 << endl;
        }
        ofs2 << "#normal " << normal_lb.lower_bound << endl;
        for (size_t th=0;th<haps.size();th++) {
            const Haplotype & hap=haps[th];
            ofs2 << "hap[" << th << "] ";
            ofs2 << normal_hap_freqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs2 << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs2 << endl;
        }
    }

    void computeBayesFactors(const vector<Read> & normalReads, const vector<Read> & tumorReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, OutputData & oData, Parameters params, string refSeq, string refSeqForAlign, vector<AlignedVariant> & close_somatic, vector<AlignedVariant>& close_germline, int index)
    {
        LOG(logDEBUG) << "computeBayesFactor: " << index << endl;
        vector<Haplotype> haps;
        vector<vector<MLAlignment> > liks;
        vector<HapPairLik> likPairs;

        LOG(logDEBUG) << "check for long indel" << endl;
        if(candidateVariants.variants[0].getString().length() > 50) { // TODO: magic number
            throw string("too long indel." + candidateVariants.variants[0].getString());
        }

        typedef map<int, AlignedVariant>::const_iterator It;

        // NOTE leftPos will be the left position of the reference sequence each haplotype will be aligned to
        vector<vector<MLAlignment> > normal_liks, normal_liks2;
        vector<vector<MLAlignment> > tumor_liks, tumor_liks2;
        int refSeqPos=leftPos;
        vector<double> normalHapFreqs, non_normalHapFreqs, tumorHapFreqs, non_tumorHapFreqs;
        vector<HapEstResult> normal_her, non_normal_her;
        vector<HapEstResult> tumor_her, non_tumor_her;
        vector<HapEstResult> tumor_her_hap2, normal_her_hap2, merged_her_hap2;

        lower_bound_t mm_lb, nmm_lb;
        map<AlignedVariant, double> normal_vpp, non_normal_vpp;
        map<AlignedVariant, double> tumor_vpp, non_tumor_vpp;
        int normal_count = (int)normalReads.size();
        int tumor_count = (int)tumorReads.size();
        double bf2;
        LOG(logDEBUG) <<  "centerPos=" << candidateVariants.centerPos << endl;
        candidateVariants.printAll();
        double hap2_bf = -100;
        bool hap2_flag = true, hap3_flag = true;
        bool haps_size_error = false;
        string closest_germline = "-";
        int distance = -1;
        try {
            hap2_bf = calc_hap2_bf_with_hap3(normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, params, refSeq, refSeqForAlign, 1.0, 1.0, 1.0, normal_her_hap2, tumor_her_hap2, merged_her_hap2, index);
        } catch (string s) {
            LOG(logERROR) << "error_" << s << endl;
            hap2_flag = false;
        } catch (std::exception& e) {
            string message = string("error_exception_").append(e.what());
            hap2_flag = false;
        }
        try {
            if (close_germline.size() != 0) {
                LOG(logINFO) << string(50, '-') << endl;
                LOG(logINFO) << "### calculate the Bayes factor with SNP information nearby:" << endl;
                //close_somaticは空
                //close_germlineは、targetに最も近いものを採用する
                LOG(logWARNING) << "choose the closest SNP" << endl;
                vector<AlignedVariant> new_close_somatic;
                vector<AlignedVariant> new_close_germline;
                int target_pos = candidateVariants.variants[0].getStartHap();
                LOG(logDEBUG) << "target_pos " << target_pos << endl;
                int min_dis= 100000, min_index=-1;
                for(int i=0;i<close_germline.size();i++) {
                    if(close_germline[i].isSNP()) {
                        int dist = abs(close_germline[i].getStartHap() - target_pos);
                        if(min_dis > dist) {
                            min_index = i;min_dis = dist;
                        }
                    }
                }
                if (min_index == -1) {
                    LOG(logWARNING) << "no germline SNP (!= indel); skip!" << endl;
                    hap3_flag = false;
                } else {
                    LOG(logDEBUG) << "min " << min_dis << " " << min_index << endl;
                    stringstream ss;
                    ss << close_germline[min_index].getStartHap() << ":" << close_germline[min_index].getString();
                    closest_germline = ss.str();
                    LOG(logWARNING) << "the closest SNP is " << closest_germline << endl;
                    distance = min_dis;
                    new_close_germline.push_back(close_germline[min_index]);
                    bool flag = Haps2::getHaplotypes(haps, normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, params, refSeq, refSeqForAlign, new_close_somatic, new_close_germline, normal_liks, tumor_liks);
                    if(haps.size()==4) {
                        LOG(logDEBUG) << "************ hap4 **************" << endl << "haps: " << haps.size() << endl;
                        int nrt = tumorReads.size(), nrn=normalReads.size();
                        vector<double> rlt(nrt*4), rln(nrn*4);
                        int idx=0;
                        for (size_t r=0;r<nrt;r++) {
                            for (size_t h=0;h<4;h++) { rlt[idx] = tumor_liks[h][r].ll;idx++; }
                        }
                        idx=0;
                        for (size_t r=0;r<nrn;r++) {
                            for (size_t h=0;h<4;h++) { rln[idx] = normal_liks[h][r].ll;idx++; }
                        }

                        output_alignments(index, "hap3_tumor", tumorReads, tumor_liks, params);
                        output_alignments(index, "hap3_normal", normalReads, normal_liks, params);

                        vector<Read> f_normal_reads, f_tumor_reads;
                        f_normal_reads.insert(f_normal_reads.end(), normalReads.begin(), normalReads.end());
                        f_tumor_reads.insert(f_tumor_reads.end(), tumorReads.begin(), tumorReads.end());
                        LOG(logDEBUG) << "before " << f_normal_reads.size() << " " << f_tumor_reads.size() << " " << rln.size()  << " " << rlt.size() << endl;
                        LOG(logINFO) << "filter uninformative reads. remaining ";
                        filter_reads(f_normal_reads, rln);
                        filter_reads(f_tumor_reads, rlt);
                        LOGP(logINFO) << "tumor reads: " << f_tumor_reads.size() << " of " << nrt << ", normal reads: " << f_normal_reads.size() << " of " << nrn << endl;
                        LOG(logDEBUG) << " " << rln.size() << " " << rlt.size() << endl;

                        output_liktable(index, "hap3_tumor", rlt, params);
                        output_liktable(index, "hap3_normal", rln, params);
                        LOG(logINFO) << "calculating the marginal likelihood based on the MUTATION model...";
                        MutationModel::estimate(haps, f_tumor_reads, f_normal_reads, rlt, rln, tumorHapFreqs, normalHapFreqs, tumor_her, normal_her, pos, leftPos, rightPos, candidateVariants, mm_lb, params.hap3_params, "mutation", params.fileName+"logs/"+itos(index)+".hap3.mutation");
                        LOG(logINFO) << "calculating the marginal likelihood based on the ERROR model...";
                        MutationModel::estimate(haps, f_tumor_reads, f_normal_reads, rlt, rln, non_tumorHapFreqs, non_normalHapFreqs, non_tumor_her, non_normal_her, pos, leftPos, rightPos, candidateVariants, nmm_lb, params.hap3_params, "error",params.fileName+"logs/"+itos(index)+".hap3.err");
                        bf2 = mm_lb.lower_bound - nmm_lb.lower_bound;
                        LOG(logDEBUG) << "**** hap4 done ***" << endl;
                    } else {
                        hap3_flag = false;
                        haps_size_error = true;
                        LOG(logERROR) << "error_" << "haps.size() != 4 in hap3; the region is a bit suspicious." << endl;
                    }
                }
            } else {
                hap3_flag = false;
            }
        } catch (string s) {
            LOG(logERROR) << "error_" << s << endl;
            hap3_flag = false;
        } catch (std::exception& e) {
            string message = string("error_exception_").append(e.what());
            LOG(logERROR) << "error_" << message << endl;
            hap3_flag = false;
        }
        if(tumor_her.empty()) {
            tumor_her.swap(tumor_her_hap2);
            normal_her.swap(normal_her_hap2);
        }
        LOG(logDEBUG) << "output" << endl;
        if(tumor_her.empty()) {
            const VariantInfo &v = candidateVariants.variants[0].info;
            int pos = v.start;
            int end = v.end;
            OutputData::Line line(oData);
            line.set("chr", params.tid);
            line.set("start", pos);
            line.set("end", end);
            line.set("ref", v.ref);
            line.set("obs", v.obs);
            line.set("ref_count_tumor", v.ref_count_tumor);
            line.set("obs_count_tumor", v.obs_count_tumor);
            line.set("ref_count_normal", v.ref_count_normal);
            line.set("obs_count_normal", v.obs_count_normal);
            line.set("missrate_tumor", v.missrate_tumor);
            line.set("strandrate_tumor", v.strandrate_tumor);
            line.set("missrate_normal", v.missrate_normal);
            line.set("strandrate_normal", v.strandrate_normal);
            line.set("fisher", v.fisher_score);
            line.set("ref_bq_tumor", v.ref_bq_tumor);
            line.set("obs_bq_tumor", v.obs_bq_tumor);
            line.set("ref_bq_normal", v.ref_bq_normal);
            line.set("obs_bq_normal", v.obs_bq_normal);
            line.set("bf_without_snp", "-");
            line.set("bf_with_snp", "-");
            line.set("germline_snp_nearby", "-");
            line.set("distance", "-");
            oData.output(line);
        } else {
            for(int i=0;i<tumor_her.size();i++) {
                LOG(logDEBUG) << i << endl;
                HapEstResult &n_result = normal_her[i];
                HapEstResult &t_result = tumor_her[i];
                string chr = params.tid;
                VariantInfo &v = t_result.info;
                int pos = v.start;
                int end = v.end;
                OutputData::Line line(oData);
                line.set("chr", params.tid);
                line.set("start", pos);
                line.set("end", end);
                line.set("ref", v.ref);
                line.set("obs", v.obs);
                line.set("ref_count_tumor", v.ref_count_tumor);
                line.set("obs_count_tumor", v.obs_count_tumor);
                line.set("ref_count_normal", v.ref_count_normal);
                line.set("obs_count_normal", v.obs_count_normal);
                line.set("missrate_tumor", v.missrate_tumor);
                line.set("strandrate_tumor", v.strandrate_tumor);
                line.set("missrate_normal", v.missrate_normal);
                line.set("strandrate_normal", v.strandrate_normal);
                line.set("fisher", v.fisher_score);
                line.set("ref_bq_tumor", v.ref_bq_tumor);
                line.set("obs_bq_tumor", v.obs_bq_tumor);
                line.set("ref_bq_normal", v.ref_bq_normal);
                line.set("obs_bq_normal", v.obs_bq_normal);
                if(hap2_flag) {
                    line.set("bf_without_snp", hap2_bf);
                } else {
                    line.set("bf_with_snp", "-");
                }
                if(hap3_flag){
                    line.set("bf_with_snp", bf2);
                } else {
                    if(haps_size_error) {
                        line.set("bf_with_snp", "-");
                    } else {
                        line.set("bf_with_snp", "-");
                    }
                }
                if(closest_germline == "-") {
                    line.set("germline_snp_nearby", "-");
                    line.set("distance", "-");
                } else {
                    line.set("germline_snp_nearby", closest_germline);
                    line.set("distance", distance);
                }
                oData.output(line);
            }
        }
    }

    double calc_bf(const vector<Haplotype> &haps, const vector<vector<MLAlignment> > &normal_liks, const vector<vector<MLAlignment> > tumor_liks, const vector<Read> & normalReads, const vector<Read> & tumorReads, const vector<Read> & mergedReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, double a0, double b0, double free_a0, vector<HapEstResult> &normal_her, vector<HapEstResult> &tumor_her, vector<HapEstResult> &merged_her) {
        // algorithm of undergrad-thesis
        LOG(logDEBUG) << "********** calc_bf ***********" << endl;
        vector<vector<MLAlignment> > liks;
        liks.insert(liks.end(), normal_liks.begin(), normal_liks.end());
        for(int i = 0;i < haps.size();i++) {
            liks[i].insert(liks[i].end(), tumor_liks[i].begin(), tumor_liks[i].end());
        }
        vector<double> mergedHapFreqs, normalHapFreqs, tumorHapFreqs;
        lower_bound_t normal_lb;
        lower_bound_t tumor_lb;
        lower_bound_t merged_lb;
        map<AlignedVariant, double> merged_vpp;
        map<AlignedVariant, double> normal_vpp;
        map<AlignedVariant, double> tumor_vpp;
        EMBasic::estimate(haps, normalReads, normal_liks, normalHapFreqs, normal_her, pos, leftPos, rightPos, candidateVariants, normal_lb, normal_vpp,  1.0, "all", params);
        EMBasic::estimate(haps, tumorReads, tumor_liks, tumorHapFreqs, tumor_her, pos, leftPos, rightPos, candidateVariants, tumor_lb, tumor_vpp, 1.0, "all", params);
        EMBasic::estimate(haps, mergedReads, liks, mergedHapFreqs, merged_her, pos, leftPos, rightPos, candidateVariants, merged_lb, merged_vpp,  1.0, "all", params);
        double bf = normal_lb.lower_bound + tumor_lb.lower_bound - merged_lb.lower_bound;
        outputLowerBounds(haps, (params.fileName+".basic_haplotypes.txt"), leftPos, rightPos, bf, normal_lb, tumor_lb, merged_lb, normalHapFreqs, tumorHapFreqs, mergedHapFreqs);
        return bf;
    }

    double calc_hap2_bf_with_hap3(const vector<Read> & normalReads, const vector<Read> & tumorReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, double a0, double b0, double free_a0, vector<HapEstResult> &normal_her, vector<HapEstResult> &tumor_her, vector<HapEstResult> &merged_her, int index) {
        LOG(logDEBUG) << "calc_hap2_bf_with_hap3" << endl;
        LOG(logINFO) << "### calculate the Bayes factor without SNP information nearby:" << endl;
        // hap1とhap2、hap3とhap4は同じとする。 (つまり、届く範囲にhetero germline snpはない場合)
        vector<Haplotype> haps;
        vector<vector<MLAlignment> > normal_liks;
        vector<vector<MLAlignment> > tumor_liks;
        Haps2::getHaplotypesBasic(haps, normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, params, refSeq, refSeqForAlign, normal_liks, tumor_liks);
        int nrt = tumorReads.size(), nrn=normalReads.size();
        vector<double> rlt(nrt*4), rln(nrn*4);
        int idx=0;
        for (size_t r=0;r<nrt;r++) {
            for (size_t h=0;h<2;h++) {
                rlt[idx] = tumor_liks[h][r].ll;idx++;
                rlt[idx] = tumor_liks[h][r].ll;idx++;
            }
        }
        idx=0;
        for (size_t r=0;r<nrn;r++) {
            for (size_t h=0;h<2;h++) {
                rln[idx] = normal_liks[h][r].ll;idx++;
                rln[idx] = normal_liks[h][r].ll;idx++;
            }
        }
        output_alignments(index, "hap2_tumor", tumorReads, tumor_liks, params);
        output_alignments(index, "hap2_normal", normalReads, normal_liks, params);
        vector<Read> f_normal_reads, f_tumor_reads;
        f_normal_reads.insert(f_normal_reads.end(), normalReads.begin(), normalReads.end());
        f_tumor_reads.insert(f_tumor_reads.end(), tumorReads.begin(), tumorReads.end());
        LOG(logINFO) << "filter uninformative reads. remaining ";
        filter_reads(f_normal_reads, rln);
        filter_reads(f_tumor_reads, rlt);
        LOGP(logINFO) << "tumor reads: " << f_tumor_reads.size() << " of " << nrt << ", normal reads: " << f_normal_reads.size() << " of " << nrn << endl;
#ifdef LOGDEBUG
        output_liktable(index, "hap2_tumor", rlt, params);
        output_liktable(index, "hap2_normal", rln, params);
#endif
        //4本用意する
        vector<Haplotype>::iterator it = haps.begin();
        haps.insert(it, haps[0]);
        haps.push_back(haps.back());
        typedef map<int, AlignedVariant>::const_iterator It;
        LOG(logINFO) << "prepared haplotypes:" << endl;
        for (size_t th=0;th<haps.size();th++) {
            const Haplotype & hap=haps[th];
            LOG(logINFO) << "hap[" << th << "] " << hap.seq << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    LOGP(logINFO) << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            LOGP(logINFO) << endl;
        }
        vector<double> normalHapFreqs, tumorHapFreqs, non_tumorHapFreqs, non_normalHapFreqs;
        vector<HapEstResult> non_tumor_her, non_normal_her;
        lower_bound_t mm_lb, nmm_lb;
        LOG(logINFO) << "calculating the marginal likelihood based on the MUTATION model...";
        MutationModel::estimate(haps, tumorReads, normalReads, rlt, rln, tumorHapFreqs, normalHapFreqs, tumor_her, normal_her, pos, leftPos, rightPos, candidateVariants, mm_lb, params.hap2_params, "mutation", params.fileName+"logs/"+itos(index)+".hap2.mutation");
        //output_mm(haps, (params.fileName+".mm.txt"), leftPos, rightPos, 0.0, mm_lb, normalHapFreqs, tumorHapFreqs);
        LOG(logINFO) << "calculating the marginal likelihood based on the ERROR model...";
        MutationModel::estimate(haps, tumorReads, normalReads, rlt, rln, non_tumorHapFreqs, non_normalHapFreqs, non_tumor_her, non_normal_her, pos, leftPos, rightPos, candidateVariants, nmm_lb, params.hap2_params, "error",params.fileName+"logs/"+itos(index)+".hap2.err");
        //output_mm(haps, (params.fileName+".non-mm.txt"), leftPos, rightPos, 0.0, nmm_lb, non_normalHapFreqs, non_tumorHapFreqs);
        double bf = mm_lb.lower_bound - nmm_lb.lower_bound;
        LOG(logDEBUG) << "********** calc_hap2_bf_with_hap3 done ***********" << endl;
        return bf;
    }



    double calc_hap2_bf(const vector<Read> & normalReads, const vector<Read> & tumorReads, const vector<Read> & mergedReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, double a0, double b0, double free_a0, vector<HapEstResult> &normal_her, vector<HapEstResult> &tumor_her, vector<HapEstResult> &merged_her) {
        LOG(logDEBUG) << "********** calc_hap2_bf ***********" << endl;
        vector<Haplotype> haps;
        vector<vector<MLAlignment> > liks;
        vector<vector<MLAlignment> > normal_liks;
        vector<vector<MLAlignment> > tumor_liks;
        Haps2::getHaplotypesBasic(haps, normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, params, refSeq, refSeqForAlign, normal_liks, tumor_liks);
        liks.insert(liks.end(), normal_liks.begin(), normal_liks.end());
        for(int i = 0;i < haps.size();i++) {
            liks[i].insert(liks[i].end(), tumor_liks[i].begin(), tumor_liks[i].end());
        }
        vector<double> mergedHapFreqs, normalHapFreqs, tumorHapFreqs;

        lower_bound_t normal_lb;
        lower_bound_t tumor_lb;
        lower_bound_t merged_lb;
        map<AlignedVariant, double> merged_vpp;
        map<AlignedVariant, double> normal_vpp;
        map<AlignedVariant, double> tumor_vpp;
        EMfor2::estimate_basic(haps, normalReads, normal_liks, normalHapFreqs, normal_her, pos, leftPos, rightPos, candidateVariants, normal_lb, normal_vpp, params, 1.0, 1.0);
        EMfor2::estimate_basic(haps, tumorReads, tumor_liks, tumorHapFreqs, tumor_her, pos, leftPos, rightPos, candidateVariants, tumor_lb, tumor_vpp, params, 1.0, 1.0);
        EMfor2::estimate_basic(haps, mergedReads, liks, mergedHapFreqs, merged_her, pos, leftPos, rightPos, candidateVariants, merged_lb, merged_vpp, params, 1.0, 1.0);
        double bf = normal_lb.lower_bound + tumor_lb.lower_bound - merged_lb.lower_bound;
        outputLowerBounds(haps, (params.fileName+".hap2_haplotypes.txt"), leftPos, rightPos, bf, normal_lb, tumor_lb, merged_lb, normalHapFreqs, tumorHapFreqs, mergedHapFreqs);
        return bf;
    }

    void outputLowerBounds(const vector<Haplotype> &haps, string fname, uint32_t leftPos, uint32_t rightPos, double bf, lower_bound_t& normal_lb, lower_bound_t &tumor_lb, lower_bound_t &merged_lb, vector<double> normalHapFreqs, vector<double> tumorHapFreqs, vector<double> mergedHapFreqs) {
#ifndef LOGDEBUG
        return;
#endif
        typedef map<int, AlignedVariant>::const_iterator It;
        std::ofstream ofs(fname.c_str(), std::ios::out | std::ios::app);
        ofs << "-------- window[" << leftPos << "-" << rightPos <<"] ---------------" << endl;
        ofs << "bayes_factor=" << bf << endl;
        ofs << "#tumor " << tumor_lb.lower_bound << endl;
        for (size_t th=0;th<tumorHapFreqs.size();th++) {
            const Haplotype & hap=haps[th];
            ofs << "hap[" << th << "] ";
            ofs << tumorHapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs << endl;
        }
        ofs << "#normal " << normal_lb.lower_bound << endl;
        for (size_t th=0;th<normalHapFreqs.size();th++) {
            const Haplotype & hap=haps[th];
            ofs << "hap[" << th << "] ";
            ofs << normalHapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs << endl;
        }
        ofs << "#merged " << merged_lb.lower_bound << endl;
        for (size_t th=0;th<mergedHapFreqs.size();th++) {
            const Haplotype & hap=haps[th];
            ofs << "hap[" << th << "] ";
            ofs << mergedHapFreqs[th] << " ";
            for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
                if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
                    ofs << "[" << it->second.getString() << " " << (it->first) << "]";
                }
            }
            ofs << endl;
        }
    }

}
