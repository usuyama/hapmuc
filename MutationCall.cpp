//
//  MutationCall.cpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/10/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//

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

namespace MutationCall
{
		string itos(int num) {
				std::string str = boost::lexical_cast<std::string>(num);
				return(str);
		}

		void output_liktable(int index, string id, const vector<double> & liks, Parameters params) {
				//hap size=4
				std::ofstream ofs((params.fileName+"logs/"+itos(index)+"."+id+".liks").c_str(), std::ios::out | std::ios::app);
				int nr = liks.size() / 4;
				for(int i=0;i<nr;i++) {
						ofs << liks[i*4];
						for(int j=1;j<4;j++)
								ofs << "\t" << liks[i*4+j];
						ofs << endl;
				}
		}


		void filter_reads(vector<Read> &reads, vector<double> &liks) {
				cout << "filter_reads: " << index << endl;
				//hap4本の場合に、どのハプロタイプにも張り付かないリードをフィルターする
				for(int i=0;i<reads.size();i++) {
						double max=-DBL_MAX, min=DBL_MAX;
						for(int j=0;j<4;j++) {
								if(liks[i*4+j]>max) max = liks[i*4+j];
								if(liks[i*4+j]<min) min = liks[i*4+j];
						}
						cout << max << " " << min << endl;
						if((max-min)<0.001) {
								//filter this read
								cout << "filter read: " << reads[i].seq_name << " " << i << " " << max << " " << min << endl;
								reads.erase(reads.begin()+i);
								liks.erase(liks.begin()+i*4, liks.begin()+i*4+4);
								i--;
						}
				}
		}


		void output_mm(vector<Haplotype> &haps, string fname, uint32_t leftPos, uint32_t rightPos, double bf, lower_bound_t& lb, vector<double> normalHapFreqs, vector<double> tumorHapFreqs) {
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
				cout << "computeBayesFactor: " << index << endl;
				vector<Haplotype> haps;
				vector<vector<MLAlignment> > liks;
				vector<HapPairLik> likPairs;

            cout << "check for long indel" << endl;
            if(candidateVariants.variants[0].getString().length() > 20) { // TODO: magic number
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
				cout <<  "centerPos=" << candidateVariants.centerPos << endl;
				candidateVariants.printAll();
				cout.flush();
				double hap2_bf = -100;
				bool hap2_flag = true, hap3_flag = true;
            string closest_germline = "-";
            int distance = -1;
				try {
						hap2_bf = calc_hap2_bf_with_hap3(normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, params, refSeq, refSeqForAlign, 1.0, 1.0, 1.0, normal_her_hap2, tumor_her_hap2, merged_her_hap2, index);
				} catch (string s) {
						cout << "error_" << s << endl;
						hap2_flag = false;
				} catch (std::exception& e) {
						string message = string("error_exception_").append(e.what());
						hap2_flag = false;
				}
				try {
						if (close_germline.size() != 0) {
								//close_somaticは空
								//close_germlineは、targetに最も近いものを採用する
								cout << "choose closest germline snp" << endl;
								vector<AlignedVariant> new_close_somatic;
								vector<AlignedVariant> new_close_germline;
								int target_pos = candidateVariants.variants[0].getStartHap();
								cout << "target_pos " << target_pos << endl;
								int min_dis= 100000, min_index=-1;
								for(int i=0;i<close_germline.size();i++) {
										int dist = abs(close_germline[i].getStartHap() - target_pos);
										if(min_dis > dist) {
												min_index = i;min_dis = dist;
										}
								}
								cout << "min " << min_dis << " " << min_index << endl;
                            stringstream ss;
                            ss << close_germline[min_index].getStartHap() << ":" << close_germline[min_index].getString();
                            closest_germline = ss.str();
                            distance = min_dis;
								new_close_germline.push_back(close_germline[min_index]);
								bool flag = Haps2::getHaplotypes(haps, normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, params, refSeq, refSeqForAlign, new_close_somatic, new_close_germline, normal_liks, tumor_liks);
								if(haps.size()==4) {
										cout << "************ hap4 **************" << endl << "haps: " << haps.size() <<  endl;
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
										vector<Read> f_normal_reads, f_tumor_reads;
										f_normal_reads.insert(f_normal_reads.end(), normalReads.begin(), normalReads.end());
										f_tumor_reads.insert(f_tumor_reads.end(), tumorReads.begin(), tumorReads.end());
										cout << "before" << f_normal_reads.size() << f_tumor_reads.size() << rln.size() << rlt.size() << endl;
										filter_reads(f_normal_reads, rln);
										filter_reads(f_tumor_reads, rlt);
										cout << "after" << f_normal_reads.size() << f_tumor_reads.size() << rln.size() << rlt.size() << endl;
										output_liktable(index, "hap3_tumor", rlt, params);
										output_liktable(index, "hap3_normal", rln, params);
										MutationModel::estimate(haps, f_tumor_reads, f_normal_reads, rlt, rln, tumorHapFreqs, normalHapFreqs, tumor_her, normal_her, pos, leftPos, rightPos, candidateVariants, mm_lb, params.hap3_params, "mutation", params.fileName+"logs/"+itos(index)+".hap3.mutation");
										//output_mm(haps, (params.fileName+".mm.txt"), leftPos, rightPos, 0.0, mm_lb, normalHapFreqs, tumorHapFreqs);
										MutationModel::estimate(haps, f_tumor_reads, f_normal_reads, rlt, rln, non_tumorHapFreqs, non_normalHapFreqs, non_tumor_her, non_normal_her, pos, leftPos, rightPos, candidateVariants, nmm_lb, params.hap3_params, "error",params.fileName+"logs/"+itos(index)+".hap3.err");
										//output_mm(haps, (params.fileName+".non-mm.txt"), leftPos, rightPos, 0.0, nmm_lb, non_normalHapFreqs, non_tumorHapFreqs);
										bf2 = mm_lb.lower_bound - nmm_lb.lower_bound;
										//cout << "**** hap4 done ***" << endl;
								} else {
                                    hap3_flag = false;
                                    cout << "error_" << "haps.size() != 4 in hap3" << endl;
                                }
						}
				} catch (string s) {
						cout << "error_" << s << endl;
						hap3_flag = false;
				} catch (std::exception& e) {
						string message = string("error_exception_").append(e.what());
						cout << "error_" << message << endl;
						hap3_flag = false;
				}
				if(tumor_her.empty()) {            
						tumor_her.swap(tumor_her_hap2);
						normal_her.swap(normal_her_hap2);
				}
            cout << "output" << endl;
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
                line.set("NN1", "-");
                line.set("NN2", "-");
                line.set("NN3", "-");
                line.set("NN4", "-");
                line.set("TN1", "-");
                line.set("TN2", "-");
                line.set("TN3", "-");
                line.set("TN4", "-");
                line.set("hap2_bf", "-");
                line.set("bf2", "-");
                line.set("closest_germline", "-");
                line.set("distance", "-");
                oData.output(line);
            } else {
				for(int i=0;i<tumor_her.size();i++) {
						cout << i;
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
						if(hap2_flag==false && hap3_flag==false) {
								line.set("NN1", "-");
								line.set("NN2", "-");
								line.set("NN3", "-");
								line.set("NN4", "-");
								line.set("TN1", "-");
								line.set("TN2", "-");
								line.set("TN3", "-");
								line.set("TN4", "-");
						} else {
								line.set("NN1", n_result.N1);
								line.set("NN2", n_result.N2);
								line.set("NN3", n_result.N3);
								line.set("NN4", n_result.N4);
								line.set("TN1", t_result.N1);
								line.set("TN2", t_result.N2);
								line.set("TN3", t_result.N3);
								line.set("TN4", t_result.N4);
						}
						if(hap2_flag) {
								line.set("hap2_bf", hap2_bf);
						} else {
								line.set("hap2_bf", "-");
						}
						if(hap3_flag){
								line.set("bf2", bf2);
						} else {
								line.set("bf2", "-");
						}
                    if(closest_germline == "-") {
                        line.set("closest_germline", "-");
                        line.set("distance", "-");
                    } else {
                        line.set("closest_germline", closest_germline);
                        line.set("distance", distance);
                    }
						oData.output(line);
				}
            }
		}

		double calc_bf(const vector<Haplotype> &haps, const vector<vector<MLAlignment> > &normal_liks, const vector<vector<MLAlignment> > tumor_liks, const vector<Read> & normalReads, const vector<Read> & tumorReads, const vector<Read> & mergedReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, double a0, double b0, double free_a0, vector<HapEstResult> &normal_her, vector<HapEstResult> &tumor_her, vector<HapEstResult> &merged_her) {
				// algorithm of undergrad-thesis
				cout << "********** calc_bf ***********" << endl;
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
				cout << "********** calc_hap2_bf_with_hap3 ***********" << endl;
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
				vector<Read> f_normal_reads, f_tumor_reads;
				f_normal_reads.insert(f_normal_reads.end(), normalReads.begin(), normalReads.end());
				f_tumor_reads.insert(f_tumor_reads.end(), tumorReads.begin(), tumorReads.end());
				filter_reads(f_normal_reads, rln);
				filter_reads(f_tumor_reads, rlt);
				output_liktable(index, "hap2_tumor", rlt, params);
				output_liktable(index, "hap2_normal", rln, params);
				//4本用意する
				vector<Haplotype>::iterator it = haps.begin();
				haps.insert(it, haps[0]);
				haps.push_back(haps.back());

				vector<double> normalHapFreqs, tumorHapFreqs, non_tumorHapFreqs, non_normalHapFreqs;
				vector<HapEstResult> non_tumor_her, non_normal_her;
				lower_bound_t mm_lb, nmm_lb;
				MutationModel::estimate(haps, tumorReads, normalReads, rlt, rln, tumorHapFreqs, normalHapFreqs, tumor_her, normal_her, pos, leftPos, rightPos, candidateVariants, mm_lb, params.hap2_params, "mutation", params.fileName+"logs/"+itos(index)+".hap2.mutation");
				//output_mm(haps, (params.fileName+".mm.txt"), leftPos, rightPos, 0.0, mm_lb, normalHapFreqs, tumorHapFreqs);
				MutationModel::estimate(haps, tumorReads, normalReads, rlt, rln, non_tumorHapFreqs, non_normalHapFreqs, non_tumor_her, non_normal_her, pos, leftPos, rightPos, candidateVariants, nmm_lb, params.hap2_params, "error",params.fileName+"logs/"+itos(index)+".hap2.err");
				//output_mm(haps, (params.fileName+".non-mm.txt"), leftPos, rightPos, 0.0, nmm_lb, non_normalHapFreqs, non_tumorHapFreqs);
				double bf = mm_lb.lower_bound - nmm_lb.lower_bound;
				return bf;
		}



		double calc_hap2_bf(const vector<Read> & normalReads, const vector<Read> & tumorReads, const vector<Read> & mergedReads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, Parameters params, string refSeq, string refSeqForAlign, double a0, double b0, double free_a0, vector<HapEstResult> &normal_her, vector<HapEstResult> &tumor_her, vector<HapEstResult> &merged_her) {
				cout << "********** calc_hap2_bf ***********" << endl;
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
