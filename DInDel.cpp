/*    
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
#include "fisher_test.cpp"
#include <boost/algorithm/string.hpp>
#include "MutationCall.hpp"
using namespace boost;
const int USECALLWINDOW=0;
using namespace seqan;
namespace po = boost::program_options;

using namespace std;
//using namespace fasta;


DetInDel::DetInDel(const string & bfName, const Parameters & _params, int multipleFiles) : params(_params)
{
	fai=NULL;
	if (params.alignAgainstReference) {
		fai = fai_load(params.refFileName.c_str());
		if (!fai) {
			cerr << "Cannot open reference sequence file." << endl;
			params.alignAgainstReference=false;
			exit(1);
		}
	}
    
	if (multipleFiles==0) {
		myBams.push_back(new MyBam(bfName));
	} else {
        
		ifstream file(bfName.c_str());
		if (!file.is_open()) {
            cout << "Cannot open file with BAM files:  " << bfName << endl;
            throw string("File open error.");
		}
		while (!file.eof()) {
			string line;
			getline(file, line);
			if (!line.empty()) {
				istringstream is(line);
				string fname;
				is >> fname;
				if (!fname.empty()) {
					cout << "Reading BAM file " << fname << endl;
					myBams.push_back(new MyBam(fname));
					myBamsFileNames.push_back(fname);
				}
			}
		}
		file.close();
	}
}



DetInDel::DetInDel(const string & normalBF, const string & tumorBF, const Parameters & _params) : params(_params)
{
	fai=NULL;
	if (params.alignAgainstReference) {
		fai = fai_load(params.refFileName.c_str());
		if (!fai) {
			cerr << "Cannot open reference sequence file." << endl;
			params.alignAgainstReference=false;
			exit(1);
		}
	}
    
    string fname = normalBF;
    if (!fname.empty()) {
        cout << "Reading BAM file " << fname << endl;
        normalBams.push_back(new MyBam(fname));
        normalBamsFileNames.push_back(fname);
        myBams.push_back(new MyBam(fname));
        myBamsFileNames.push_back(fname);
    }
    
    fname = tumorBF;
    if (!fname.empty()) {
        cout << "Reading BAM file " << fname << endl;
        tumorBams.push_back(new MyBam(fname));
        tumorBamsFileNames.push_back(fname);
        myBams.push_back(new MyBam(fname));
        myBamsFileNames.push_back(fname);
    }
}




DetInDel::~DetInDel()
{
	if (params.alignAgainstReference && fai) {
		fai_destroy(fai);
	}
	for (size_t b=0;b<myBams.size();b++) delete myBams[b];
}

void DetInDel::analyzeDifference(const pair<Haplotype, Haplotype> & hp1, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos)
{
	cout << "Inference results" << endl;
    
	if (params.analyzeLowFreqDiffThreshold<-100.0) {
		size_t offset=50;
		cout << "Haplotype pair: " << endl;
		cout << hp1.first << endl << hp1.second << endl;
        
        
		cout << "h.1 alignment: " << endl;
		cout << string(offset,' ') << hp1.first.seq << endl;
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
			MLAlignment ml=om.calcLikelihood();
            
			double lm=ml.ll;
			ObservationModelFBMax op(hp1.second, reads[r], leftPos, params.obsParams);
			MLAlignment ml2=op.calcLikelihood();
			double lp=ml2.ll;
            
			if (lm>=lp) om.printAlignment(offset);
		}
        
		cout << "h.2 alignment: " << endl;
		cout << string(offset,' ') << hp1.second.seq << endl;
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
			MLAlignment ml=om.calcLikelihood();
            
			double lm=ml.ll;
			ObservationModelFBMax op(hp1.second, reads[r], leftPos, params.obsParams);
			MLAlignment ml2=op.calcLikelihood();
			double lp=ml2.ll;
            
			if (lp>=lm) op.printAlignment(offset);
            
            
		}
        
        
        
	} else {
        
		double ll=0.0,l1=0.0, l2=0.0;
		vector<size_t > show;
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
			MLAlignment ml=om.calcLikelihood();
            
			double lm=ml.ll;
			ObservationModelFBMax op(hp1.second, reads[r], leftPos, params.obsParams);
			MLAlignment ml2=op.calcLikelihood();
			double lp=ml2.ll;
			double dll=log(exp(lp)+exp(lm))+log(.5);
			l1+=(addLogs(lm,lm)+log(.5));
			l2+=(addLogs(lp,lp)+log(.5));
			ll+=dll;
            
			if (lp-lm>params.analyzeLowFreqDiffThreshold) {
				show.push_back(r);
			}
            //		cout << "read[" << r <<"]: 1-mq: " << 1.0-reads[r].mapQual << " first hap lik: " << lm << " second hap lik: " << lp << " combined: " << dll << " lp+lm: " << ll << " lm+lm: " << l1 << " lp+lp: " << l2 << endl;
            
		}
        
		size_t offset=50;
		if (show.size()) {
			cout << "Haplotype pair: " << endl;
			cout << hp1.first << endl << hp1.second << endl;
            
            
			cout << "h.1 alignment: " << endl;
			cout << string(offset,' ') << hp1.first.seq << endl;
            
			for (size_t i=0;i<show.size();i++) {
				Read rr=reads[show[i]];
				ObservationModelFBMax om2(hp1.first, rr, leftPos, params.obsParams);
				om2.calcLikelihood();
                
                //		om2.computeMarginals();
				om2.printAlignment(offset);
			}
            
			cout << endl << endl;
			cout << "h.2 alignment: " << endl;
			cout << string(offset,' ') << hp1.second.seq << endl;
            
			for (size_t i=0;i<show.size();i++) {
				Read rr=reads[show[i]];
				ObservationModelFBMax om2(hp1.second, rr, leftPos, params.obsParams);
				om2.calcLikelihood();
                //		om2.computeMarginals();
				om2.printAlignment(offset);
			}
		}
		else { cout << "No differences in log-likelihoods over threshold." << endl; };
	}
    
}

void DetInDel::showAlignments(const pair<Haplotype, Haplotype> & hp1, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos)
{
	cout << "Inference results" << endl;
    
	double ll=0.0;
	int offset=50;
	vector <double> lf(reads.size(),0.0), ls(reads.size(),0);
	cout << "h.1 alignment: " << endl;
	cout << string(offset,' ') << hp1.first.seq << endl;
	for (size_t r=0;r<reads.size();r++) {
		ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
		double lm=om.getLogLikelihood();
		lf[r]=lm;
		if (lm<params.analyzeLowFreqDiffThreshold) {
			om.printAlignment(offset);
		}
	}
    
	cout << "h.2 alignment: " << endl;
	cout << string(offset,' ') << hp1.second.seq << endl;
	for (size_t r=0;r<reads.size();r++) {
		ObservationModelFBMax om(hp1.second, reads[r], leftPos,params.obsParams);
		double lm=om.getLogLikelihood();
		ls[r]=lm;
		if (lm<params.analyzeLowFreqDiffThreshold) {
			om.printAlignment(offset);
		}
		ll+=addLogs(lf[r],ls[r])+log(.5);
	}
	cout << "Total loglikelihood: " << ll << endl;
    
    
}

void DetInDel::showAlignmentsPerHaplotype(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, uint32_t candPos, uint32_t leftPos)
{
	cout << "ALIGNMENTS" << endl;
    
	vector<std::set<size_t> > maxHap(haps.size());
	for (size_t r=0;r<reads.size();r++) {
		size_t idx=0;
		double ml=-HUGE_VAL;
		for (size_t h=0;h<haps.size();h++) {
			if (liks[h][r].ll>ml) {
				ml=liks[h][r].ll;
				idx=h;
			}
		}
		maxHap[idx].insert(r);
	}
    
	int offset=50;
	for (size_t h=0;h<haps.size();h++) {
		cout << "*******************************************" << endl;
		cout << endl << "HAPLOTYPE " << h << endl << endl;
		cout << string(offset,' ') << haps[h].seq << endl;
		BOOST_FOREACH(size_t r, maxHap[h]) {
			ObservationModelFBMax om(haps[h], reads[r], leftPos,params.obsParams);
			om.calcLikelihood();
			om.printAlignment(offset);
		}
	}
    
    
}




string DetInDel::getRefSeq(uint32_t lpos, uint32_t rpos)
{
	if (!fai) throw string("FAI error.");
    
	char *str;
	char *ref;
    
	str = (char*)calloc(strlen(params.tid.c_str()) + 30, 1);
	sprintf(str, "%s:%d-%d", params.tid.c_str(), lpos, rpos);
	int len;
	ref = fai_fetch(fai, str, &len);
	if (len==0) throw string("faidx error: len==0");
	free(str);
	string res(ref);
	free(ref);
    
	transform(res.begin(), res.end(), res.begin(), ::toupper);
	return res;
}


double DetInDel::getMaxHap(Haplotype & h1, Haplotype &h2, HapPairLik & hpl, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs)
{
    
	size_t idx=0, midx;
	double maxll=-HUGE_VAL;
	for (idx=0;idx<likPairs.size();idx++) {
        double ll=likPairs[idx].ll;
        if (ll>maxll) {
            maxll=ll;
            midx=idx;
        }
	}
	h1=haps[likPairs[midx].h1];
	h2=haps[likPairs[midx].h2];
    
	/*
     cout << "getMaxHap: " << midx <<  " h1: " << likPairs[midx].h1 << " h2: " << likPairs[midx].h2 << endl;
     cout << "indelcoverage h1: ";
     for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage1.begin();it!=likPairs[midx].hapIndelCoverage1.end();it++) {
     cout << "[" << it->second.nf << "," << it->second.nr << "]";
     }
     cout << endl;
     cout << "indelcoverage h2: ";
     for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage2.begin();it!=likPairs[midx].hapIndelCoverage2.end();it++) {
     cout << "[" << it->second.nf << "," << it->second.nr << "]";
     }
     cout << endl;
     */
    
	hpl=likPairs[midx];
	return maxll;
}

void DetInDel::outputMaxHap(ostream *output, const string & prefix, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs)
{
    
	Haplotype h1, h2;
	HapPairLik hpl;
	getMaxHap(h1,h2, hpl, haps, likPairs);
	*output << prefix << " " << hpl.ll << " " << hpl.numFirst << " " << hpl.numSecond << " " << hpl.numIndFirst << " " << hpl.numIndSecond << " " << hpl.numOffBoth << " " << h1.seq << " " << h2.seq << " ";
	for (map<int, AlignedVariant>::const_iterator it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << "!";
	for (map<int, AlignedVariant>::const_iterator it=h2.indels.begin();it!=h2.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << "!";
	for (map<int, AlignedVariant>::const_iterator it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << "!";
	for (map<int, AlignedVariant>::const_iterator it=h2.snps.begin();it!=h2.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << endl;
    
    
}

void DetInDel::outputTopHaps(ostream *output, const string & prefix, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs, int n)
{
	// output n most likely haplotype pairs
    
	for (int ns=0;ns<n && ns<int(likPairs.size());ns++) {
		const Haplotype & h1 = haps[likPairs[ns].h1];
		const Haplotype & h2 = haps[likPairs[ns].h2];
		const HapPairLik & hpl = likPairs[ns];
		*output << prefix << " " << ns+1 << " " << hpl.ll << " " << hpl.numFirst << " " << hpl.numSecond << " " << hpl.numIndFirst << " " << hpl.numIndSecond << " " << hpl.numOffBoth << " " << h1.seq << " " << h2.seq << " ";
		for (map<int, AlignedVariant>::const_iterator it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << "!";
		for (map<int, AlignedVariant>::const_iterator it=h2.indels.begin();it!=h2.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << "!";
		for (map<int, AlignedVariant>::const_iterator it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << "!";
		for (map<int, AlignedVariant>::const_iterator it=h2.snps.begin();it!=h2.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << endl;
	}
    
}

void DetInDel::outputHapsAndFreqs(ostream *output, const string & prefix, const vector<Haplotype> & haps, const vector<double> & freqs, uint32_t leftPos)
{
	// output n most likely haplotype pairs
    
	for (size_t h=0;h<haps.size();h++) {
		const Haplotype & h1 = haps[h];
		*output << prefix << " " << h+1 << " " << freqs[h] << " ";
		for (map<int, AlignedVariant>::const_iterator it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString()!="*REF") *output << leftPos+it->first << "," << it->second.getString() << "|";
		//*output << "!";
		//for (map<int, AlignedVariant>::const_iterator it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString()!="*REF") *output << leftPos+it->first << "," << it->second.getString() << "|";
		*output << endl;
	}
    
}



#include <fstream>
#include <sstream>


DetInDel::CIGAR DetInDel::getCIGAR(const Haplotype & hap, const Read & read, const MLAlignment & ml, int refSeqStart)
{
	if (hap.ml.hpos.size()!=hap.size()) throw string("Haplotype has not been aligned!");
	if (ml.hpos.size()!=read.size()) throw string("Read is not properly aligned!");
	const MLAlignment & hml=hap.ml; // alignment of haplotype to reference
    
	//string qname = bam1_qname(read.getBam());
	const int debug = 0;
	/*
     if (qname == "IL8_4337:8:102:11530:1494") {
     cout << "YES" << endl;
     cout << qname << endl;
     debug = 1;
     }
     */
    
	vector<int> npos(read.size()); // npos records position of read base on the reference sequence
	for (int b=0;b<int(read.size());b++) {
		if (ml.hpos[b]>=0) npos[b]=hml.hpos[ml.hpos[b]]; else npos[b]=ml.hpos[b];
	}
    
	if (debug) {
		for (size_t h=0;h<npos.size();h++) {
			cout << "[" << h << "," << npos[h] << "]";
		}
		cout << endl;
		cout << endl;
        
		for (size_t h=0;h<ml.hpos.size();h++) {
			cout << "[" << h << "," << ml.hpos[h] << "]";
		}
		cout << endl;
		cout << endl;
		for (size_t h=0;h<hml.hpos.size();h++) {
            cout << "[" << h << "," << hml.hpos[h] << "]";
        }
        cout << endl;
	}
	CIGAR cig;
    
	int b=0, prevponr=0; // position on reference previous base aligned on the reference (ie no deletion/insertion/LO/RO)
    
	// determine last base in read aligned to the haplotype
	b=read.size()-1;
	while (npos[b]<0) b--;
	int lastbonh=b;
    
	if (lastbonh<0) {
		// clip the whole read, read is de facto off haplotype
		cig.push_back(CIGAR::CIGOp(BAM_CSOFT_CLIP, read.size()));
		return cig;
	}
    
    
	if (debug) {
		cout << "lastbonh: " << lastbonh << endl;
	}
	// find first base in read that is aligned to haplotype and to reference
	// all sequence before that is considered 'soft clipped', ie not aligned. This may include sequence that matches perfectly to the reference
	b=0;
	while (npos[b]<0) b++;
 	if (b>0) cig.push_back(CIGAR::CIGOp(BAM_CSOFT_CLIP, b));
 	prevponr=npos[b];
 	cig.refPos=refSeqStart+prevponr;
    
 	int curr_cop=BAM_CMATCH;
 	int len_curr_cop=1;
    
    
    
 	while (b<lastbonh) {
        
 		int chp=npos[b]; // position on reference of current base in read
 		int nhp=npos[b+1];
        
 		if (nhp==MLAlignment::INS) {
 			if (chp==MLAlignment::INS) {
 				// stay on inserted sequence
 				if (curr_cop!=BAM_CINS) throw string("Error(1)!");
 				len_curr_cop++;
 			} else {
 				if (chp>=0) {
 					// going from on reference to insertions
 					if (curr_cop!=BAM_CMATCH) throw string("Error(2)!");
 					// write CIGAR
					cig.push_back(CIGAR::CIGOp(BAM_CMATCH,len_curr_cop));
                    
					// update current CIGAR operation
					len_curr_cop=1;
					curr_cop=BAM_CINS;
                    
					prevponr=chp;
 				} else throw string("How is this possible? (1)");
 			}
            
 		} else if (chp>=0 && nhp>=0 && nhp-chp==1) {
            // MATCH to MATCH
 			if (curr_cop!=BAM_CMATCH) {
 				cout << "b: " << b << " chp: " << chp << " nhp: " << nhp << endl;
 				throw string("Error(3)!");
 			}
 			len_curr_cop++;
 			prevponr=nhp;
 		} else if (chp>=0 && nhp>=0 && nhp-chp>1) {
 			// deletion
 			if (curr_cop!=BAM_CMATCH) throw string("Error(4)!");
 			// write CIGAR
 			cig.push_back(CIGAR::CIGOp(BAM_CMATCH,len_curr_cop));
            
 			// write deletion CIGAR
 			cig.push_back(CIGAR::CIGOp(BAM_CDEL,nhp-chp-1));
            
 			curr_cop=BAM_CMATCH;
 			len_curr_cop=1;
            
 			prevponr=nhp;
 		} else if (chp==MLAlignment::INS && nhp-prevponr==1) {
 			// from inserted bases to bases matched to the reference
 			cig.push_back(CIGAR::CIGOp(BAM_CINS,len_curr_cop));
            
 			//
 			curr_cop=BAM_CMATCH;
 			len_curr_cop=1;
            
 			prevponr=nhp;
 		} else if (chp==MLAlignment::INS && nhp-prevponr>1) {
 			// next base is again on reference but more than 1 reference base from the last read base aligned to the haplotype
 			cig.push_back(CIGAR::CIGOp(BAM_CINS,len_curr_cop));
 			cig.push_back(CIGAR::CIGOp(BAM_CDEL,nhp-prevponr-1));
            
 			curr_cop=BAM_CMATCH;
 			len_curr_cop=1;
            
 			prevponr=nhp;
 		}
 		b++;
 	}
    
 	// write last cigar
 	cig.push_back(CIGAR::CIGOp(curr_cop,len_curr_cop));
    
 	// write soft_clip at the end
 	if (read.size()-1 - lastbonh>0) {
 		cig.push_back(CIGAR::CIGOp(BAM_CSOFT_CLIP,read.size()-1 - lastbonh));
 	}
    
 	/*
     cout << "cig: ";
     BOOST_FOREACH(CIGAR::CIGOp cop, cig) {
     cout << "(" << cop.first << "," << cop.second << ")" ;
     }
     cout << endl;
     */
 	return cig;
}


void DetInDel::getReadsFromBams(vector<MyBam *> & Bams, uint32_t leftPos, uint32_t rightPos, vector<Read> & reads, uint32_t & oldLeftPos, uint32_t & oldRightFetchReadPos, vector<Read *> & readBuffer, const bool reset)
{
    cout << "   getReadsFromBams()"<< endl;
    BOOST_FOREACH(MyBam * b, Bams) {
        cout << " " << b->fileName;
    }
    cout << endl;
	// filter using map quality
	class SortFunc {
	public:
		static bool sortFunc(const Read & r1, const Read & r2)
		{
			// sort in decreasing order
			if (r1.mapQual>r2.mapQual) return true; else return false;
		}
	};
    
	if (leftPos<oldLeftPos) {
		cerr << "Windows are not sorted!" << endl;
		exit(3);
	}
    
	reads.clear();
    
	//if (int(rightPos-leftPos)<3*params.minReadOverlap) throw string("Choose a larger width or a smaller minReadOverlap.");
    
    
	int maxDev = int ( libraries.getMaxInsertSize());
    
	//maxDev = 100;
	//cerr << "CHANGE THIS CHANGE THIS" << endl;
    
	string_hash< list<int> > mapped_name_to_idx, unmapped_name_to_idx; // query name to read idx
	string_hash< list<int> >::const_iterator hash_it;
    
	int numUnknownLib = 0;
	string_hash <int> unknownLib;
	const int LEFTPAD = 200;
	// note the idea is to get only consider reads starting at position leftMostReadPos (and not ones merely overlapping)
	// LEFTPAD should take care of overlap effects (note that leftMostReadPos is already generous, based on library insert size)
    
    
	uint32_t rightFetchReadPos = rightPos+maxDev;
	uint32_t rightMostReadPos = rightPos+maxDev;
    
	uint32_t leftFetchReadPos = leftPos-maxDev-LEFTPAD;
	uint32_t leftMostReadPos = leftPos-maxDev-LEFTPAD; // left most position of reads we want to seriously consider
    
	// reset indicates whether we want to remake the readBuffer
    
	vector<Read*> newReadBuffer;
	bool leftOverlapsPrevious = false;
	if (reset) {
		for (size_t r=0;r<readBuffer.size();r++) {
			if (readBuffer[r]!=NULL) delete readBuffer[r];
		}
		readBuffer.clear();
		oldRightFetchReadPos = rightFetchReadPos;
	} else {
		// clear reads that do not overlap with new window [leftMostReadPos, rightMostReadPos]
        
		for (size_t r=0;r<readBuffer.size();r++) {
			uint32_t rend = readBuffer[r]->getEndPos();
			uint32_t rbeg = readBuffer[r]->getBam()->core.pos;
			if (rbeg<leftMostReadPos) {
				delete readBuffer[r];
				readBuffer[r]=NULL;
			} else {
				newReadBuffer.push_back(readBuffer[r]);
			}
		}
        
		// note that it is required that the new leftPos of the window >= the old leftPos
		// therefore if leftFetchReadPos<=oldRightFetchReadPos the new w
		if (leftMostReadPos<oldRightFetchReadPos) {
			leftFetchReadPos = oldRightFetchReadPos;
			leftOverlapsPrevious = true;
		}
	}
    
    
    
    
	// cout << "leftFetchReadPos: " << leftFetchReadPos << " rightFetchReadPos: " << rightFetchReadPos << " oldRightFetchReadPos: " << oldRightFetchReadPos << endl;
	// cout << "leftMostReadPos: " << leftMostReadPos << " rightMostReadPos: " << rightMostReadPos << " leftOverlapsPrevious: " << int(leftOverlapsPrevious) << endl;
	// store updated readBuffer
	readBuffer.swap(newReadBuffer);
    
    //	cout << "leftPos : " << leftPos << " rightPos: " << rightPos << " maxDev: " << maxDev << endl;
    
    
	// first clean readbuffer
    
    
	int numReads = readBuffer.size();
    
	vector<Read> newReads;
	if (leftFetchReadPos<=rightFetchReadPos) {
		cout << "Fetching reads...." << endl;
		for (size_t b=0;b<Bams.size();b++) {
			Read::FetchReadData data(&newReads, int(b), &(this->libraries), &Bams, numReads, params.maxReads*100);
			bam_fetch(Bams[b]->bf, Bams[b]->idx, Bams[b]->getTID(params.tid), leftFetchReadPos , rightFetchReadPos,&data , &Read::fetchFuncVectorPooled);
			numUnknownLib += data.numUnknownLib;
			numReads = data.numReads;
		}
		oldRightFetchReadPos = rightFetchReadPos;
	}
    
	// add new reads to readBuffer
    
	for (size_t r=0;r<newReads.size();r++) {
		if (newReads[r].getBam()->core.pos>=leftFetchReadPos) {
			// only store reads that do not overlap with the boundary;
			// reads overlapping with boundary will have been picked up before.
			readBuffer.push_back(new Read(newReads[r]));
		}
	}
    
	if (0) {
		// check with regular fetch
		vector<Read> tmpReads;
		for (size_t b=0;b<Bams.size();b++) {
			Read::FetchReadData data(&tmpReads, int(b), &(this->libraries), &Bams, numReads, params.maxReads*100);
			bam_fetch(Bams[b]->bf, Bams[b]->idx, Bams[b]->getTID(params.tid), leftMostReadPos , rightMostReadPos,&data , &Read::fetchFuncVectorPooled);
		}
		for (size_t r=0;r<tmpReads.size();r++) {
			if (tmpReads[r].getBam()->core.pos>=leftMostReadPos) {
				string qname = string(bam1_qname(tmpReads[r].getBam()));
				cout << "glp: "  << leftPos << " qname: " << qname << " pos: " << tmpReads[r].pos << " end: " << tmpReads[r].getEndPos() << endl;
			}
		}
	}
    
    
    
	// check readbuffer for duplicates (debugging)
	if (0) {
		string_hash <int> qnameCount;
		for (size_t r=0;r<readBuffer.size();r++) {
			string qname = string(bam1_qname(readBuffer[r]->getBam()));
			//cout << "lp: "  << leftPos << " qname: " << qname << " pos: " << readBuffer[r]->pos << " end: " << readBuffer[r]->getEndPos() << endl;
			string_hash<int>::iterator it = qnameCount.find(qname);
			if (it == qnameCount.end()) {
				qnameCount[qname]=1;
			} else {
				qnameCount[qname]++;
				if (qnameCount[qname]>2) {
					cerr << "Duplicate reads: readbuffer problem!" << endl;
					throw string("duplicate reads!");
				}
			}
		}
	}
    
    
    
	newReads.clear();
    
	size_t oldNumReads=readBuffer.size();
    
    
	// copy readBuffer to reads
    
	for (size_t r=0;r<readBuffer.size();r++) {
		reads.push_back(Read(*readBuffer[r]));
	}
    
    
	// get query names
    
	vector<int> unmapped;
	for (size_t r=0; r<reads.size();r++) {
		if (reads[r].isUnmapped())  {
			unmapped.push_back(r);
			unmapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
            //	 cout << " __reads[" << r  << "]: " << bam1_qname(reads[r].getBam()) << " UNMAPPED" << endl;
		} else {
			mapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
            //	cout << " __reads[" << r  << "]: " << bam1_qname(reads[r].getBam()) << " pos: " << reads[r].pos << " " << reads[r].getBam()->core.pos << " mpos: " << reads[r].getBam()->core.mpos << " mu: " << reads[r].mateIsUnmapped() << endl;
		}
	}
    
    
	// filter reads based on haplotype window, the minimum read overlap for mapped reads, minimum mapping quality and maximum read length
    
	/*
     cout << "name_to_idx.size: " << mapped_name_to_idx.size() << endl;
     for (hash_it = mapped_name_to_idx.begin();hash_it!=mapped_name_to_idx.end();hash_it++) {
     cout << "hit: " << hash_it->first;
     BOOST_FOREACH(int x, hash_it->second) {
     cout << " " << x;
     }
     cout << endl;
     }
     */
	int numTIDmismatch = 0, numOrphan =0, numOrphanUnmapped = 0, numInRegion = 0;
    
	// reads are filtered by setting mapping quality to -1
	vector<Read> filteredReads;
	double minMapQual = params.mapQualThreshold;
	if (minMapQual<0.0) minMapQual=0.0;
	for (int r=0;r<int(reads.size());r++) {
		//cout << "***" << endl;
		//cout << "reads.mapQual " << reads[r].mapQual <<  " bq: " << reads[r].getBam()->core.qual << endl;
		bool filter = false;
		int tf = 0;
		if (reads[r].size()>params.maxReadLength) filter=true;
        
		if (reads[r].getEndPos()<leftMostReadPos || reads[r].pos>rightMostReadPos) filter=true;
        
		if (!reads[r].isUnmapped()) {
            //		cout << "mapped" << endl;
			if (reads[r].pos+int(reads[r].size())<int(leftPos)+params.minReadOverlap || reads[r].pos>int(rightPos)-params.minReadOverlap) {
                //		cout << " { " << reads[r].pos+reads[r].size() << " " << leftPos+params.minReadOverlap << " " << reads[r].pos << " " << rightPos-params.minReadOverlap << " } " << endl;
				filter=true;
				tf = 1;
			} else if (reads[r].mateIsUnmapped() == false ){
				if (reads[r].getBam()->core.mtid != reads[r].getBam()->core.tid) {
					// filter = true;
					// cout << "TIDERR: reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << endl;
					numTIDmismatch++;
				} else {
					// find mate of mapped read
					// filter if we cannot find it (mapped to another chromosome, those are a bit suspicious)
                    
					// lookup mapped read
					hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
					if (hash_it == mapped_name_to_idx.end()) { numOrphan++; filter=true; } else {
						if (hash_it->second.size()>2) cerr << "HUH? DUPLICATE READ LABELS???" << endl;
						if (reads[r].mateIsUnmapped() == false) {
							filter = true;
						}
                        
						BOOST_FOREACH(int idx, hash_it->second) {
							if (idx != r) {
								reads[r].mateLen = reads[idx].size();
								reads[r].matePos = reads[idx].pos;
								filter = false;
								if (reads[r].matePos != reads[r].getBAMMatePos()) {
									cerr << "matepos inconsistency!" << endl;
									cerr << reads[r].matePos << " " << reads[r].getBAMMatePos() << endl;
									exit(1);
								}
							}
						}
                        
						if (filter == true) {
							numOrphan++;
							tf = 2;
						}
						//cout << "mapped read: " << r << " " << qname << " pos: " << reads[r].pos << " " << reads[r].getBam()->core.mtid << " " <<  reads[r].getBam()->core.mpos << " mateunmap: " << reads[r].mateIsUnmapped() << endl;
					}
				}
			} else if (reads[r].mateIsUnmapped() == true) {
				reads[r].matePos=reads[r].pos;
				hash_it = unmapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
				if (hash_it == unmapped_name_to_idx.end()) { filter=true; } else {
					filter = true;
					if (hash_it->second.size()>2) cerr << "HUH? DUPLICATE READ LABELS???" << endl;
					BOOST_FOREACH(int idx, hash_it->second) {
						if (idx != r) {
							reads[r].mateLen = reads[idx].size();
							filter = false;
						}
					}
                    
				}
				if (filter==true) {
					numOrphan++;
					tf = 3;
				}
			}
			if (filter == false) numInRegion++;
            
		} else {
            //		cout << " unmapped" << endl;
			// read is unmapped
			if (params.mapUnmappedReads) {
                //			cout << " unmapped " << qname << " ; ";
                
				// lookup mapped read
				hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
				int idx;
				if (hash_it == mapped_name_to_idx.end()) { numOrphanUnmapped++; filter=true; tf = 4;} else {
					if (hash_it->second.size()!=1) {
						cerr << "UNMAPPED READ HAS MORE THAN ONE MATE!" << endl;
						exit(1);
					}
					idx = *(hash_it->second.begin());
					//cout << " FOUND " << idx << endl ;
					// check if mate will overlap with haplotype
					uint32_t range_l, range_r; // range of mate
                    
					int maxInsert = (int) reads[idx].getLibrary().getMaxInsertSize();
					int minInsert = 0;
					uint32_t rpos = reads[idx].pos;
                    
                    //				cout << "idx: " << idx << " unmapped: " << reads[idx].isUnmapped() << " rpos: " << rpos << " isreverse: " << reads[idx].isReverse() << endl;
					if (reads[idx].isReverse()) {
						range_l = rpos-maxInsert;
						range_r = rpos-minInsert;
					} else {
						range_l = rpos+minInsert;
						range_r = rpos+maxInsert;
					}
                    
                    //				cout << "insert: " << insert << " std: " << std << " range_l : " << range_l << " range_r: " << range_r << " leftPos: " << leftPos << "  rightPos: " << rightPos << endl;
                    
					if (range_r>leftPos && range_l<rightPos) {
						numInRegion++;
						filter=false;
						reads[r].mapQual = reads[idx].mapQual;
						reads[r].matePos = reads[idx].pos;
						reads[r].mateLen = reads[idx].size();
						if (reads[r].isReverse() == reads[idx].isReverse()) {
							reads[r].reverse();
							reads[r].complement();
						}
					} else {
						filter=true;
						tf = 5;
					}
				}
			} else {
				filter = true;
			}
            
		}
		if (filter == true) reads[r].mapQual = -1.0;
        
        //		cout << "reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << " Filter: " << tf << " filter: " << filter <<  " mq: " << reads[r].mapQual << endl;
	}
    
    
	int nUnmapped = 0;
	int nMateposError = 0;
	sort(reads.begin(), reads.end(), SortFunc::sortFunc);
	size_t max; for (max=0;max<params.maxReads && max<reads.size();max++) if (!(reads[max].mapQual<minMapQual)) {
		if (reads[max].matePos==-1 && reads[max].isPaired() && !reads[max].mateIsUnmapped() ) {
			nMateposError++;
			reads[max].matePos = reads[max].pos;
		};
		filteredReads.push_back(Read(reads[max]));
		if (reads[max].isUnmapped()) nUnmapped++;
	} else break;
    
    
	filteredReads.swap(reads);
	filteredReads.clear();
	//Read::filterReads(reads, params.maxReads, params.mapQualThreshold, params.maxReadLength, params.minReadOverlap, leftPos, rightPos);
    
	if (params.filterReadAux!="") {
		if (params.filterReadAux.size()>1) {
			size_t before=reads.size();
			int exclude=1;
			if (params.filterReadAux[0]=='+') exclude=0;
			string match=params.filterReadAux.substr(1,params.filterReadAux.size());
			Read::filterReads(reads, exclude, match);
			size_t after=reads.size();
			if (!params.quiet) cout << "filterAux: " << before-after << " reads were filtered based on match string " << params.filterReadAux << endl;
		}
	}
    
	if (!params.quiet) cout << "Number of reads: " << reads.size() << " out of " << oldNumReads << " # unmapped reads: " << nUnmapped << " numReadsUnknownLib: " << numUnknownLib << " numChrMismatch: " << numTIDmismatch << " numMappedWithoutMate: " << numOrphan << " numUnmappedWithoutMate: " << numOrphanUnmapped << endl;
	if (nMateposError) {
		cerr << "The mate position of " << nMateposError << " reads was recorded as -1 in the BAM file" << endl;
	}
    
	if (params.showReads) {
		for (size_t r=0;r<reads.size();r++) {
			cout << "read[" << r << "]: " << reads[r] << endl;
		}
	}
    
	if (reads.size()<2) {
        cout << "size=" << reads.size() << endl;
		throw string("too_few_reads");
	} else if (reads.size()>=params.maxReads) {
        cout << "size=" << reads.size() << endl;
		throw string("above_read_count_threshold");
	}
    
}



void DetInDel::getReads(uint32_t leftPos, uint32_t rightPos, vector<Read> & reads, uint32_t & oldLeftPos, uint32_t & oldRightFetchReadPos, vector<Read *> & readBuffer, bool reset)
{
    cout << "   getReads()" << endl;
	// filter using map quality
	class SortFunc {
	public:
		static bool sortFunc(const Read & r1, const Read & r2)
		{
			// sort in decreasing order
			if (r1.mapQual>r2.mapQual) return true; else return false;
		}
	};
    
	if (leftPos<oldLeftPos) {
		cerr << "Windows are not sorted!" << endl;
		exit(3);
	}
    
	reads.clear();
    
	//if (int(rightPos-leftPos)<3*params.minReadOverlap) throw string("Choose a larger width or a smaller minReadOverlap.");
    
    
	int maxDev = int ( libraries.getMaxInsertSize());
    
	//maxDev = 100;
	//cerr << "CHANGE THIS CHANGE THIS" << endl;
    
	string_hash< list<int> > mapped_name_to_idx, unmapped_name_to_idx; // query name to read idx
	string_hash< list<int> >::const_iterator hash_it;
    
	int numUnknownLib = 0;
	string_hash <int> unknownLib;
	const int LEFTPAD = 200;
	// note the idea is to get only consider reads starting at position leftMostReadPos (and not ones merely overlapping)
	// LEFTPAD should take care of overlap effects (note that leftMostReadPos is already generous, based on library insert size)
    
    
	uint32_t rightFetchReadPos = rightPos+maxDev;
	uint32_t rightMostReadPos = rightPos+maxDev;
    
	uint32_t leftFetchReadPos = leftPos-maxDev-LEFTPAD;
	uint32_t leftMostReadPos = leftPos-maxDev-LEFTPAD; // left most position of reads we want to seriously consider
    
	// reset indicates whether we want to remake the readBuffer
    
	vector<Read*> newReadBuffer;
	bool leftOverlapsPrevious = false;
	if (reset) {
		for (size_t r=0;r<readBuffer.size();r++) {
			if (readBuffer[r]!=NULL) delete readBuffer[r];
		}
		readBuffer.clear();
		oldRightFetchReadPos = rightFetchReadPos;
	} else {
		// clear reads that do not overlap with new window [leftMostReadPos, rightMostReadPos]
        
		for (size_t r=0;r<readBuffer.size();r++) {
			uint32_t rend = readBuffer[r]->getEndPos();
			uint32_t rbeg = readBuffer[r]->getBam()->core.pos;
			if (rbeg<leftMostReadPos) {
				delete readBuffer[r];
				readBuffer[r]=NULL;
			} else {
				newReadBuffer.push_back(readBuffer[r]);
			}
		}
        
		// note that it is required that the new leftPos of the window >= the old leftPos
		// therefore if leftFetchReadPos<=oldRightFetchReadPos the new w
		if (leftMostReadPos<oldRightFetchReadPos) {
			leftFetchReadPos = oldRightFetchReadPos;
			leftOverlapsPrevious = true;
		}
	}
    
    
    
    
	// cout << "leftFetchReadPos: " << leftFetchReadPos << " rightFetchReadPos: " << rightFetchReadPos << " oldRightFetchReadPos: " << oldRightFetchReadPos << endl;
	// cout << "leftMostReadPos: " << leftMostReadPos << " rightMostReadPos: " << rightMostReadPos << " leftOverlapsPrevious: " << int(leftOverlapsPrevious) << endl;
	// store updated readBuffer
	readBuffer.swap(newReadBuffer);
    
    //	cout << "leftPos : " << leftPos << " rightPos: " << rightPos << " maxDev: " << maxDev << endl;
    
    
	// first clean readbuffer
    
    
	int numReads = readBuffer.size();
    
	vector<Read> newReads;
	if (leftFetchReadPos<=rightFetchReadPos) {
		cout << "Fetching reads...." << endl;
		for (size_t b=0;b<myBams.size();b++) {
			//bam_fetch(myBams[b].bf, myBams[b].idx, myBams[b].getTID(params.tid), leftPos+params.minReadOverlap, rightPos-params.minReadOverlap, &reads, &Read::fetchFuncVector);
			Read::FetchReadData data(&newReads, int(b), &(this->libraries), &myBams, numReads, params.maxReads*100);
			bam_fetch(myBams[b]->bf, myBams[b]->idx, myBams[b]->getTID(params.tid), leftFetchReadPos , rightFetchReadPos,&data , &Read::fetchFuncVectorPooled);
			numUnknownLib += data.numUnknownLib;
			numReads = data.numReads;
		}
		oldRightFetchReadPos = rightFetchReadPos;
	}
    
	// add new reads to readBuffer
    
	for (size_t r=0;r<newReads.size();r++) {
		if (newReads[r].getBam()->core.pos>=leftFetchReadPos) {
			// only store reads that do not overlap with the boundary;
			// reads overlapping with boundary will have been picked up before.
			readBuffer.push_back(new Read(newReads[r]));
		}
	}
    
	if (0) {
		// check with regular fetch
		vector<Read> tmpReads;
		for (size_t b=0;b<myBams.size();b++) {
			//bam_fetch(myBams[b].bf, myBams[b].idx, myBams[b].getTID(params.tid), leftPos+params.minReadOverlap, rightPos-params.minReadOverlap, &reads, &Read::fetchFuncVector);
			Read::FetchReadData data(&tmpReads, int(b), &(this->libraries), &myBams, numReads, params.maxReads*100);
			bam_fetch(myBams[b]->bf, myBams[b]->idx, myBams[b]->getTID(params.tid), leftMostReadPos , rightMostReadPos,&data , &Read::fetchFuncVectorPooled);
		}
		for (size_t r=0;r<tmpReads.size();r++) {
			if (tmpReads[r].getBam()->core.pos>=leftMostReadPos) {
				string qname = string(bam1_qname(tmpReads[r].getBam()));
				cout << "glp: "  << leftPos << " qname: " << qname << " pos: " << tmpReads[r].pos << " end: " << tmpReads[r].getEndPos() << endl;
			}
		}
	}
    
    
    
	// check readbuffer for duplicates (debugging)
	if (1) {
		string_hash <int> qnameCount;
		for (size_t r=0;r<readBuffer.size();r++) {
			string qname = string(bam1_qname(readBuffer[r]->getBam()));
			//cout << "lp: "  << leftPos << " qname: " << qname << " pos: " << readBuffer[r]->pos << " end: " << readBuffer[r]->getEndPos() << endl;
			string_hash<int>::iterator it = qnameCount.find(qname);
			if (it == qnameCount.end()) {
				qnameCount[qname]=1;
			} else {
				qnameCount[qname]++;
				if (qnameCount[qname]>2) {
					cerr << "Duplicate reads: readbuffer problem!" << endl;
					throw string("duplicate reads!");
				}
			}
		}
	}
    
    
    
	newReads.clear();
    
	size_t oldNumReads=readBuffer.size();
    
    
	// copy readBuffer to reads
    
	for (size_t r=0;r<readBuffer.size();r++) {
		reads.push_back(Read(*readBuffer[r]));
	}
    
    
	// get query names
    
	vector<int> unmapped;
	for (size_t r=0; r<reads.size();r++) {
		if (reads[r].isUnmapped())  {
			unmapped.push_back(r);
			unmapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
            //	 cout << " __reads[" << r  << "]: " << bam1_qname(reads[r].getBam()) << " UNMAPPED" << endl;
		} else {
			mapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
            //	cout << " __reads[" << r  << "]: " << bam1_qname(reads[r].getBam()) << " pos: " << reads[r].pos << " " << reads[r].getBam()->core.pos << " mpos: " << reads[r].getBam()->core.mpos << " mu: " << reads[r].mateIsUnmapped() << endl;
		}
	}
    
    
	// filter reads based on haplotype window, the minimum read overlap for mapped reads, minimum mapping quality and maximum read length
    
	/*
     cout << "name_to_idx.size: " << mapped_name_to_idx.size() << endl;
     for (hash_it = mapped_name_to_idx.begin();hash_it!=mapped_name_to_idx.end();hash_it++) {
     cout << "hit: " << hash_it->first;
     BOOST_FOREACH(int x, hash_it->second) {
     cout << " " << x;
     }
     cout << endl;
     }
     */
	int numTIDmismatch = 0, numOrphan =0, numOrphanUnmapped = 0, numInRegion = 0;
    
	// reads are filtered by setting mapping quality to -1
	vector<Read> filteredReads;
	double minMapQual = params.mapQualThreshold;
	if (minMapQual<0.0) minMapQual=0.0;
	for (int r=0;r<int(reads.size());r++) {
		//cout << "***" << endl;
		//cout << "reads.mapQual " << reads[r].mapQual <<  " bq: " << reads[r].getBam()->core.qual << endl;
		bool filter = false;
		int tf = 0;
		if (reads[r].size()>params.maxReadLength) filter=true;
        
		if (reads[r].getEndPos()<leftMostReadPos || reads[r].pos>rightMostReadPos) filter=true;
        
		if (!reads[r].isUnmapped()) {
            //		cout << "mapped" << endl;
			if (reads[r].pos+int(reads[r].size())<int(leftPos)+params.minReadOverlap || reads[r].pos>int(rightPos)-params.minReadOverlap) {
                //		cout << " { " << reads[r].pos+reads[r].size() << " " << leftPos+params.minReadOverlap << " " << reads[r].pos << " " << rightPos-params.minReadOverlap << " } " << endl;
				filter=true;
				tf = 1;
			} else if (reads[r].mateIsUnmapped() == false ){
				if (reads[r].getBam()->core.mtid != reads[r].getBam()->core.tid) {
					// filter = true;
					// cout << "TIDERR: reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << endl;
					numTIDmismatch++;
				} else {
					// find mate of mapped read
					// filter if we cannot find it (mapped to another chromosome, those are a bit suspicious)
                    
					// lookup mapped read
					hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
					if (hash_it == mapped_name_to_idx.end()) { numOrphan++; filter=true; } else {
						if (hash_it->second.size()>2) cerr << "HUH? DUPLICATE READ LABELS???" << endl;
						if (reads[r].mateIsUnmapped() == false) {
							filter = true;
						}
                        
						BOOST_FOREACH(int idx, hash_it->second) {
							if (idx != r) {
								reads[r].mateLen = reads[idx].size();
								reads[r].matePos = reads[idx].pos;
								filter = false;
								if (reads[r].matePos != reads[r].getBAMMatePos()) {
									cerr << "matepos inconsistency!" << endl;
									cerr << reads[r].matePos << " " << reads[r].getBAMMatePos() << endl;
									exit(1);
								}
							}
						}
                        
						if (filter == true) {
							numOrphan++;
							tf = 2;
						}
						//cout << "mapped read: " << r << " " << qname << " pos: " << reads[r].pos << " " << reads[r].getBam()->core.mtid << " " <<  reads[r].getBam()->core.mpos << " mateunmap: " << reads[r].mateIsUnmapped() << endl;
					}
				}
			} else if (reads[r].mateIsUnmapped() == true) {
				reads[r].matePos=reads[r].pos;
				hash_it = unmapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
				if (hash_it == unmapped_name_to_idx.end()) { filter=true; } else {
					filter = true;
					if (hash_it->second.size()>2) cerr << "HUH? DUPLICATE READ LABELS???" << endl;
					BOOST_FOREACH(int idx, hash_it->second) {
						if (idx != r) {
							reads[r].mateLen = reads[idx].size();
							filter = false;
						}
					}
                    
				}
				if (filter==true) {
					numOrphan++;
					tf = 3;
				}
			}
			if (filter == false) numInRegion++;
            
		} else {
            //		cout << " unmapped" << endl;
			// read is unmapped
			if (params.mapUnmappedReads) {
                //			cout << " unmapped " << qname << " ; ";
                
				// lookup mapped read
				hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
				int idx;
				if (hash_it == mapped_name_to_idx.end()) { numOrphanUnmapped++; filter=true; tf = 4;} else {
					if (hash_it->second.size()!=1) {
						cerr << "UNMAPPED READ HAS MORE THAN ONE MATE!" << endl;
						exit(1);
					}
					idx = *(hash_it->second.begin());
					//cout << " FOUND " << idx << endl ;
					// check if mate will overlap with haplotype
					uint32_t range_l, range_r; // range of mate
                    
					int maxInsert = (int) reads[idx].getLibrary().getMaxInsertSize();
					int minInsert = 0;
					uint32_t rpos = reads[idx].pos;
                    
                    //				cout << "idx: " << idx << " unmapped: " << reads[idx].isUnmapped() << " rpos: " << rpos << " isreverse: " << reads[idx].isReverse() << endl;
					if (reads[idx].isReverse()) {
						range_l = rpos-maxInsert;
						range_r = rpos-minInsert;
					} else {
						range_l = rpos+minInsert;
						range_r = rpos+maxInsert;
					}
                    
                    //				cout << "insert: " << insert << " std: " << std << " range_l : " << range_l << " range_r: " << range_r << " leftPos: " << leftPos << "  rightPos: " << rightPos << endl;
                    
					if (range_r>leftPos && range_l<rightPos) {
						numInRegion++;
						filter=false;
						reads[r].mapQual = reads[idx].mapQual;
						reads[r].matePos = reads[idx].pos;
						reads[r].mateLen = reads[idx].size();
						if (reads[r].isReverse() == reads[idx].isReverse()) {
							reads[r].reverse();
							reads[r].complement();
						}
					} else {
						filter=true;
						tf = 5;
					}
				}
			} else {
				filter = true;
			}
            
		}
		if (filter == true) reads[r].mapQual = -1.0;
        
        //		cout << "reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << " Filter: " << tf << " filter: " << filter <<  " mq: " << reads[r].mapQual << endl;
	}
    
    
	int nUnmapped = 0;
	int nMateposError = 0;
	sort(reads.begin(), reads.end(), SortFunc::sortFunc);
	size_t max; for (max=0;max<params.maxReads && max<reads.size();max++) if (!(reads[max].mapQual<minMapQual)) {
		if (reads[max].matePos==-1 && reads[max].isPaired() && !reads[max].mateIsUnmapped() ) {
			nMateposError++;
			reads[max].matePos = reads[max].pos;
		};
		filteredReads.push_back(Read(reads[max]));
		if (reads[max].isUnmapped()) nUnmapped++;
	} else break;
    
    
	filteredReads.swap(reads);
	filteredReads.clear();
	//Read::filterReads(reads, params.maxReads, params.mapQualThreshold, params.maxReadLength, params.minReadOverlap, leftPos, rightPos);
    
	if (params.filterReadAux!="") {
		if (params.filterReadAux.size()>1) {
			size_t before=reads.size();
			int exclude=1;
			if (params.filterReadAux[0]=='+') exclude=0;
			string match=params.filterReadAux.substr(1,params.filterReadAux.size());
			Read::filterReads(reads, exclude, match);
			size_t after=reads.size();
			if (!params.quiet) cout << "filterAux: " << before-after << " reads were filtered based on match string " << params.filterReadAux << endl;
		}
	}
    
	if (!params.quiet) cout << "Number of reads: " << reads.size() << " out of " << oldNumReads << " # unmapped reads: " << nUnmapped << " numReadsUnknownLib: " << numUnknownLib << " numChrMismatch: " << numTIDmismatch << " numMappedWithoutMate: " << numOrphan << " numUnmappedWithoutMate: " << numOrphanUnmapped << endl;
	if (nMateposError) {
		cerr << "The mate position of " << nMateposError << " reads was recorded as -1 in the BAM file" << endl;
	}
    
	if (params.showReads) {
		for (size_t r=0;r<reads.size();r++) {
			cout << "read[" << r << "]: " << reads[r] << endl;
		}
	}
    
	if (reads.size()<2) {
		throw string("too_few_reads");
	} else if (reads.size()>=params.maxReads) {
		throw string("above_read_count_threshold");
	}
    
}

string DetInDel::get_var_symbol(string ref, string obs) {
    if(ref == "-") return("+" + obs);
    else if(obs == "-") return("-" + ref);
    else return(ref + "=>" + obs);
}

vector<AlignedVariant> DetInDel::parse_close_vars(string s) {
    vector<AlignedVariant> variants;
    if(s == "-") return variants;
    vector<string> str_vars;
    vector<string> pos_var;
    split(str_vars, s, is_any_of(","));
    BOOST_FOREACH(string x, str_vars) {
        split(pos_var, x, is_any_of(":"));
        AlignedVariant v(pos_var[1], atoi(pos_var[0].c_str()), -1, 0);
        variants.push_back(v);
    }
    return variants;
}

AlignedCandidates DetInDel::getCandidateVariants(string line, vector<AlignedVariant>& close_somatic, vector<AlignedVariant>& close_germline) {
    int leftPos, rightPos;
    stringstream  linestream(line);
    string cand_somatic, cand_germline;
    VariantInfo vi;
    linestream >> vi.chr >> vi.start >> vi.end >> vi.ref >> vi.obs >> vi.ref_count_tumor >> vi.obs_count_tumor >> vi.ref_count_normal >> vi.obs_count_normal >> vi.missrate_tumor >> vi.strandrate_tumor >> vi.missrate_normal >> vi.strandrate_normal >> vi.fisher_score >> cand_somatic >> cand_germline;

    leftPos = vi.start - 90;//window size
    rightPos = vi.start + 90;
    vector<AlignedVariant> variants;
    AlignedVariant variant(get_var_symbol(vi.ref, vi.obs), vi.start, -1, 0);
    variant.info = vi;
    variants.push_back(variant);
    close_somatic = parse_close_vars(cand_somatic);
    close_germline = parse_close_vars(cand_germline);
    return AlignedCandidates(vi.chr, variants, leftPos, rightPos);
}

void DetInDel::mutationCall(const string & variantsFileName)
{
	cout << "mutationCall " << variantsFileName;
	
    ofstream output;
    ofstream glfOutput;
    string callsFile=params.fileName; callsFile.append(".calls.txt");
	string glfFile=params.fileName; glfFile.append(".glf.txt");
    std::ofstream ofs((params.fileName+".haplotypes.txt").c_str());
    ofs << "haplotypes" << endl;
    std::ofstream ofs_lb((params.fileName+".lower_bounds.txt").c_str());
	OutputData oData=params.makeMutationData(output);
    output.open(callsFile.c_str());
    if (!output.is_open()) {
        throw(string("Cannot open file ").append(callsFile).append(" for writing."));
    }
    oData.outputLine(oData.headerString());
	glfOutput.open(glfFile.c_str());
	if (!glfOutput.is_open()) {
		throw(string("Cannot open file ").append(glfFile).append(" for writing."));
	}
	OutputData glfData=params.makeGLFOutputData(glfOutput);
	glfData.outputLine(glfData.headerString());
    
	ifstream varfile(variantsFileName.c_str());
    string line;
	int index=0;
	//for (map<uint32_t,InDel>::const_iterator it=indels.begin();it!=indels.end();it++, cnt++) {
    
	vector<Read *> normalReadBuf;
    vector<Read *> tumorReadBuf;
	uint32_t oldLeftPosForN=0, oldRightFetchReadPosForN=0;
    uint32_t oldLeftPosForT=0, oldRightFetchReadPosForT=0;
    
	string oldTid("-1");
    
	// NOTE ReadBuffer should be reset on first usage or on chromosome change!
	bool resetReadBuffer = true;
    
	while(getline(varfile, line)) {
        vector<AlignedVariant> close_somatic_vars, close_germline_vars;
		AlignedCandidates candidateVariants = getCandidateVariants(line, close_somatic_vars, close_germline_vars);
		if (candidateVariants.variants.size()==0) continue;
        cout << "for each candidate variants(size=" << candidateVariants.variants.size() << ")" << endl;
		vector<Read> normalReads;
        vector<Read> tumorReads;
        vector<Read> mergedReads;
		uint32_t pos, leftPos, rightPos;
		// get lowest and highest position
        leftPos = candidateVariants.leftPos;
		rightPos = candidateVariants.rightPos;
        
		pos = candidateVariants.centerPos;
		params.tid=candidateVariants.tid;
        
		if (params.tid!=oldTid) {
			// reinit
			resetReadBuffer = true;
			oldTid = params.tid;
			oldLeftPosForT = 0;
            oldLeftPosForN = 0;
		}
        ///oldがnormalのみになってる。
		if (leftPos < oldLeftPosForN) {
			cerr << "leftPos: " << leftPos << " oldLeftPosForN: " << oldLeftPosForN << endl;
			cerr << "Candidate variant files must be sorted on left position of window!" << endl;
			exit(1);
		}
        
        
        // TODO either add tid to AlignedVariant or infer it from the vector of aligned variants
        // change alige
		index++;
		bool skipped = false;
        
		if (!params.quiet) cout << "****" << endl << " tid: " << params.tid << " pos: " << pos << " leftPos: " << leftPos << " " << " rightPos: " << rightPos << endl;
        
		string message="ok";
		try {
            normalReads.clear();
            tumorReads.clear();
			getReadsFromBams(normalBams, leftPos, rightPos, normalReads, oldLeftPosForN, oldRightFetchReadPosForN, normalReadBuf, resetReadBuffer);
            getReadsFromBams(tumorBams, leftPos, rightPos, tumorReads, oldLeftPosForT, oldRightFetchReadPosForT, tumorReadBuf, resetReadBuffer);
            cout << "read size: n, t =" << normalReads.size() << " " << tumorReads.size() << endl;
            uint32_t rs=(int(leftPos)>params.minReadOverlap)?(leftPos-params.minReadOverlap):0;
            uint32_t re=rightPos+params.minReadOverlap;
            string refSeq=getRefSeq(rs+1, re+1);
            string refSeqForAlign=getRefSeq(leftPos+1, rightPos+1);
            MutationCall::computeBayesFactors(index, normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, oData, glfData, params, refSeq, refSeqForAlign, close_somatic_vars, close_germline_vars);
            cout << "after read size: n, t =" << normalReads.size() << " " << tumorReads.size() << endl;
            
		}
		catch (string s) {
			for (size_t x=0;x<s.size();x++) if (s[x]==' ') s[x]='_';
			message=string("error_").append(s);
			skipped=true;
			goto _end;
		}
		catch (std::bad_alloc) {
			message = string("error_bad_alloc");
			skipped = true;
			goto _end;
		}
		catch (std::exception& e) {
			message = string("error_exception_").append(e.what());
			skipped = true;
			goto _end;
		}
        
    _end:
        
		if (skipped) {
			cerr << "skipped " << params.tid << " " << pos << " reason: " << message << endl;
			//OutputData::Line line(oData);
			//line.set("msg", message);
			//line.set("index", index);
			//oData.output(line);
            
			OutputData::Line gline(glfData);
			gline.set("msg", message);
			gline.set("index", index);
			gline.set("tid", params.tid);
			gline.set("lpos", leftPos);
			gline.set("rpos", rightPos);
			glfData.output(gline);
            
			// reset read buffer: all reads will be fetched again
			resetReadBuffer = true;
		} else {
			resetReadBuffer = false;
		}
        
		oldLeftPosForT = leftPos;
        oldLeftPosForN = leftPos;
	}
    
	//output.close();
	glfOutput.close();
    
	// clean up read buffer
	for (size_t r=0;r<normalReadBuf.size();r++) {
		if (normalReadBuf[r]!=NULL) delete normalReadBuf[r];
	}
    for (size_t r=0;r<tumorReadBuf.size();r++) {
		if (tumorReadBuf[r]!=NULL) delete tumorReadBuf[r];
	}
    
    
}




bool DetInDel::alignHaplotypes(vector<Haplotype> & haps,  uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants)
{
    cout << "alignHaplotypes" << endl;
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
        
		string refSeq=getRefSeq(start+1, end+1);
        
        
        
        
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
    /*   cout << "#haplotype list[debug4]" << endl;
     for (size_t th=0;th<haps.size();th++) {
     const Haplotype & hap=haps[th];
     cout << "hap[" << th << "] " << hap.seq << endl;
     for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
     if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
     cout << "[" << it->second.getString() << " " << (it->first) << "]";
     }
     }
     cout << endl;
     }*/
    
    
	// add REF allele as variant to each haplotype in order to compute coverage statistics
	for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
		for (size_t h=0;h<haps.size();h++) haps[h].addRefVariant(it->first);
	}
    
    /*  cout << "#haplotype list[debug5]" << endl;
     for (size_t th=0;th<haps.size();th++) {
     const Haplotype & hap=haps[th];
     cout << "hap[" << th << "] " << hap.seq << endl;
     for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
     if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
     cout << "[" << it->second.getString() << " " << (it->first) << "]";
     }
     }
     cout << endl;
     }*/
    
	if (!params.quiet) {
		for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
			cout << "aligned_var@pos " << pos << " " << leftPos+it->first;
			BOOST_FOREACH(AlignedVariant av, it->second) {
				cout << " " << av.getString();
			}
			cout << endl;
		}
	}
    
    
	return true;
}


void DetInDel::computeHapPosition(const Haplotype & hap, const Read & read, vector<int> & alPos, int leftPos)
{
	// get position on haplotype of read alignment to reference from the aligned position of first and last base in the read
    
	const bam1_t *b=read.getBam();
	const bam1_core_t *c=&b->core;
	uint32_t* cigar=bam1_cigar(b);
	int k, end, start;
	end = c->pos;
    
	int offs=0, l=0, lend; // offset due to SOFT_SKIP at beginning
    // lend is base for which end is computed (there might be a SOFT_CLIP at the end of the read)
    
	bool al=false;
	for (k = 0; k < (int) c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
        
		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CREF_SKIP) al=true;
        
		if (!al && op == BAM_CSOFT_CLIP) offs += cigar[k] >> BAM_CIGAR_SHIFT;
        
		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP)
            l += cigar[k] >> BAM_CIGAR_SHIFT;
        
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP) {
			end += cigar[k] >> BAM_CIGAR_SHIFT;
			lend=l;
		}
	}
    
	start=c->pos-leftPos; // make relative to alignment of haplotypes to reference
	end-=leftPos;
    
	// lookup start and end in alignment of haplotype to reference
    
	for (int x=0;x<int(hap.ml.hpos.size());x++) if (hap.ml.hpos[x]==start) {
		alPos.push_back(hap.ml.hpos[x]-offs);
		break;
	}
    
	for (int x=int(hap.ml.hpos.size())-1;x>=0;x--) if (hap.ml.hpos[x]==end) {
		alPos.push_back(hap.ml.hpos[x]-lend);
		break;
	}
    
    
}

void DetInDel::computeLikelihoodsFaster(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap)
{
	liks.clear();
	liks=vector<vector<MLAlignment> >(haps.size(),vector<MLAlignment>(reads.size()));
	onHap = vector<int>(reads.size(),0); // records whether a read was aligned onto at least one haplotype
    
	const unsigned int kmer=4;
    
	for (size_t h=0;h<haps.size();h++) {
		const Haplotype & hap = haps[h];
		//cout << "Haplotype[" << h << "]: " << endl;
		HapHash hash(kmer, hap);
		for (size_t r=0;r<reads.size();r++) {
			// given BWA alignment of read to reference, estimate a number of likely alignments to the haplotype
			vector<int> alPos;
			computeHapPosition(hap, reads[r], alPos, leftPos);
            //	cout << "[" << r << ": " << alPos.size() ;
            
            
			ObservationModelS om(hap, reads[r], leftPos, params.obsParams);
            
			// align using guessed alignments and the haplotype hash
			liks[h][r]=om.align(hash);
            //	cout << "," << liks[h][r].ll << "] ";
			if (!liks[h][r].offHapHMQ) onHap[r]=1; // if on-haplotype with artificial high mapping quality
            
			/*
             seqan::Score<int> score(-1, -460, -100,-960);
             ObservationModelSeqAn om2(hap, reads[r], leftPos, params.obsParams, score);
             om2.align();
             */
            
 		}
		//cout << "done" << endl;
	}
    
	// check HMQ off-haplotype reads
    
	// realign a couple of high-mapping quality reads to obtain new candidate indels
	// propose new set of candidate haplotypes by realigning all reads to the new set of candidate haplotypes
    
    
	// --bamFiles /Users/caa/Documents/workspace/DInDelFastProb/bamfiles.txt --output test --region 12036338-12036340 --maxReadLength 60 --tid 17 --maxHap 8 --showEmpirical  --minReadOverlap 20 --width 60 --maxLengthIndel 10 --ref /Users/caa/data/human_b36_female.Y.fa --pError 0.0005  --maxRead 2000 --computeMAP
}

double DetInDel::getPairPrior(const AlignedVariant & av1, const AlignedVariant & av2, int leftPos, const AlignedCandidates & candidateVariants)
{
	std::set<AlignedVariant> vars;  vars.insert(av1); vars.insert(av2);
	double ll = 0.0;
	BOOST_FOREACH(AlignedVariant avar, vars) {
		double lnf = 0.0;
		int type = avar.getType();
		const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
        
		if (type==Variant::SNP) lnf = log(params.priorSNP); else if (type==Variant::DEL || type==Variant::INS) lnf = log(params.priorIndel);
		if (av==NULL) {
			ll += lnf;
		} else {
			double prior = av->getFreq();
			if (prior<0.0) ll += lnf; else ll+=log(prior);
		}
	}
    
	return ll;
    
}

double DetInDel::getHaplotypePrior(const Haplotype & h1, const Haplotype & h2, int leftPos, const AlignedCandidates & candidateVariants)
{
	// returns log prior probability of the pair of haplotypes, based on settings in params
	// one day maybe change prior for known SNPs
	double ll=0.0;
    
	// count indels
	typedef  map <int, AlignedVariant>::const_iterator AVIt;
	//hash_map <int, int> indels, snps;
	std::set <AlignedVariant> indels, snps;
	for (AVIt it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>")==string::npos ) {
        //indels[it->first]=1;
        indels.insert(it->second);
	}
	for (AVIt it=h2.indels.begin();it!=h2.indels.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>")==string::npos ) {
		//indels[it->first]=1;
		indels.insert(it->second);
	}
    
	for (AVIt it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>D")==string::npos) {
		snps.insert(it->second);
		//snps[it->first]=1;
	}
	for (AVIt it=h2.snps.begin();it!=h2.snps.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>D")==string::npos) {
		snps.insert(it->second);
		//snps[it->first]=1;
	}
    
	BOOST_FOREACH(AlignedVariant indel, indels) {
        //	cout << "indel: " << indel.getString() << " " << indel.getStartHap();
		const AlignedVariant *av = candidateVariants.findVariant(indel.getStartHap()+leftPos, indel.getType(), indel.getString());
		if (av==NULL) {
            //	cout << " not found. " << endl;
			ll += log(params.priorIndel);
		} else {
			double prior = av->getFreq();
            //	cout << " found: " << prior << " " << av->getStartHap() << " " << av->getFreq() << endl;
			if (prior<0.0) ll += log(params.priorIndel); else ll+=log(prior);
		}
	}
	BOOST_FOREACH(AlignedVariant indel, snps) {
        //	cout << "snp: " << indel.getString() << " " << indel.getStartHap();
		const AlignedVariant *av = candidateVariants.findVariant(indel.getStartHap()+leftPos, indel.getType(), indel.getString());
		if (av==NULL) {
            //	cout << " not found. " << endl;
			ll += log(params.priorIndel);
		} else {
			double prior = av->getFreq();
            //	cout << " found: " << prior << " " << av->getStartHap() << " " << av->getFreq() << endl;
			if (prior<0.0) ll += log(params.priorIndel); else ll+=log(prior);
		}
	}
	/*
     BOOST_FOREACH(AlignedVariant snp, snps) {
     const AlignedVariant *av = candidateVariants.findVariant(snp.getStartHap()+leftPos, snp.getString());
     if (av==NULL) ll += log(params.priorSNP); else {
     double prior = av->getFreq();
     if (prior<0.0) ll += log(params.priorSNP); else ll+=log(prior);
     }
     }
     */
	/*
     int numIndels=int(indels.size());
     int numSNPs=int(snps.size());
     
     ll+=double(numIndels)*log(params.priorIndel);
     ll+=double(numSNPs)*log(params.priorSNP);
     */
    //	cout << "ll: " << ll << endl;
	return ll;
}
#include <fstream>


void DetInDel::computeAnotherLowerBound(const vector<double> & resp, const double a0, const vector<double> & ak, const vector<double> & ln_p_x_given_h, const vector<double> & lpi, const vector<int> & compatible, lower_bound_t & lb) {
    cout << endl << "###computeAnotherLowerBound###" << endl;
    int nh = ak.size(); //number of haps 
    int nr = resp.size() / nh; //number of reads
    cout << "###compatible" << endl;
    for(int i=0;i < nh;i++) cout << " " << compatible[i];
    cout << endl;
    cout << "lpi = "; 
    for(int i=0;i < nh;i++) cout << " " << lpi[i];
    for(int i=0;i < nh;i++) cout << " " << ak[i];
    cout << endl;
    //ln_p_x_given_z
    double tmp = 0.0;
    for(int i=0;i < nr;i++) {
        for(int j=0;j < nh;j++) {
            if(!compatible[j]) continue;
            tmp += resp[i*nh+j] * ln_p_x_given_h[i*nh+j];
        }
    }
    lb.ln_p_z_given_pi = tmp;
    //entropy
    double entropy = 0.0;
    for(int i=0;i < nr;i++) {
        for(int j=0;j < nh;j++) {
            if(!compatible[j]) continue;
            double t = resp[i*nh+j] * log(resp[i*nh+j]); 
            if(!isnan(t)) entropy -= t; //check NaN
        }
    }
    //ln_c
    tmp = 0.0;
    //lgamma(sum(a0)) - sum(lgamma(a0))
    double ahat = 0.0;
    for(size_t h=0;h<nh;h++) if (compatible[h]) {
        ahat += a0;
        tmp += lgamma(a0);
    }
    double ln_c_a0 = lgamma(ahat) - tmp;
    tmp = 0.0;
    ahat = 0.0;
    for(size_t h=0;h<nh;h++) if (compatible[h]) {
        ahat += ak[h];
        tmp += lgamma(ak[h]);
    }
    double ln_c_a = lgamma(ahat) - tmp;
    lb.lower_bound = lb.ln_p_z_given_pi + entropy + ln_c_a0 - ln_c_a;
    cout << "another lb=" << lb.lower_bound << ", lb.ln_p_z_given_pi=" << lb.ln_p_z_given_pi << ", entropy=" << entropy << ", ln_c_a0(" << ln_c_a0 << ") - ln_c_a(" << ln_c_a << ")=" << ln_c_a0 - ln_c_a << endl;
}


double DetInDel::computeModelEvidence(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, const vector<double> & hapFreqs) {
    double me = 0.0;
    for(int n=0;n<reads.size();n++) {
        double l = 0.0;
        // cout << "#" << n << " ";
        for(int k=0;k<haps.size();k++) {
            //   cout << "[" << hapFreqs[k] << " " << exp(liks[k][n].ll) << "]";
            l += hapFreqs[k] * exp(liks[k][n].ll);
        }
        double t = log(l);
        if(!isnan(t)) {
            me += t;
        }
        //cout << " " << t << endl;
    }
    //cout << endl;
    return me;
}


void DetInDel::estimateHaplotypeFrequencies(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs)
{
    
	// estimate haplotype frequencies using EM
	hapFreqs.clear();
    
	size_t nh=haps.size();
	size_t nr=reads.size();
    
    
    
    
	vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods
    
	vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
	vector<double> pi(nh); // log of haplotype frequencies
	vector<double> nk(nh,0.0); // counts for each haplotype
    
	hapFreqs=nk;
    
	// initialize frequencies
	for (size_t h=0;h<nh;h++) pi[h]=1.0/double(nh);
    
	for (size_t h=0;h<nh;h++) for (size_t r=0;r<nr;r++) {
		// initialize read-haplotype likelihoods
		rl[h*nr+r]=liks[h][r].ll;
        
		// initialize expectations of indicator variables
		z[h*nr+r]=0.5;
	}
    
    
	bool converged=false;
	double tol=params.EMtol;
    
	double eOld=-HUGE_VAL, eNew;
    
	cout << "EM HapFreqs:";
    
	int iter=0;
	while (!converged) {
        
		// compute expectation of indicator variables
		for (size_t h=0;h<nh;h++) nk[h]=0.0;
        
		int idx=0;
		for (size_t r=0;r<nr;r++) {
			double lognorm=-HUGE_VAL;
			// compute responsibilities
			for (size_t h=0;h<nh;h++) {
				z[h*nr+r]=pi[h]+(rl[h*nr+r]);
				lognorm=addLogs(lognorm, z[h*nr+r]);
			}
			// normalize and compute counts
			for (size_t h=0;h<nh;h++) {
				z[nr*h+r]-=lognorm;
				z[nr*h+r]=exp(z[nr*h+r]);
                
				nk[h]+=z[nr*h+r];
			}
		}
        
		// compute frequencies
        
		for (size_t h=0;h<nh;h++) {
			double nph=nk[h]/nr;
			pi[h]=log(nph);
		}
        
        
		idx=0;
		eNew=0.0;
		for (size_t h=0;h<nh;h++) {
            
            for (size_t r=0;r<nr;r++) {
                // compute responsibilities
				eNew+=z[idx]*( pi[h]+rl[idx]);
				idx++;
			}
		}
		//cout << "[" << eNew << "]" << endl;
		//
		if (eOld>eNew) throw string("EM Error in estimateHapFreqs");
		converged=(fabs(eOld-eNew))<tol || iter>25;
        
		eOld=eNew;
        
        
		iter++;
	}
    
	for (size_t h=0;h<nh;h++) { cout << " " << exp(pi[h]); }
	cout << endl;
    
	// output haplotypes and estimated frequencies
    
	for (size_t h=0;h<nh;h++) hapFreqs[h]=exp(pi[h]);
}


void DetInDel::computePairLikelihoods(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<HapPairLik> & likPairs, bool usePrior, const AlignedCandidates & candidateVariants, int leftPos)
{
    //	cout << "Computing pair likelihoods...\n";
	likPairs.clear();
	size_t lh=haps.size();
	likPairs.reserve(lh*(lh/2));
	//const double log10=log(10);
	double maxLL=-HUGE_VAL; int hpm, hmm;
	size_t midx;
	for (size_t hp=0;hp<lh;hp++) for (size_t hm=hp;hm<lh;hm++) {
		double ll=0.0;
		HapPairLik hpl;
		hpl.numIndFirst=0;
		hpl.numIndSecond=0;
		hpl.numOffBoth=0;
		hpl.numOffBothError=0.0;
		hpl.numFirst=0;
		hpl.numSecond=0;
		hpl.h1=hp;
		hpl.h2=hm;
		for (size_t r=0;r<reads.size();r++) {
			ll+=(addLogs(liks[hp][r].ll,liks[hm][r].ll)+log(.5));
			const MLAlignment & ml1=liks[hp][r];
			const MLAlignment & ml2=liks[hm][r];
            
			// record which indel was on which haplotype
			if (!ml1.offHap && ml1.ll>ml2.ll && ml1.indels.size()!=0) hpl.numIndFirst++;
			else if (!ml2.offHap && ml2.ll>ml1.ll && ml2.indels.size()!=0) hpl.numIndSecond++;
			if (ml1.offHapHMQ && ml2.offHapHMQ) { hpl.numOffBoth++; hpl.numOffBothError+=reads[r].mapQual; };
			if (ml1.ll>=ml2.ll) hpl.numFirst++;
			if (ml2.ll>=ml1.ll) hpl.numSecond++;
            
            //		determine read coverage of all the variants
            //		TODO this part is really slow!
			for (map<int, bool>::const_iterator it=ml1.hapIndelCovered.begin();it!=ml1.hapIndelCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage1[it->first].nr++; else hpl.hapIndelCoverage1[it->first].nf++; }
			for (map<int, bool>::const_iterator it=ml2.hapIndelCovered.begin();it!=ml2.hapIndelCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage2[it->first].nr++; else hpl.hapIndelCoverage2[it->first].nf++; }
			for (map<int, bool>::const_iterator it=ml1.hapSNPCovered.begin();it!=ml1.hapSNPCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage1[it->first].nr++; else hpl.hapSNPCoverage1[it->first].nf++; }
			for (map<int, bool>::const_iterator it=ml2.hapSNPCovered.begin();it!=ml2.hapSNPCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage2[it->first].nr++; else hpl.hapSNPCoverage2[it->first].nf++; }
            
            
		}
		if (usePrior) ll+=getHaplotypePrior(haps[hp], haps[hm], leftPos, candidateVariants);
        
		hpl.ll=ll;
		if (ll>maxLL) {
			maxLL=ll;
			hpm=hp;
			hmm=hm;
			midx=likPairs.size();
		}
        //	cout << "hp: " << hp << " hm: " << hm << " ll: " << ll << endl;
		likPairs.push_back(hpl);
	}
    
    
    //	cout << "ML hap: " << hpm << " " << hmm << " midx: " << midx << endl;
	/*
     cout << haps[hpm] << endl;
     cout << "indels: "; for (map<int, AlignedVariant>::const_iterator it=haps[hpm].indels.begin();it!=haps[hpm].indels.end();it++) {
     cout << "[" << it->first << "," << it->second.getString() << "]";
     }
     cout << endl;
     cout << "coverage: ";
     for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage1.begin();it!=likPairs[midx].hapIndelCoverage1.end();it++) {
     cout << "[" << it->first << "," << it->second.nf << "," << it->second.nr << "]";
     }
     cout << endl;
     cout << haps[hmm] << endl;
     cout << "indels: "; for (map<int, AlignedVariant>::const_iterator it=haps[hmm].indels.begin();it!=haps[hmm].indels.end();it++) {
     cout << "[" << it->first << "," << it->second.getString() << "]";
     }
     cout << endl;
     for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage2.begin();it!=likPairs[midx].hapIndelCoverage2.end();it++) {
     cout << "[" << it->first << "," << it->second.nf << "," << it->second.nr << "]";
     }
     
     cout << endl;
     */
    
	class SortFunc {
	public:
		static bool sortFunc(const HapPairLik & hpl1, const HapPairLik & hpl2)
		{
			// sort in decreasing order
			return hpl1.ll>hpl2.ll;
		}
	};
	sort(likPairs.begin(), likPairs.end(),SortFunc::sortFunc);
}

void DetInDel::statisticsHaplotypePair(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, HapPairLik & hpl, OutputData::Line & line)
{
	hpl.numIndFirst=0;
	hpl.numIndSecond=0;
	hpl.numOffBoth=0;
	hpl.numOffBothError=0.0;
	hpl.numFirst=0;
	hpl.numSecond=0;
	int hp=hpl.h1;
	int hm=hpl.h2;
	for (size_t r=0;r<reads.size();r++) {
		const MLAlignment & ml1=liks[hp][r];
		const MLAlignment & ml2=liks[hm][r];
        
		if (!ml1.offHap && ml1.ll>ml2.ll && ml1.indels.size()!=0) hpl.numIndFirst++;
		else if (!ml2.offHap && ml2.ll>ml1.ll && ml2.indels.size()!=0) hpl.numIndSecond++;
		if (ml1.offHapHMQ && ml2.offHapHMQ) { hpl.numOffBoth++; hpl.numOffBothError+=reads[r].mapQual; };
		if (ml1.ll>=ml2.ll) hpl.numFirst++;
		if (ml2.ll>=ml1.ll) hpl.numSecond++;
        
        //		determine read coverage of all the variants
        //		TODO this part is really slow!
		for (map<int, bool>::const_iterator it=ml1.hapIndelCovered.begin();it!=ml1.hapIndelCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage1[it->first].nr++; else hpl.hapIndelCoverage1[it->first].nf++; }
		for (map<int, bool>::const_iterator it=ml2.hapIndelCovered.begin();it!=ml2.hapIndelCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage2[it->first].nr++; else hpl.hapIndelCoverage2[it->first].nf++; }
		for (map<int, bool>::const_iterator it=ml1.hapSNPCovered.begin();it!=ml1.hapSNPCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage1[it->first].nr++; else hpl.hapSNPCoverage1[it->first].nf++; }
		for (map<int, bool>::const_iterator it=ml2.hapSNPCovered.begin();it!=ml2.hapSNPCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage2[it->first].nr++; else hpl.hapSNPCoverage2[it->first].nf++; }
        
		// record which indel was on which haplotype
	}
    
	line.set("num_off_hap", hpl.numOffBoth);
	line.set("num_mapped_to_first",hpl.numFirst);
	line.set("num_mapped_to_second",hpl.numSecond);
}



void parseRegionString(const string & region, int & start, int & end)
{
	string filtered;
	for(size_t x=0;x<region.size();x++) {
		char c=region[x];
		if (c=='-') filtered+=' ';
		else if (c!=',') filtered+=c;
	}
	istringstream is(filtered);
	string e; is >> e;
	if (!from_string(start,e,std::dec)) throw string("Cannot parse region for start!");
	is >> e;
	if (!from_string(end,e,std::dec)) throw string("Cannot parse region end!");
}

void getParameters(po::variables_map & vm, Parameters & params)
{
	params.maxHap=vm["maxHap"].as<uint32_t> ();
	params.maxReads=vm["maxRead"].as<uint32_t> ();
	params.width=vm["width"].as<uint32_t> ();
	params.mapQualThreshold=vm["mapQualThreshold"].as<double>();
	params.skipMaxHap=vm["skipMaxHap"].as<uint32_t>();
	//params.glfNumHap=vm["glfNumHap"].as<uint32_t>();
	params.inferenceMethod=vm["inferenceMethod"].as<string>();
	params.minReadOverlap=vm["minReadOverlap"].as<uint32_t>();
	params.maxReadLength=vm["maxReadLength"].as<uint32_t>();
	//params.scaleErr=vm["mapScaleError"].as<double>();
	//params.minCount=vm["minCount"].as<uint32_t>();
	params.maxHapReadProd=vm["maxHapReadProd"].as<uint32_t>();
    
	params.priorSNP=vm["priorSNP"].as<double>();
	params.priorIndel=vm["priorIndel"].as<double>();
    params.a0[0]=vm["a0"].as<double>();
    params.a0[1]=vm["a1"].as<double>();
    params.a0[2]=vm["a2"].as<double>();
    params.b0[0]=vm["b0"].as<double>();
    params.b0[1]=vm["b1"].as<double>();
    params.c0[0]=vm["c0"].as<double>();
    params.c0[1]=vm["c1"].as<double>();
	params.bayesType=vm["bayesType"].as<string>();

	if (vm.count("ref")) {
		params.alignAgainstReference=true;
		params.refFileName=vm["ref"].as<string>();
	} else {
		params.alignAgainstReference=false;
	}
    
	params.obsParams.pError=vm["pError"].as<double>();
	params.obsParams.pMut=vm["pMut"].as<double>();
    
	//params.obsParams.baseQualThreshold=vm["baseQualThreshold"].as<double>();
	//	params.obsParams.fixedBaseQual=vm["fixedBaseQual"].as<double>();
	params.obsParams.maxLengthIndel=vm["maxLengthIndel"].as<int>();
	params.obsParams.maxLengthDel=params.obsParams.maxLengthIndel;
	params.obsParams.mapQualThreshold=vm["capMapQualThreshold"].as<double>();
	params.obsParams.capMapQualFast=vm["capMapQualFast"].as<double>();
	//params.obsParams.scaleErr=vm["obsScaleError"].as<double>();
	//params.obsParams.numE=vm["numE"].as<int>();
	params.obsParams.padCover = vm["flankRefSeq"].as<int>();
	params.obsParams.maxMismatch = vm["flankMaxMismatch"].as<int>();
	params.checkAllCIGARs=vm["checkAllCIGARs"].as<int>();
    
	params.varFileIsOneBased=vm.count("varFileIsOneBased")?true:false;
	params.outputRealignedBAM=vm.count("outputRealignedBAM")?true:false;
	params.analyzeLowFreq=vm.count("compareReadHap")?true:false;
	params.analyzeLowFreqDiffThreshold=vm["compareReadHapThreshold"].as<double>();
	params.showHapDist=vm.count("showEmpirical")?true:false;
	params.showCandHap=vm.count("showCandHap")?true:false;
	params.showReads=vm.count("showReads")?true:false;
	params.quiet=vm.count("quiet")?true:false;
	params.computeML=vm.count("computeML")?true:false;
	params.computeMAP=vm.count("computeMAP")?true:false;
	params.doDiploid=vm.count("doDiploid")?true:false;
    
	params.filterHaplotypes=vm.count("filterHaplotypes")?true:false;
    
	params.printCallsOnly=vm.count("printCallsOnly")?true:false;
	params.estimateHapFreqs=vm.count("doPooled")?true:false;
	params.outputPooledLikelihoods=vm.count("opl")?true:false;
	params.showHapAlignments=vm.count("showHapAlignments")?true:false;
	if (vm.count("filterReadAux")) params.filterReadAux=vm["filterReadAux"].as<string>();
	if (vm.count("processRealignedBAM")) params.processRealignedBAM=vm["processRealignedBAM"].as<string>();
    
	params.slower=vm.count("faster")?false:true;
	params.changeINStoN=vm.count("changeINStoN")?true:false;
	params.outputGLF=true;
    
    
	// removed options
    /*
     params.outputRealignedBAM=vm.count("outputRealignedBAM")?true:false;
     params.obsParams.modelType=vm["modelType"].as<string>();
     params.mapUnmappedReads=vm.count("mapUnmapped")?true:false;
     params.obsParams.mapUnmappedReads=vm.count("mapUnmapped")?true:false;
     params.obsParams.pFirstgLO=vm["pFirstgLO"].as<double>();
     //params.numOutputTopHap=vm["numOutputTopHap"].as<int>();
     
     */
    
}



int smain(int argc, char *argv[])
{
	if (1) {
		Haplotype hap;
		Read read;
        
		//hap.seq=          "ATCGTGTAGCTCTCTGGCTGGCTAGCTGATTGGCTCTTGCC";
		//read.seq.seq=              "CTCTCTGGCTGGCTAGCGAT";
		Haplotype ref;
		//                 012345678901234567890123456789
		hap.seq=          "ATCGATTCGTGATATATATATTCAATGTAGTCGCTAG";
		read.seq.seq=     "ATCGATTCGTGATAATATTCAATGTAGTCGCTAG";
        
        
		//hap.seq=          "ATCGATTCGTGATATATATATTCAATGTAGTCGCTAG";
		//read.seq.seq=     "ATCGATTCGTGATATATATATAATTCAATGTAGTCGCTAG";
        
		//                 012345678901234567890123456789
		//hap.seq=          "ATCGATTCGTGTTTTTTCAATGTAGTCGCTAG";
		//read.seq.seq=     "ATCGATTCGTGTTTTTCAATGTAGTCGCTAG";
        
		read.mapQual=1-1e-16;
        
		ObservationModelParameters obsParams;
		read.setAllQual(0.99);
        
		ObservationModelFBMaxErr omfbe(hap, read, 0, obsParams);
		/*
         ObservationModelS oms(hap, read, 0, obsParams);
         HapHash hash(4, hap);
         
         oms.align(hash);
         */
        
		double ll= omfbe.calcLikelihood();
        
		cout << "ll: " << ll << endl;
		cout << string(50,' ') << hap.seq << endl;
		omfbe.printAlignment(50);
	}
	if (0) {
		Haplotype hap;
		Read read;
        
		//hap.seq=          "ATCGTGTAGCTCTCTGGCTGGCTAGCTGATTGGCTCTTGCC";
		//read.seq.seq=              "CTCTCTGGCTGGCTAGCGAT";
        
		hap.seq=          "AAAATCACCAACACTTCATAATCTATTTTTTCCCCTGAGGAACTTCCTAAAATGAATAAAAAAAAACCCCAGCCACATCTGCATTTGCAAACAGGAAACTCTGCAAGCCATACTAAGACCAAAGCTTAGTT";
		read.seq.seq=     "CAAACAGGAAACTCTGCAAGCCATACTAAGACCAAAGCTTAGTTA";
        
        
		read.mapQual=1-1e-16;
        
		ObservationModelParameters obsParams;
		read.setAllQual(0.99);
        
		ObservationModelFBMaxErr omfbe(hap, read, 0, obsParams);
		/*
         ObservationModelS oms(hap, read, 0, obsParams);
         HapHash hash(4, hap);
         
         oms.align(hash);
         */
        
		double ll= omfbe.calcLikelihood();
        
		cout << "ll: " << ll << endl;
		cout << string(50,' ') << hap.seq << endl;
		omfbe.printAlignment(50);
	}
    
	return 0;
    
}






#ifdef DINDEL
int main(int argc, char *argv[])
{
	po::options_description which("[Required] Program option");
	which.add_options()
	("analysis", po::value<string>()->default_value("mutationCall"),"Analysis type:\n"
     "getCIGARindels:  Extract indels from CIGARs of mapped reads, and infer libary insert size distributions\n"
     "mutationCall\n");
    
	po::options_description required("[Required] ");
	required.add_options()
	("ref", po::value<string>(),"fasta reference sequence (should be indexed with .fai file)")
	("outputFile", po::value<string>(),"file-prefix for output results");
    
	po::options_description baminput("[Required] BAM input. Choose one of the following");
	baminput.add_options()
	("bamFile",po::value<string>(), "read alignment file (should be indexed)")
	("bamFiles",po::value<string>(), "file containing filepaths for BAMs to be jointly analysed (not possible for --analysis==indels");
    
    po::options_description bams_tn("[Required for analysis = mutationCall] :");
    bams_tn.add_options()
    ("tumorBamFile", po::value<string>(), "bam files of tumor samples")
    ("normalBamFile", po::value<string>(), "bam files of normal samples");
    
	po::options_description regioninput("[Required for analysis == getCIGARindels]: \nRegion to be considered for extraction of candidate indels.");
	regioninput.add_options()
	("region", po::value<string>(),"region to be analysed in format start-end, eg. 1000-2000")
	("tid", po::value<string>(),"target sequence (eg 'X') ");
    
	po::options_description varfileinput("[Required for analysis == indels, mutationCall]");
	varfileinput.add_options()
	("varFile", po::value<string>(), "file with candidate variants to be tested.")
	("varFileIsOneBased", "coordinates in varFile are one-based");
    
	po::options_description output_options("Output options");
	output_options.add_options()
	("outputRealignedBAM", "output BAM file with realigned reads")
	("processRealignedBAM", po::value<string>(),"ABSOLUTE path to script to process realigned BAM file")
	//("outputGLF", "outputGLF for individuals in each bam file")
	("quiet", "quiet output");
	//("printCallsOnly", "print only genotypes where call_lik_ref>0.0001 (only affects --single)");
    
	po::options_description single_analysis("parameters for analysis==indels option");
	single_analysis.add_options()
	("doDiploid", "analyze data assuming a diploid sequence")
	("doPooled", "estimate haplotype frequencies using Bayesian EM algorithm.\nMay be applied to single individual and pools.");
    
	po::options_description analysis_opt("General algorithm parameters");
	analysis_opt.add_options()
	//("mapUnmapped", "remap unmapped reads for which mate is mapped")
	("faster","use faster but less accurate ungapped read-haplotype alignment model")
	("filterHaplotypes","prefilter haplotypes based on coverage")
	("flankRefSeq",po::value<int>()->default_value(2),"#bases of reference sequence of indel region")
	("flankMaxMismatch",po::value<int>()->default_value(2),"max number of mismatches in indel region")
	("priorSNP", po::value<double>()->default_value(1.0/1000.0), "prior probability of a SNP site")
	("priorIndel", po::value<double>()->default_value(1.0/10000.0), "prior probability of a detected indel not being a sequencing error")
	("width", po::value<uint32_t>()->default_value(80), "number of bases to left and right of indel")
	("maxHap", po::value<uint32_t>()->default_value(8), "maximum number of haplotypes in likelihood computation")
	("maxRead", po::value<uint32_t>()->default_value(50000), "maximum number of reads in likelihood computation")
	("mapQualThreshold", po::value<double>()->default_value(0.99), "lower limit for read mapping quality")
	("capMapQualThreshold", po::value<double>()->default_value(100.0), "upper limit for read mapping quality in observationmodel_old (phred units)")
	("capMapQualFast", po::value<double>()->default_value(45.0), "cap mapping quality in alignment using fast ungapped method\n (WARNING: setting it too high (>50) might result in significant overcalling!)")
	("skipMaxHap", po::value<uint32_t>()->default_value(200), "skip computation if number of haplotypes exceeds this number")
	//("glfNumHap", po::value<uint32_t>()->default_value(5), "number of haplotypes per glf-class")
	//("numOutputTopHap", po::value<int>()->default_value(5), "number of haplotype pairs output to haplotype file")
	("minReadOverlap", po::value<uint32_t>()->default_value(40),"minimum overlap between read and haplotype")
	("maxReadLength", po::value<uint32_t>()->default_value(500),"maximum length of reads")
	("minCount", po::value<uint32_t>()->default_value(1), "minimum number of WS observations of indel")
	("maxHapReadProd",po::value<uint32_t>()->default_value(10000000), "skip if product of number of reads and haplotypes exceeds this value")
	("changeINStoN", "change sequence of inserted sequence to 'N', so that no penalty is incurred if a read mismatches the inserted sequence");
	po::options_description pooled_analysis("parameters for --pooled option");
	pooled_analysis.add_options()
	("bayesa0", po::value<double>()->default_value(0.001), "Dirichlet a0 parameter haplotype frequency prior")
    ("a0", po::value<double>()->default_value(1.0), "beta a0 prior for normal model")
    ("a1", po::value<double>()->default_value(1.0), "beta a0 prior for normal model")
    ("a2", po::value<double>()->default_value(1.0), "beta a0 prior for normal model")
    ("b0", po::value<double>()->default_value(1.0), "beta b0 prior for normal model")
    ("b1", po::value<double>()->default_value(1.0), "beta b0 prior for normal model")
    ("c0", po::value<double>()->default_value(1.0), "beta b0 prior for normal model")
    ("c1", po::value<double>()->default_value(1.0), "beta b0 prior for normal model")
	("bayesType",po::value<string>()->default_value("singlevariant"), "Bayesian EM program type (all or singlevariant or priorpersite)");
    
    
	po::options_description option_filter("General algorithm filtering options");
	option_filter.add_options()
	("checkAllCIGARs",po::value<int>()->default_value(1),"include all indels at the position of the call site")
	("filterReadAux", po::value<string>(), "match string for exclusion of reads based on auxilary information");
    
    
	po::options_description obsModel("Observation model parameters");
	obsModel.add_options()
	("pError", po::value<double>()->default_value(5e-4), "probability of a read indel")
	//("modelType", po::value<string>()->default_value("probabilistic"), "probabilistic/threshold")
	("pMut", po::value<double>()->default_value(1e-5), "probability of a mutation in the read")
	("maxLengthIndel", po::value<int>()->default_value(5), "maximum length of a _sequencing error_ indel in read [not for --faster option]");
	//("pFirstgLO",po::value<double>()->default_value(0.01),"probability of transition from off the haplotype to on the haplotype");
    
	po::options_description libParams("Library options");
	libParams.add_options()
	("libFile", po::value<string>(), "file with library insert histograms (as generated by --analysis getCIGARindels)");
    
    
	po::options_description miscAnalysis("Misc results analysis options");
	miscAnalysis.add_options()
	("compareReadHap",  "compare likelihood differences in reads against haplotypes")
	("compareReadHapThreshold", po::value<double>()->default_value(0.5), "difference threshold for viewing")
	("showEmpirical", "show empirical distribution over nucleotides")
	("showCandHap", "show candidate haplotypes for fast method")
	("showHapAlignments","show for each haplotype which reads map to it")
	("showReads","show reads")
	("inferenceMethod",po::value<string>()->default_value("empirical"), "inference method")
	("opl","output likelihoods for every read and haplotype");
    
	required.add(which).add(baminput).add(regioninput).add(varfileinput).add(output_options).add(single_analysis).add(analysis_opt).add(pooled_analysis).add(option_filter).add(obsModel).add(libParams).add(miscAnalysis).add(bams_tn);
    
	po::variables_map vm;
    
	try {
        po::store(po::parse_command_line(argc, argv, required), vm);
	} catch (boost::program_options::error) {
        cout << "Error parsing input options. Usage: \n\n" << required <<"\n";
        exit(1);
	}
	po::notify(vm);
    
	// analysis
	if (!(vm.count("analysis"))) {
		cerr << "Error: Specify which analysis (--analysis) is required." << endl;
		exit(1);
	}
    
    // required
	if (!(vm.count("ref") && vm.count("outputFile"))) {
		cerr << "Error: One of the following options was not specified:  --ref --tid or --outputFile" << endl;
		exit(1);
	}
    
	if (vm.count("getCIGARindels") && vm.count("region") && !vm.count("tid")) {
		cerr << "--tid must be specified if analysis==getCIGARindels and --region option is used. " << endl;
		exit(1);
	}
#define DEBUGGING
#ifndef DEBUGGING
	try {
#endif
		// extract required parameters
		string file;
		int multipleFiles=0;
		string analysis=vm["analysis"].as<string>();
        
		// baminput
		if (analysis=="indels" || analysis=="getCIGARindels") {
			if (!(vm.count("bamFile") || vm.count("bamFiles"))) {
                cerr << "Error: Specify either --bamFile or --bamFiles." << endl;
                exit(1);
			}
            
			if (vm.count("bamFile")) {
                file=vm["bamFile"].as< string >();
                cout << "Reading BAM file: " << file << endl;
			} else if (vm.count("bamFiles")) {
                file=vm["bamFiles"].as<string>();
                multipleFiles=1;
			}
		}
        
		string faFile=vm["ref"].as<string>();
		string outputFile=vm["outputFile"].as< string >();
        
		string modelType="probabilistic"; //vm["modelType"].as< string >();
		Parameters params(string("1"), outputFile, modelType);
		getParameters(vm, params);
        if (analysis=="mutationCall") {
            cout << "mutationCall" << endl;
            string varFile = vm["varFile"].as<string>();
            string tumorBF = vm["tumorBamFile"].as<string>();
            string normalBF = vm["normalBamFile"].as<string>();
            DetInDel detInDel(normalBF, tumorBF, params);
            if (vm.count("libFile")) {
				cout << "Detected library file..." << endl;
				detInDel.params.mapUnmappedReads=true;
				detInDel.params.obsParams.mapUnmappedReads=true;
				detInDel.addLibrary(vm["libFile"].as<string>());
			}
			detInDel.params.print();
			detInDel.mutationCall(varFile);
        } else if (analysis=="getCIGARindels") {
			GetCandidatesFromCIGAR gcfc;
			string outputFile=vm["outputFile"].as< string >();// outputFile.append(".variants.txt");
			fasta::Fasta fa(faFile);
			if (vm.count("region")) {
				string tid=vm["tid"].as<string>();
				string region=vm["region"].as<string>();
				int start, end;
				parseRegionString(region, start, end);
				DetInDel detInDel(file, params, multipleFiles);
				const vector<MyBam *> &  bams = detInDel.getMyBams();
                
				cout << "Getting indels from CIGARs in mapped reads from region " << tid << ":" << start << "-" << end << endl;
				gcfc.getIndelFromCIGARRegion(bams,tid, start, end, outputFile, fa);
                
			} else {
				if (multipleFiles) {
					cerr << "Can extract the full set of indels from only BAM file at a time." << endl;
					exit(1);
				}
				gcfc.get(file, outputFile, faFile);
			}
		/*} else if (analysis=="indels") {
			if (!vm.count("varFile")) {
				cerr << "Please specify the file with the candidate variants." << endl;
				exit(1);
			}
            
			string varFile = vm["varFile"].as<string>();
            
			DetInDel detInDel(file, params, multipleFiles);
            
			if (vm.count("libFile")) {
				cout << "Detected library file..." << endl;
				detInDel.params.mapUnmappedReads=true;
				detInDel.params.obsParams.mapUnmappedReads=true;
				detInDel.addLibrary(vm["libFile"].as<string>());
			}
			detInDel.params.print();
            
            
			detInDel.detectIndels(varFile);
		} else if (analysis == "realignCandidates") {
			GetCandidatesFromCIGAR gcfc;
			string outputFile=vm["outputFile"].as< string >(); outputFile.append(".variants.txt");
            
			if (!vm.count("varFile")) {
				cerr << "Please specify the file with the candidate variants." << endl;
				exit(1);
			}
            
			string varFile = vm["varFile"].as<string>();
            
			if (varFile == outputFile) {
				cerr << "outputFile is same as variant file used for input!" << endl;
				exit(1);
			}
            
			gcfc.realignCandidateFile(varFile, params.varFileIsOneBased,outputFile, faFile);*/
		} else {
			cerr << "Unrecognized --analysis option." << endl;
			exit(1);
		}
#ifndef DEBUGGING
    }
    catch (string s) {
        cout << "Exception: " << s << endl;
        exit(1);
    }
#endif
    return 0;
}
#endif

