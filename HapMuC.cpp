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


HapMuC::HapMuC(const string & normalBF, const string & tumorBF, const Parameters & _params) : params(_params)
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


HapMuC::~HapMuC()
{
	if (params.alignAgainstReference && fai) {
		fai_destroy(fai);
	}
	for (size_t b=0;b<myBams.size();b++) delete myBams[b];
}


string HapMuC::getRefSeq(uint32_t lpos, uint32_t rpos)
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

void HapMuC::getReadsFromBams(vector<MyBam *> & Bams, uint32_t leftPos, uint32_t rightPos, vector<Read> & reads, uint32_t & oldLeftPos, uint32_t & oldRightFetchReadPos, vector<Read *> & readBuffer, const bool reset)
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
    
    
	int maxDev = int (libraries.getMaxInsertSize());
    
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
    
    
    
    
	 cout << "leftFetchReadPos: " << leftFetchReadPos << " rightFetchReadPos: " << rightFetchReadPos << " oldRightFetchReadPos: " << oldRightFetchReadPos << endl;
	 cout << "leftMostReadPos: " << leftMostReadPos << " rightMostReadPos: " << rightMostReadPos << " leftOverlapsPrevious: " << int(leftOverlapsPrevious) << endl;
	// store updated readBuffer
	readBuffer.swap(newReadBuffer);
    
    	cout << "leftPos : " << leftPos << " rightPos: " << rightPos << " maxDev: " << maxDev << endl;
    
    
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
            cout << b << " " << numReads << endl;
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


string HapMuC::get_var_symbol(string ref, string obs) {
    if(ref == "-") return("+" + obs);
    else if(obs == "-") return("-" + ref);
    else return(ref + "=>" + obs);
}

vector<AlignedVariant> HapMuC::parse_close_vars(string s) {
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

AlignedCandidates HapMuC::getCandidateVariants(string line, vector<AlignedVariant>& close_somatic, vector<AlignedVariant>& close_germline) {
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

void HapMuC::mutationCall(const string & variantsFileName)
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


#ifdef HAPMUC
int main(int argc, char *argv[])
{
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
    
	required.add(baminput).add(varfileinput).add(output_options).add(analysis_opt).add(pooled_analysis).add(option_filter).add(obsModel).add(libParams).add(miscAnalysis).add(bams_tn);
    
	po::variables_map vm;
    
	try {
        po::store(po::parse_command_line(argc, argv, required), vm);
	} catch (boost::program_options::error) {
        cout << "Error parsing input options. Usage: \n\n" << required <<"\n";
        exit(1);
	}
	po::notify(vm);
    
    // required
	if (!(vm.count("ref") && vm.count("outputFile"))) {
		cerr << "Error: One of the following options was not specified:  --ref --tid or --outputFile" << endl;
		exit(1);
	}
#define DEBUGGING
#ifndef DEBUGGING
	try {
#endif
		// extract required parameters
		string file;
		int multipleFiles=0;
        
		string faFile=vm["ref"].as<string>();
		string outputFile=vm["outputFile"].as< string >();
        
		string modelType="probabilistic"; //vm["modelType"].as< string >();
		Parameters params(string("1"), outputFile, modelType);
		getParameters(vm, params);
        
            string varFile = vm["varFile"].as<string>();
            string tumorBF = vm["tumorBamFile"].as<string>();
            string normalBF = vm["normalBamFile"].as<string>();
            HapMuC hapmuc(normalBF, tumorBF, params);
            if (vm.count("libFile")) {
				cout << "Detected library file..." << endl;
				hapmuc.params.mapUnmappedReads=true;
				hapmuc.params.obsParams.mapUnmappedReads=true;
				hapmuc.addLibrary(vm["libFile"].as<string>());
			}
			hapmuc.params.print();
			hapmuc.mutationCall(varFile);
        
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
