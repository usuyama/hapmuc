/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
/*
 * Dindel
 * http://www.sanger.ac.uk/resources/software/dindel/
 *
 * Copyright 2010, Kees Albers
 * Released under the GPL Version 3.
 */
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
#include <boost/foreach.hpp>
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
#include "MutationCall.hpp"
#include "log.h"
#include "cmdline.h"
#include <boost/algorithm/string.hpp>

const int USECALLWINDOW=0;
using namespace seqan;

using namespace std;
//using namespace fasta;


HapMuC::HapMuC(const string & normalBF, const string & tumorBF, const Parameters & _params) : params(_params)
{
	fai=NULL;
	fai = fai_load(params.refFileName.c_str());
	if (!fai) {
		LOG(logERROR) << "Cannot open reference sequence file." << endl;
		exit(1);
	}

    string fname = normalBF;
    if (!fname.empty()) {
        LOG(logINFO) << "Reading BAM file " << fname << endl;
        normalBams.push_back(new MyBam(fname));
        normalBamsFileNames.push_back(fname);
        myBams.push_back(new MyBam(fname));
        myBamsFileNames.push_back(fname);
    }

    fname = tumorBF;
    if (!fname.empty()) {
        LOG(logINFO) << "Reading BAM file " << fname << endl;
        tumorBams.push_back(new MyBam(fname));
        tumorBamsFileNames.push_back(fname);
        myBams.push_back(new MyBam(fname));
        myBamsFileNames.push_back(fname);
    }
}


HapMuC::~HapMuC()
{
	if (fai) {
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
    LOG(logDEBUG) << "getReadsFromBams()" << endl;
    BOOST_FOREACH(MyBam * b, Bams) {
        LOG(logDEBUG) << " " << b->fileName << endl;
    }
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
		LOG(logERROR) << "Windows are not sorted!" << endl;
		exit(3);
	}

	reads.clear();

	if (int(rightPos-leftPos)<params.minReadOverlap) throw string("You should choose a larger maxReadLength or a smaller minReadOverlap.");


	int maxDev = 800; //int (libraries.getMaxInsertSize());

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

    LOG(logDEBUG) << "leftFetchReadPos: " << leftFetchReadPos << " rightFetchReadPos: " << rightFetchReadPos << " oldRightFetchReadPos: " << oldRightFetchReadPos << endl;
    LOG(logDEBUG) << "leftMostReadPos: " << leftMostReadPos << " rightMostReadPos: " << rightMostReadPos << " leftOverlapsPrevious: " << int(leftOverlapsPrevious) << endl;
	// store updated readBuffer
	readBuffer.swap(newReadBuffer);
    LOG(logDEBUG) << "leftPos : " << leftPos << " rightPos: " << rightPos << " maxDev: " << maxDev << endl;

	// first clean readbuffer
	int numReads = readBuffer.size();

	vector<Read> newReads;
	if (leftFetchReadPos<=rightFetchReadPos) {
		LOG(logINFO) << "Fetching reads...." << endl;
		for (size_t b=0;b<Bams.size();b++) {
			Read::FetchReadData data(&newReads, int(b), &Bams, numReads, params.maxReads*100);
			bam_fetch(Bams[b]->bf, Bams[b]->idx, Bams[b]->getTID(params.tid), leftFetchReadPos , rightFetchReadPos,&data , &Read::fetchFuncVectorPooled);
			numReads = data.numReads;
            LOG(logDEBUG) << b << " " << numReads << endl;
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
		} else {
			mapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
		}
	}
	int numTIDmismatch = 0, numOrphan =0, numOrphanUnmapped = 0, numInRegion = 0;

	// reads are filtered by setting mapping quality to -1
	vector<Read> filteredReads;
	double minMapQual = params.mapQualThreshold;
	if (minMapQual<0.0) minMapQual=0.0;
	for (int r=0;r<int(reads.size());r++) {
		string qname = string(bam1_qname(reads[r].getBam()));
        LOG(logDEBUGREADS) << "qname: " << qname << " ";
		bool filter = false;
		int tf = 0;
        LOGP(logDEBUGREADS) << "reads[r].pos: " << reads[r].pos << " reads[r].getEndPos(): " << reads[r].getEndPos() << endl;
		if (reads[r].getEndPos()<leftMostReadPos || reads[r].pos>rightMostReadPos) {
            filter=true;
            LOG(logDEBUGREADS) << "filtered by pos ";
        } else if (reads[r].mapQual < minMapQual) {
            filter = true;
            LOG(logDEBUGREADS) << "filter by mapping quality: " << bam1_qname(reads[r].getBam()) << " mq: " << reads[r].mapQual << endl;
        }
        if (!reads[r].isUnmapped()) {
            LOG(logDEBUGREADS) << endl << "is mapped" << endl;
            if (reads[r].pos+int(reads[r].size())<=int(leftPos)+params.minReadOverlap || reads[r].pos>=int(rightPos)-params.minReadOverlap) {
                LOG(logDEBUGREADS) << "filtered by minOverlap" << endl;
                filter=true;
                tf = 1;
            } else if (reads[r].mateIsUnmapped() == false ){
                LOG(logDEBUGREADS) << "mate is mapped" << endl;
                LOG(logDEBUGREADS) << "lookup mate and filter if we cannot find it (mapped to another chromosome, those are a bit suspicious)" << endl;
                if (reads[r].getBam()->core.mtid != reads[r].getBam()->core.tid) {
                    LOG(logWARNING) << "TIDERR: reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << endl;
                    numTIDmismatch++; filter = true;
                } else {
                    hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
                    if (hash_it == mapped_name_to_idx.end()) {
                        LOG(logDEBUGREADS) << "the mate is not found" << endl;
                        numOrphan++; filter=true;
                    } else {
                        if (hash_it->second.size()>2) LOG(logERROR) << "HUH? DUPLICATE READ LABELS???" << endl;
                        filter = true;
                        LOG(logDEBUGREADS) << "found the mate" << endl;
                        BOOST_FOREACH(int idx, hash_it->second) {
                            if (idx != r) {
                                reads[r].mateLen = reads[idx].size();
                                reads[r].matePos = reads[idx].pos;
                                filter = false;
                                if (reads[r].matePos != reads[r].getBAMMatePos()) {
                                    LOG(logINFO) << "matepos inconsistency!" << endl;
                                    LOG(logINFO) << reads[r].matePos << " " << reads[r].getBAMMatePos() << endl;
                                    exit(1);
                                }
                            }
                        }
                        if (filter == true) {
                            LOG(logWARNING) << "filter: reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " because we could not find the mate read nearby" << endl;
                            numOrphan++;
                            tf = 2;
                        }
                    }
                }
            } else if (reads[r].mateIsUnmapped() == true) {
                LOG(logDEBUGREADS) << "mate is UNmapped; filter" << endl;
                filter = true;
                /*
                 reads[r].matePos=reads[r].pos;
                 hash_it = unmapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
                 if (hash_it == unmapped_name_to_idx.end()) { filter=true; } else {
                 filter = true;
                 if (hash_it->second.size()>2) LOG(logERROR) << "HUH? DUPLICATE READ LABELS???" << endl;
                 BOOST_FOREACH(int idx, hash_it->second) {
                 if (idx != r) {
                 reads[r].mateLen = reads[idx].size();
                 filter = false;
                 }
                 }

                 }
                 */
                if (filter==true) {
                    numOrphan++;
                    tf = 3;
                }
            }
            if (filter == false) numInRegion++;

        } else {
            // read is unmapped
            LOG(logDEBUGREADS) << endl << "is UNmapped; filter" << endl;
            filter = true;
        }
        LOG(logDEBUGREADS) << "reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << " Filter: " << tf << " filter: " << filter <<  " mq: " << reads[r].mapQual << endl;
        if (filter == true) reads[r].mapQual = -1.0;
    }

    int nUnmapped = 0;
    int nMateposError = 0;
    sort(reads.begin(), reads.end(), SortFunc::sortFunc);
    LOG(logDEBUG) << "fetched reads:" << endl;
    size_t max; for (max=0;max<params.maxReads && max<reads.size();max++) if (!(reads[max].mapQual<minMapQual)) {
        if (reads[max].matePos==-1 && reads[max].isPaired() && !reads[max].mateIsUnmapped() ) {
            nMateposError++;
            reads[max].matePos = reads[max].pos;
        };
        LOG(logDEBUG) << "reads[" << max << "]: " << reads[max].seq_name << " pos: " << reads[max].pos << " matePos: " << reads[max].matePos << " mateLen: " << reads[max].mateLen  <<  " mq: " << reads[max].mapQual << endl;
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
			if (!params.quiet) LOG(logINFO) << "filterAux: " << before-after << " reads were filtered based on match string " << params.filterReadAux << endl;
		}
	}

	LOG(logINFO) << "Number of reads: " << reads.size() << " out of " << oldNumReads << " # unmapped reads: " << nUnmapped << " numReadsUnknownLib: " << numUnknownLib << " numChrMismatch: " << numTIDmismatch << " numMappedWithoutMate: " << numOrphan << " numUnmappedWithoutMate: " << numOrphanUnmapped << endl;
	if (nMateposError) {
		LOG(logERROR) << "The mate position of " << nMateposError << " reads was recorded as -1 in the BAM file" << endl;
	}

	if (params.showReads) {
		for (size_t r=0;r<reads.size();r++) {
			LOG(logINFO) << "read[" << r << "]: " << reads[r] << endl;
		}
	}

	if (reads.size()<2) {
        LOG(logERROR) << "size=" << reads.size() << endl;
		throw string("too_few_reads");
	} else if (reads.size()>=params.maxReads) {
        LOG(logERROR) << "size=" << reads.size() << endl;
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
    boost::split(str_vars, s, boost::is_any_of(","));
    BOOST_FOREACH(string x, str_vars) {
        boost::split(pos_var, x, boost::is_any_of(":"));
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
    linestream >> vi.chr >> vi.start >> vi.end >> vi.ref >> vi.obs >> vi.ref_count_tumor >> vi.obs_count_tumor >> vi.ref_count_normal >> vi.obs_count_normal >> vi.missrate_tumor >> vi.strandrate_tumor >> vi.missrate_normal >> vi.strandrate_normal >> vi.ref_bq_tumor >> vi.obs_bq_tumor >> vi.ref_bq_normal >> vi.obs_bq_normal >> vi.fisher_score >> cand_germline;
    leftPos = vi.start - 0;//window size
    rightPos = vi.start + 0;
    vector<AlignedVariant> variants;
    AlignedVariant variant(get_var_symbol(vi.ref, vi.obs), vi.start, -1, 0);
    variant.info = vi;
    variants.push_back(variant);
    close_germline = parse_close_vars(cand_germline);
    // reads must overlap the candidate mutation
    /*
    BOOST_FOREACH(AlignedVariant av, close_somatic) {
        if(leftPos > av.getStartHap())
            leftPos = av.getStartHap();
        if(rightPos < av.getStartHap())
            rightPos = av.getStartHap();
    }
    BOOST_FOREACH(AlignedVariant av, close_germline) {
        if(leftPos > av.getStartHap())
            leftPos = av.getStartHap();
        if(rightPos < av.getStartHap())
            rightPos = av.getStartHap();
    }
     */
    return AlignedCandidates(vi.chr, variants, leftPos, rightPos);
}

void HapMuC::mutationCall(const string & variantsFileName)
{
    LOG(logDEBUG) << "mutationCall " << variantsFileName << endl;
    ofstream output;
    string callsFile=params.fileName; callsFile.append(".calls.txt");
    OutputData oData=params.makeMutationData(output);
    output.open(callsFile.c_str());
    if (!output.is_open()) {
        throw(string("Cannot open file ").append(callsFile).append(" for writing."));
    }
    oData.outputLine(oData.headerString());
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
#ifdef LOGDEBUG
    system(("mkdir -p " + params.fileName + "logs").c_str());
#endif
    while(getline(varfile, line)) {
        vector<AlignedVariant> close_somatic_vars, close_germline_vars;
        AlignedCandidates candidateVariants = getCandidateVariants(line, close_somatic_vars, close_germline_vars);
        if (candidateVariants.variants.size()==0) continue;
        LOG(logDEBUG) << "for each candidate variants(size=" << candidateVariants.variants.size() << ")" << endl;
        vector<Read> normalReads;
        vector<Read> tumorReads;
        vector<Read> mergedReads;
        uint32_t pos, leftPos, rightPos;
        // get lowest and highest position
        leftPos = candidateVariants.leftPos - params.maxReadLength;
        rightPos = candidateVariants.rightPos + params.maxReadLength;

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
            LOG(logERROR) << "leftPos: " << leftPos << " oldLeftPosForN: " << oldLeftPosForN << endl;
            LOG(logERROR) << "Candidate variant files must be sorted on left position of window!" << endl;
            exit(1);
        }


        // TODO either add tid to AlignedVariant or infer it from the vector of aligned variants
        // change alige
        index++;
        bool skipped = false;

        LOG(logINFO) << "****" << endl;
        LOG(logINFO) << " tid: " << params.tid << " pos: " << pos << " leftPos: " << leftPos << " " << " rightPos: " << rightPos << endl;

        string message="ok";
        try {
            normalReads.clear();
            tumorReads.clear();
            getReadsFromBams(normalBams, leftPos, rightPos, normalReads, oldLeftPosForN, oldRightFetchReadPosForN, normalReadBuf, resetReadBuffer);
            getReadsFromBams(tumorBams, leftPos, rightPos, tumorReads, oldLeftPosForT, oldRightFetchReadPosForT, tumorReadBuf, resetReadBuffer);
            LOG(logDEBUG) << "read size: n, t =" << normalReads.size() << " " << tumorReads.size() << endl;

            uint32_t rs=(int(leftPos)>params.minReadOverlap)?(leftPos-params.minReadOverlap):0;
            uint32_t re=rightPos+params.minReadOverlap;
            string refSeq=getRefSeq(rs+1, re+1);
            string refSeqForAlign=getRefSeq(leftPos+1, rightPos+1);
            MutationCall::computeBayesFactors(normalReads, tumorReads, pos, leftPos, rightPos, candidateVariants, oData, params, refSeq, refSeqForAlign, close_somatic_vars, close_germline_vars, index);
            LOG(logDEBUG) << "after read size: n, t =" << normalReads.size() << " " << tumorReads.size() << endl;

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
            LOG(logERROR) << "skipped " << params.tid << " " << pos << " reason: " << message << endl;
            OutputData::Line line(oData);
            VariantInfo v = candidateVariants.variants[0].info;
            line.set("chr", params.tid);
            line.set("start", v.start);
            line.set("end", v.end);
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

            // reset read buffer: all reads will be fetched again
            resetReadBuffer = true;
        } else {
            resetReadBuffer = false;
        }

        oldLeftPosForT = leftPos;
        oldLeftPosForN = leftPos;
    }

    // clean up read buffer
    for (size_t r=0;r<normalReadBuf.size();r++) {
        if (normalReadBuf[r]!=NULL) delete normalReadBuf[r];
    }
    for (size_t r=0;r<tumorReadBuf.size();r++) {
        if (tumorReadBuf[r]!=NULL) delete tumorReadBuf[r];
    }
}

void getParameters(cmdline::parser & a, Parameters & params)
{
    params.maxReads=a.get<int>("maxReads");
    params.mapQualThreshold=a.get<double>("mapQualThreshold");
    params.minReadOverlap=a.get<int>("minReadOverlap");
    params.maxReadLength=a.get<int>("maxReadLength");

    params.priorSNP=a.get<double>("priorSNP");
    params.priorIndel=a.get<double>("priorIndel");

    params.refFileName=a.get<string>("ref");

    params.obsParams.pError=a.get<double>("pError");
    params.obsParams.pMut=a.get<double>("pMut");

    params.obsParams.maxLengthIndel=a.get<int>("maxLengthIndel");
    params.obsParams.maxLengthDel=params.obsParams.maxLengthIndel;
    params.obsParams.mapQualThreshold=a.get<double>("capMapQualThreshold");
    params.obsParams.padCover = a.get<int>("flankRefSeq");
    params.obsParams.maxMismatch = a.get<int>("flankMaxMismatch");

    params.showReads=a.exist("showReads")?true:false;
    params.quiet=a.exist("quiet")?true:false;

    params.showHapAlignments=a.exist("showHapAlignments")?true:false;
    params.filterReadAux=a.get<string>("filterReadAux");
}

int main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("ref", 'f', "fasta reference sequence (should be indexed with .fai file)", true, "");
    a.add<string>("tumor", 'a', "tumor bam file", true, "");
    a.add<string>("normal", 'b', "normal bam file", true, "");
    a.add<string>("out", 'o', "file-prefix for output results", true, "");
    a.add<string>("windows", 'w', "file with candidate variants to be tested.", true, "");
    
    a.add("quiet", '\0', "quiet output");
    
    a.add<int>("flankRefSeq", '\0', "#bases of reference sequence of indel region", false, 2, cmdline::range(1, 10));
    a.add<int>("flankMaxMismatch", '\0', "max number of mismatches in indel region", false, 2, cmdline::range(1, 10));
    
    a.add<double>("priorSNP", '\0', "prior probability of a SNP site", false, 1.0/1000.0);
    a.add<double>("priorIndel", '\0', "prior probability of a detected indel not being a sequencing error", false, 1.0/10000.0);
    a.add<int>("maxReads", '\0', "maximum number of reads in likelihood computation", false, 50000);
    a.add<double>("mapQualThreshold", '\0', "lower limit for read mapping quality", false, 0.9);
    a.add<double>("capMapQualThreshold", '\0', "upper limit for read mapping quality in observationmodel_old (phred units)", false, 100.0);
    a.add<int>("minReadOverlap", '\0', "minimum overlap between read and haplotype", false, 100);
    a.add<int>("maxReadLength", '\0', "maximum length of reads", false, 100);
    
    a.add<string>("filterReadAux", '\0', "match string for exclusion of reads based on auxilary information", false, "");
    
    a.add<double>("pError", '\0', "probability of a read indel", false, 5e-4);
    a.add<double>("pMut", '\0', "probability of a mutation in the read", false, 1e-5);
    a.add<int>("maxLengthIndel", '\0', "maximum length of a _sequencing error_ indel in read", false, 5);
    
    a.add("showHapAlignments", '\0', "show for each haplotype which reads map to it");
    a.add("showReads", '\0', "show reads");
    
    a.parse_check(argc, argv);
    try {
        // extract required parameters
        string file;
        int multipleFiles=0;

        string faFile=a.get<string>("ref");
        string outputFile=a.get<string>("out");

        string modelType="probabilistic";
        Parameters params(string("1"), outputFile, modelType);
        getParameters(a, params);

        FILE* pFile = fopen((outputFile + ".log").c_str(), "w");
        Output2FILE::Stream() = pFile;
        //  FILELog::ReportingLevel() = FILELog::FromString("DEBUG");
        string varFile = a.get<string>("windows");
        string tumorBF = a.get<string>("tumor");
        string normalBF = a.get<string>("normal");
        HapMuC hapmuc(normalBF, tumorBF, params);
        hapmuc.params.print();
        hapmuc.mutationCall(varFile);
    }
    catch (string s) {
        LOG(logERROR) << "Exception: " << s << endl;
        exit(1);
    }
    return 0;
}