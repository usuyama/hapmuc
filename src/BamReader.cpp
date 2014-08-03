/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include "BamReader.h"

#include <list>
#include <vector>
#include <map>

#include "log.h"
#include "MyBam.hpp"
#include "PairedRead.h"
#include "SingleRead.h"

BamReader::BamReader(const MyBam *_bam, const faidx_t *_fai,
                     int _maxReads) {
    fai = _fai;
    bam = _bam;
    readBuffer = std::vector<SamRead *>();
    maxReads = _maxReads;
}

// get reads, which cover the given chr:position
// reset indicates whether we reset the readBuffer this time
std::vector<IRead *> BamReader::getPairedReads(std::string chr,
                                               int targetPosition,
                                               int left,
                                               int right,
                                               int maxInsertSize) {
    // modifing window
    const int PAD = 10;
    left -= PAD;
    left = (left < 0) ? 0 : left;
    right += PAD;

    // readBuffer is responsible for the memory management of reads.
    for (size_t r = 0; r < readBuffer.size(); r++) {
        if (readBuffer[r] != NULL) {
            delete readBuffer[r];
            readBuffer[r] = NULL;
        }
    }
    readBuffer.clear();
    
    LOG(logDEBUG) << "leftFetchReadPos: " << left << " rightFetchReadPos: " << right << std::endl;
    
    //
    // get reads using samtools/bam_fetch function
    //
    int numReads = (int) readBuffer.size();
    std::vector<SamRead *> newReads;
    LOG(logDEBUG) << "Fetching reads...." << std::endl;
    SamRead::FetchReadData data(&newReads, fai, bam, numReads, maxReads);
    try {
        bam_fetch(bam->bf, bam->idx, bam->getTID(chr),
                  left , right, &data, &SamRead::fetchFunc);
    } catch (std::string &s) {
        LOG(logDEBUG) << s << ": free read buffer" << std::endl;
        for (size_t r = 0; r < newReads.size(); r++) {
            delete newReads[r];
            newReads[r] = NULL;
        }
        newReads.clear();
        throw s;
    }
    numReads = data.numReads;
    LOG(logDEBUG) << numReads << " reads from bam_fetch" << std::endl;
    std::swap(readBuffer, newReads);
    
    // create hash of seq name and Read instance
    // for finding mate-pair of a given read later
    std::map<std::string, std::list<int> > seqNameHash; // query name to read idx
    std::map<std::string, std::list<int> >::const_iterator hashIt;
    for (size_t r = 0; r < readBuffer.size(); r++) {
        if (!readBuffer[r]->isUnmapped()) {
            seqNameHash[std::string(bam1_qname(readBuffer[r]->getBam()))].push_back((int)r);
            if (seqNameHash[std::string(bam1_qname(readBuffer[r]->getBam()))].size() > 2) {
                LOG(logERROR) << "there are more than two reads with same seq_name..." << std::endl;
                throw std::string("duplicate_seq_names");
            }
        }
    }
    
    std::map<std::string, bool> usedSeqNameHash;
    
    //
    // Filter reads by several criteria, and if it passed
    // we make a ReadPair by finding mate-pair using seqNameHash
    //
    std::vector<IRead *> readPairs;
    for (int r = 0; r < int(readBuffer.size()); r++) {
        SamRead &read = *readBuffer[r];
        std::string qname = std::string(bam1_qname(read.getBam()));
        // LOG(logDEBUG) << "qname: " << qname << " ";
        // LOGP(logDEBUG) << "reads[r].pos: " << read.pos
        // << " reads[r].getEndPos(): " << read.getEndPos() << std::endl;

        if (read.getBam()->core.mtid == 0 || read.getBam()->core.tid == 0) {
            LOG(logWARNING) << "TID null: reads[" << r << "]: "
            << bam1_qname(read.getBam()) << ", maybe you should use `samtools fix`" << std::endl;
        }

        if (read.isUnmapped()) {
            LOG(logDEBUG) << "filtered(unmapped): reads[" << r << "]: "
            << bam1_qname(read.getBam()) << std::endl;
            continue;
        } else if (read.mateIsUnmapped()) {
            LOG(logDEBUG) << "filtered(mate is not mapped): reads[" << r << "]: "
            << bam1_qname(read.getBam()) << std::endl;
            continue;
        } else if (read.getBam()->core.mtid != 0 &&
                   read.getBam()->core.tid != 0 &&
                   read.getBam()->core.mtid != read.getBam()->core.tid) {
            // chromosomes of this read and its pair are different
            LOG(logWARNING) << "TIDERR: reads[" << r << "]: "
            << bam1_qname(read.getBam()) << " matePos: "
            << read.matePos << " mateLen: " << read.mateLen
            << " (" << read.getBam()->core.mtid << " vs "
            << read.getBam()->core.tid << ")" << std::endl;
            continue;
        } else if (std::abs(read.pos - read.matePos) > maxInsertSize) {
            LOG(logWARNING) << "filtered(mate_read_too_far): reads[" << r << "]: "
            << bam1_qname(read.getBam()) << " pos: " << read.pos << " matePos: "
            << read.matePos << std::endl;
            continue;
        } else if (usedSeqNameHash.find(qname) != usedSeqNameHash.end()) {
            // LOG(logDEBUG) << "we already used this read as mate pair."  << std::endl;
            continue;
        } else {
            hashIt = seqNameHash.find(std::string(bam1_qname(read.getBam())));
            if (hashIt->second.size() < 2) {
                // mate-read not found
                readPairs.push_back(new PairedRead(&read, NULL));
            } else {
                for (std::list<int>::const_iterator it = hashIt->second.begin(); it != hashIt->second.end(); it++) {
                    // skip for the read itself
                    if (*it != r) {
                        SamRead &mateRead = *readBuffer[*it];
                        read.mateLen = (int) mateRead.seq.length();
                        read.matePos = mateRead.pos;
                        if (read.matePos != read.getBAMMatePos()) {
                            LOG(logWARNING) << "matepos inconsistency!" << std::endl;
                            LOG(logWARNING) << read.matePos
                            << " " << read.getBAMMatePos() << std::endl;
                        } else {
                            if (read.cover(targetPosition)) {
                                readPairs.push_back(new PairedRead(&read, &mateRead));
                            } else {
                                readPairs.push_back(new PairedRead(&mateRead, &read));
                            }
                        }
                    }
                }
            }
            usedSeqNameHash[qname] = true;
        }
    }
    
    LOG(logINFO) << "Number of read pairs: " << readPairs.size() << std::endl;
    if (readPairs.size() < 5) {
        LOG(logERROR) << "size=" << readPairs.size() << std::endl;
        throw std::string("too_few_reads");
    } else if (readPairs.size() >= maxReads) {
        LOG(logERROR) << "size=" << readPairs.size() << std::endl;
        throw std::string("too_many_reads_in_window");
    }
    
    return readPairs;
}



// get reads, which cover the given chr:position
// reset indicates whether we reset the readBuffer this time
std::vector<IRead *> BamReader::getSingleReads(std::string chr,
                                               int targetPosition,
                                               int left,
                                               int right) {
    // modifing window
    const int PAD = 10;
    left -= PAD;
    left = (left < 0) ? 0 : left;
    right += PAD;

    // readBuffer is responsible for the memory management of reads.
    for (size_t r = 0; r < readBuffer.size(); r++) {
        if (readBuffer[r] != NULL) {
            delete readBuffer[r];
            readBuffer[r] = NULL;
        }
    }
    readBuffer.clear();

    LOG(logDEBUG) << "leftFetchReadPos: " << left << " rightFetchReadPos: " << right << std::endl;

    //
    // get reads using samtools/bam_fetch function
    //
    int numReads = (int) readBuffer.size();
    std::vector<SamRead *> newReads;
    LOG(logDEBUG) << "Fetching reads...." << std::endl;
    SamRead::FetchReadData data(&newReads, fai, bam, numReads, maxReads);
    try {
        bam_fetch(bam->bf, bam->idx, bam->getTID(chr),
                  left , right, &data, &SamRead::fetchFunc);
    } catch (std::string &s) {
        LOG(logDEBUG) << s << ": free read buffer" << std::endl;
        for (size_t r = 0; r < newReads.size(); r++) {
            delete newReads[r];
            newReads[r] = NULL;
        }
        newReads.clear();
        throw s;
    }

    numReads = data.numReads;
    LOG(logDEBUG) << numReads << " reads from bam_fetch" << std::endl;
    std::swap(readBuffer, newReads);

    //
    // Filter reads by several criteria, and if it passed
    //
    std::vector<IRead *> singleReads;
    for (int r = 0; r < int(readBuffer.size()); r++) {
        SamRead &read = *readBuffer[r];
        
        // LOG(logDEBUG) << "qname: " << read.seq_name << " ";
        // LOGP(logDEBUG) << "reads[r].pos: " << read.pos
        // << " reads[r].getEndPos(): " << read.getEndPos() <<
        // " mapqual: " << read.mapQual << std::endl;

        if (read.isUnmapped()) {
            LOG(logDEBUG) << "filtered(unmapped): reads[" << r << "]: "
            << bam1_qname(read.getBam()) << std::endl;
            continue;
        } else {
            singleReads.push_back(new SingleRead(&read));
        }
    }

    LOG(logINFO) << "Number of single-reads: " << singleReads.size() << std::endl;
    if (singleReads.size() < 5) {
        LOG(logERROR) << "size=" << singleReads.size() << std::endl;
        throw std::string("too_few_reads");
    } else if (singleReads.size() >= maxReads) {
        LOG(logERROR) << "size=" << singleReads.size() << std::endl;
        throw std::string("too_many_reads_in_window");
    }

    return singleReads;
}

BamReader::~BamReader() {
    for (int i = 0; i < readBuffer.size(); i++) {
        delete readBuffer[i];
        readBuffer[i] = NULL;
    }
}

