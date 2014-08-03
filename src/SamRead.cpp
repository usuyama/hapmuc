#include "SamRead.h"

#include <algorithm>
#include <cmath>

#include "bam.h"

#include "log.h"
#include "MD.h"
#include "CIGAR.h"
#include "VariantFromSAM.h"

int SamRead::fetchFunc(const bam1_t *b, void *data) {
    FetchReadData *ptr = (FetchReadData *) data;
    if (!((b->core.flag & BAM_FDUP) || (b->core.flag & BAM_FQCFAIL))) {
        try {
            ptr->reads->push_back(new SamRead(b, ptr->myBam->bh, ptr->fai));
            ptr->numReads++;
        } catch (std::string s) {
            ptr->numUnusualReads++;
            LOG(logERROR) << "fetchFuncVectorPooled[" << (std::string)reinterpret_cast<char *>(bam1_qname(b))
            << "]: " << s << std::endl;
        }
    }
    if (ptr->numReads > ptr->maxNumReads) {
        LOG(logWARNING) << "too_many_reads_in_window: " << ptr->numReads << std::endl;
        throw std::string("too_many_reads_in_window");
    }
    if (ptr->numReads % 10000 == 9999) {
        std::cout << "numreads: " << ptr->numReads << std::endl;
    }
    return 0;
}

bam1_t *SamRead::getBam() const {
    return bam;
}

SamRead::SamRead(const bam1_t *b, bam_header_t *_bamHeader, const faidx_t *_fai) {
    fai = _fai;
    const bam1_core_t *c = &b->core;
    uint32_t len = c->l_qseq;

    double mapPhred = (double) c->qual;
    mapQual = (1.0 - pow(10.0, -mapPhred / 10.0));

    if (mapQual < 0.0 || mapQual > 1.0 || std::isnan(mapQual) || std::isinf(mapQual)) {
        throw std::string("Phred error.");
    } else if (mapQual < 1e-16) {
        mapQual = 1e-16;
    } else if (mapQual > 1 - 1e-16) {
        mapQual = 1 - 1e-16;
    }

    pos = c->pos; // zero-based
    seq_name = (std::string)reinterpret_cast<char *>(bam1_qname(b));
    seq.reserve(len);
    qual.reserve(len);
    for (size_t x = 0; x < len; x++) {
        seq += (bam_nt16_rev_table[ bam1_seqi(bam1_seq(b), x) ]);
        // convert phred to probability
        double basePhred = (double)(((uint8_t *) bam1_qual(b))[x]);
        double q = (1.0 - pow(10.0, -basePhred / 10.0));
        if (q < 0.0 || q > 1.0 || std::isnan(q) || std::isinf(q)) {
            throw std::string("Phred error.");
        }
        if (q < 1e-16) {
            q = 1e-16;
        }
        if (q > 1.0 - 1e-16) {
            q = 1.0 - 1e-16;
        }
        qual.push_back(q);   // base quality is on log10 scale
    }

    bam = new bam1_t;
    *bam = *b;
    bam->data = new uint8_t[b->m_data];
    bam->m_data = b->m_data;
    for (int m = 0; m < b->m_data; m++) {
        bam->data[m] = b->data[m];
    }

    if (bam->core.flag & BAM_FREVERSE) {
        onReverseStrand = true;
    } else {
        onReverseStrand = false;
    }

    matePos = bam->core.mpos;
    mateLen = -1;
    this->bamHeader = _bamHeader;

    uint32_t *rawCigar = bam1_cigar(b);
    
    leftMostPos = pos;
    rightMostPos = getEndPos();

    if(!isUnmapped()) {
        rightMostPos = pos;
        hasIndel = false;
        hasHardClip = false;
        hasSoftClip = false;
        hasOtherCigarFlag = false;
        bool isLeftMost = true;
        softClipSize = 0;
        
        // set CIGAR string
        for (int k = 0; k < c->n_cigar; ++k) {
            int op = rawCigar[k] & BAM_CIGAR_MASK;
            int32_t len = rawCigar[k] >> BAM_CIGAR_SHIFT;
            // update position for the next cigar
            if (op == BAM_CMATCH) {
                cigars.push_back(CIGAR(CIGAR::MATCH, len));
                rightMostPos += len;
                isLeftMost = false;
            } else if (op == BAM_CINS) {
                cigars.push_back(CIGAR(CIGAR::INS, len));
                hasIndel = true;
                isLeftMost = false;
            } else if (op == BAM_CDEL) {
                cigars.push_back(CIGAR(CIGAR::DEL, len));
                rightMostPos += len;
                hasIndel = true;
                isLeftMost = false;
            } else if (op == BAM_CSOFT_CLIP) {
                rightMostPos += len;

                if (softClipSize < (int)len) {
                    softClipSize = (int)len;
                }
                hasSoftClip = true;
                cigars.push_back(CIGAR(CIGAR::SOFTCLIP, len));
                if (isLeftMost) {
                    leftMostPos -= len;
                    isLeftMost = false;
                }
            } else if (op == BAM_CHARD_CLIP) {
                hasHardClip = true;
            } else {
                hasOtherCigarFlag = true;
            }
        }

        parseVariants();
    }
}

void SamRead::getReadWithoutSoftClipping(std::string *_seq, std::vector<double> *_bq, int *_pos) const {
    *_pos = pos;
    if (isUnmapped() || !hasSoftClip) {
        *_seq = seq;
        *_bq = qual;
    } else {
        *_seq = "";
        _bq->clear();
        int index = 0;
        for (int k = 0;k < cigars.size();k++) {
            const CIGAR &c = cigars[k];
            if (c.type == CIGAR::MATCH || c.type == CIGAR::INS) {
                *_seq += seq.substr(index, c.length);
                _bq->insert(_bq->end(), qual.begin() + index, qual.begin() + index + c.length);
            } else if (c.type == CIGAR::SOFTCLIP) {
                // do nothing
            }
            if (c.type != CIGAR::DEL) {
                index += c.length;
            }
        }
    }
}

std::vector<VariantFromSAM> SamRead::getVariantsFromSAM() const {
    if (_isVariantsFromSAMAvailable) {
        return variantsFromSAM;
    } else {
        throw std::string("variant_from_sam_not_available");
    }
}

bool SamRead::hasVariantFromSam(const VariantFromSAM &variant) const {
    if (!_isVariantsFromSAMAvailable) {
        LOG(logWARNING) << "variant_from_sam_not_available" << std::endl;
        return false;
    } else {
        return std::find(variantsFromSAM.begin(), variantsFromSAM.end(), variant) != variantsFromSAM.end();
    }
}

bool SamRead::cover(int _pos) const {
    return leftMostPos <= _pos && rightMostPos >= _pos;
//    return pos <= _pos && getEndPos() >= _pos;
}

SamRead::~SamRead() {
    if (bam != NULL) {
        delete[] bam->data;
        delete bam;
        bam = NULL;
    }
}

void SamRead::parseVariants() {
    // set MD tag
    uint8_t *mdp = bam_aux_get(bam, "MD");
    if (mdp == NULL) {
        LOG(logWARNING) << "MD tag not found: " << seq_name << std::endl;
        LOG(logWARNING) << "please run samtools calmd" << std::endl;
        /*
// extern void bam_fillmd1_core(bam1_t *b, char *ref, int flag, int max_nm); // bam_md.c

        std::stringstream ss;
        ss << bam->core.tid << ":" << pos << "-" << getEndPos();
        char *ref;
        int len;
        ref = fai_fetch(fai, ss.str().c_str(), &len);
        if (len == 0) {
            throw std::string("faidx error: len==0");
        }
        bam_fillmd1_core(bam, ref, 0, 100);
        free(ref);
        mdp = bam_aux_get(bam, "MD");
         */
        _isVariantsFromSAMAvailable = false;
        throw std::string("md_tag_not_found");
    }

    std::string md_str((char *)mdp);
    LOG(logDEBUG2) << md_str << std::endl;
    mds = MD::from_md_string(md_str);

    // setting variants from CIGAR and MD
    if (!hasHardClip && !hasOtherCigarFlag) {
        variantsFromSAM = VariantFromSAM::from_cigar_and_md(cigars, mds, leftMostPos, seq);
        LOG(logDEBUG2) << seq_name << std::endl;
        LOG(logDEBUG2) << pos << std::endl;
        for (int i = 0;i < variantsFromSAM.size();i++) {
            LOG(logDEBUG2) << variantsFromSAM[i] << std::endl;
        }
        _isVariantsFromSAMAvailable = true;
    } else {
        _isVariantsFromSAMAvailable = false;
    }
}

bool SamRead::mateIsUnmapped() const {
    return (bam->core.flag & BAM_FMUNMAP) != 0;
}

bool SamRead::isUnmapped() const {
    return (bam->core.flag & BAM_FUNMAP) != 0;
}

int32_t SamRead::getBAMMatePos() const {
    return bam->core.mpos;
}

uint32_t SamRead::getEndPos() const {
    return bam->core.n_cigar ? bam_calend(&bam->core, bam1_cigar(bam)) : bam->core.pos + 1;
}
