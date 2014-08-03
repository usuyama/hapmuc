#ifndef __haplotype_builder__SamRead__
#define __haplotype_builder__SamRead__

#include <vector>

#include "MyBam.hpp"
#include "CIGAR.h"
#include "MD.h"
#include "VariantFromSAM.h"
#include "faidx.h"

class SamRead {
    std::vector<VariantFromSAM> variantsFromSAM;
    bool _isVariantsFromSAMAvailable;
    bool _variantsFromSAMCalcDone;
    void parseVariants();
public:
    // for bam_fetch
    class FetchReadData {
    public:
        FetchReadData(std::vector<SamRead *> *_reads, const faidx_t *_fai, const MyBam *_myBam, int _numReads = 0, int _maxNumReads = 100000) {
            reads = _reads;
            numReads = _numReads;
            myBam = _myBam;
            maxNumReads = _maxNumReads;
            numUnusualReads = 0;
            fai = _fai;
        }
        const faidx_t *fai;
        std::vector<SamRead *> *reads;
        const MyBam *myBam;
        int numReads;
        int maxNumReads;
        int numUnusualReads;
    };
    static int fetchFunc(const bam1_t *b, void *data);
    
    std::vector<CIGAR> cigars;
    std::vector<MD> mds;
    const faidx_t *fai;

    bool hasVariantFromSam(const VariantFromSAM &variant) const;
    std::vector<VariantFromSAM> getVariantsFromSAM() const;
    
    void getReadWithoutSoftClipping(std::string *seq, std::vector<double> *bq, int *pos) const;

    std::string seq_name, seq, cigar_str;
    std::vector<double> qual;
    double mapQual;
    int pos, matePos, mateLen;
    int leftMostPos, rightMostPos; // including soft-clipping bases
    int softClipSize;
    bam_header_t *bamHeader;
    bam1_t *bam;
    bam1_t *getBam() const;

    bool mateIsUnmapped() const;
    bool onReverseStrand;
    bool hasIndel;
    bool hasSoftClip, hasHardClip, hasOtherCigarFlag;
    bool isUnmapped() const;
    int32_t getBAMMatePos() const;
    uint32_t getEndPos() const;
    
    bool cover(int pos) const;
    
    SamRead(const bam1_t *b, bam_header_t *_bamHeader, const faidx_t *fai);
    ~SamRead();

    friend std::ostream& operator<<(std::ostream &stream, const SamRead &r) {
        stream << "SamRead " << r.pos << " " << r.leftMostPos << " " << r.onReverseStrand << " ";
        stream << r.mapQual << " " << r.seq.size();
        stream << r.seq << std::endl;
        for (int i = 0;i < r.qual.size();i++) {
            stream << r.qual[i] << " ";
        }
        return stream;
    }
private:
    SamRead(const SamRead&);
    SamRead operator=(const SamRead&);
};

#endif /* defined(__haplotype_builder__SamRead__) */
