#include "SingleRead.h"

SingleRead::SingleRead(SamRead *_read) {
    read = _read;
}

std::string SingleRead
::getSeqName() const {
    return read->seq_name;
}

std::vector<IRead *> SingleRead::toIReadVector(std::vector<SingleRead> &reads) {
    std::vector<IRead *> ireads;
    ireads.reserve(reads.size());
    for (std::vector<SingleRead>::iterator it = reads.begin(); it != reads.end(); it++) {
        ireads.push_back(&(*it));
    }
    return ireads;
}

bool SingleRead::coverPosition(int pos) const {
    return read->cover(pos);
}

bool SingleRead::hasVariant(const VariantFromSAM &variant) const {
    return read->hasVariantFromSam(variant);
}

double SingleRead::getMapQuality() const {
    return read->mapQual;
}

bool SingleRead::hasIndel() const {
    return read->hasIndel;
}

bool SingleRead::hasSoftClip() const {
    return read->hasSoftClip;
}

double SingleRead::getMapQualityPrimary() const {
    return read->mapQual;
}

bool SingleRead::hasIndelPrimary() const {
    return read->hasIndel;
}

bool SingleRead::hasSoftClipPrimary() const {
    return read->hasSoftClip;
}

SamRead *SingleRead::getPrimarySamRead() const {
    return read;
}
