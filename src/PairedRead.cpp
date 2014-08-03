#include "PairedRead.h"

PairedRead::PairedRead(SamRead *read, SamRead *mateRead) {
    first = read;
    second = mateRead;
}

std::string PairedRead::getSeqName() const {
    return first->seq_name;
}

std::vector<IRead *> PairedRead::toIReadVector(std::vector<PairedRead> &reads) {
    std::vector<IRead *> ireads;
    ireads.reserve(reads.size());
    for (std::vector<PairedRead>::iterator it = reads.begin(); it != reads.end(); it++) {
        ireads.push_back(&(*it));
    }
    return ireads;
}

bool PairedRead::coverPosition(int pos) const {
    if (first->cover(pos)) return true;
    if (second != NULL && second->cover(pos)) return true;
    return false;
}

bool PairedRead::hasVariant(const VariantFromSAM &variant) const {
    if (first->hasVariantFromSam(variant)) return true;
    if (second != NULL && second->hasVariantFromSam(variant)) return true;
    return false;
}

double PairedRead::getMapQuality() const {
    if (second != NULL && !second->isUnmapped()) {
        return first->mapQual * second->mapQual;
    } else {
        return first->mapQual;
    }
}

bool PairedRead::hasIndel() const {
    if (second != NULL && !second->isUnmapped()) {
        return first->hasIndel || second->hasIndel;
    } else {
        return first->hasIndel;
    }
}

bool PairedRead::hasSoftClip() const {
    if (second != NULL && !second->isUnmapped()) {
        return first->hasSoftClip || second->hasSoftClip;
    } else {
        return first->hasSoftClip;
    }
}

bool PairedRead::hasIndelPrimary() const {
    return first->hasIndel;
}

bool PairedRead::hasSoftClipPrimary() const {
    return first->hasSoftClip;
}

double PairedRead::getMapQualityPrimary() const {
    return first->mapQual;
}

SamRead *PairedRead::getPrimarySamRead() const {
    return first;
}