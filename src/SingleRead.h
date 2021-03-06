#ifndef __haplotype_builder__SingleRead__
#define __haplotype_builder__SingleRead__

#include <string>
#include <vector>

#include "IRead.h"
#include "SamRead.h"
#include "Variant.h"

class SingleRead : public IRead {
public:
    SamRead *read;

    virtual std::string getSeqName() const;
    virtual bool coverPosition(int pos) const;
    virtual bool hasVariant(const VariantFromSAM &variant) const;
    virtual double getMapQuality() const;
    virtual bool hasIndel() const;
    virtual bool hasSoftClip() const;
    virtual double getMapQualityPrimary() const;
    virtual bool hasIndelPrimary() const;
    virtual bool hasSoftClipPrimary() const;
    virtual SamRead *getPrimarySamRead() const;

    static std::vector<IRead *> toIReadVector(std::vector<SingleRead> &reads);
    SingleRead(SamRead *read);
};

#endif /* defined(__haplotype_builder__SingleRead__) */
