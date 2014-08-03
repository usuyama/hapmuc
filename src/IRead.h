#ifndef __IRead_h__
#define __IRead_h__

#include <string>

#include "VariantFromSAM.h"
#include "Variant.h"
#include "SamRead.h"

// interface
class IRead {
public:
    virtual ~IRead() {}
    virtual std::string getSeqName() const = 0;
    virtual bool coverPosition(int pos) const = 0;
    virtual bool hasVariant(const VariantFromSAM &variant) const = 0;
    virtual double getMapQuality() const = 0;
    virtual bool hasSoftClip() const = 0;
    virtual bool hasIndel() const = 0;
    virtual double getMapQualityPrimary() const = 0;
    virtual bool hasSoftClipPrimary() const = 0;
    virtual bool hasIndelPrimary() const = 0;
    virtual SamRead *getPrimarySamRead() const = 0;
};

#endif
