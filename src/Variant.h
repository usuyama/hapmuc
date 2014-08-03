/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef __Variant__
#define __Variant__

#include <string>
#include <iostream>

#include "VariantFromSAM.h"

class Variant {
public:
    int length, startInGenome, endInGenome;
    std::string seq, originalString;
    std::string ref, obs;
    typedef enum { INS, DEL, SNP } Type;
    Type type;
    bool isIndel() const;
    bool isSNP() const;
    Variant();
    Variant(const std::string &_str, int _startInGenome);
    friend std::ostream &operator<<(std::ostream &_out, const Variant &variant) {
        _out << "[" << variant.originalString << " "
        << variant.startInGenome << " "
        << variant.endInGenome << "]";
        return _out;
    }
    bool isOverlap(const Variant &var, int margin = 0) const;
    std::string getSymbol() const;
private:
    void initFromString(const std::string &str);
};


#endif /* defined(__Variant__) */
