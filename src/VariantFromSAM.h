#ifndef __VariantFromSAM__
#define __VariantFromSAM__

#include <iostream>
#include <vector>

#include "MD.h"
#include "CIGAR.h"

class VariantFromSAM {
public:
    typedef enum { MATCH, MISMATCH, DEL, INS, SOFTCLIP } Type;
    Type type;
    int left, right, length;
    std::string ref, obs;
    static std::vector<VariantFromSAM> from_cigar_and_md(const std::vector<CIGAR> &cigars,
                                                         const std::vector<MD> &mds,
                                                         int left_most_pos,
                                                         const std::string &seq);
    static VariantFromSAM from_cigar(const CIGAR &cigar);
    VariantFromSAM(Type _type, std::string _ref, std::string _obs, int _left, int _right);
    VariantFromSAM() {};
    friend bool operator==(const VariantFromSAM &v1, const VariantFromSAM &v2) {
        return v1.type == v2.type
               && v1.left == v2.left
               && v1.right == v2.right
               && ((v1.type == SOFTCLIP)
                   || (v1.ref == v2.ref && v1.obs == v2.obs));
    }
    
    friend std::ostream& operator<<(std::ostream &os, const VariantFromSAM &v) {
        std::string type_strings[] = { "MATCH", "MISMATCH", "DEL", "INS", "SOFTCLIP" };
        os << type_strings[v.type] << " " << v.left << "-"
           << v.right << " " << v.ref << " " << v.obs << " len: " << v.length;
        return os;
        }
};

#endif
