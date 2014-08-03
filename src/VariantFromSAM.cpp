#include "VariantFromSAM.h"

#include <iostream>
#include <vector>

#include "MD.h"
#include "CIGAR.h"

VariantFromSAM::VariantFromSAM(Type _type, std::string _ref, std::string _obs, int _left, int _right) {
    type = _type;
    ref = _ref;
    obs = _obs;
    left = _left;
    right = _right;
    length = right - left;
}

VariantFromSAM VariantFromSAM::from_cigar(const CIGAR &cigar) {
    VariantFromSAM var;
    var.length = cigar.length;
    if (cigar.type == CIGAR::MATCH) var.type = MATCH;
    if (cigar.type == CIGAR::DEL) {
        var.type = DEL;
        var.obs = "-";
    }
    if (cigar.type == CIGAR::INS) {
        var.type = INS;
        var.ref = "-";
    }
    if (cigar.type == CIGAR::SOFTCLIP) {
        var.type = SOFTCLIP;
        var.ref = "?";
    }
    return var;
}

template<typename T>
inline void set_next_item(T &current_item, int &current_idx, const std::vector<T> &items) {
    current_idx++;
    if (current_idx < items.size()) {
        current_item = items[current_idx];
    }
}

std::vector<VariantFromSAM> VariantFromSAM::from_cigar_and_md(const std::vector<CIGAR> &cigars,
                                                              const std::vector<MD> &mds,
                                                              int left_most_pos,
                                                              const std::string &seq) {
    std::vector<VariantFromSAM> variants;
    int cigar_idx = 0, md_idx = 0;
    CIGAR current_cigar = cigars[0];
    MD current_md = mds[0];
    while (cigar_idx < cigars.size() && md_idx < mds.size()) {
        if (current_cigar.type == CIGAR::SOFTCLIP) {
            variants.push_back(VariantFromSAM::from_cigar(current_cigar));
            set_next_item(current_cigar, cigar_idx, cigars);
        } else if (current_md.type == MD::DEL) {
            VariantFromSAM v = VariantFromSAM::from_cigar(current_cigar);
            v.ref = current_md.ref;
            variants.push_back(v);
            set_next_item(current_cigar, cigar_idx, cigars);
            set_next_item(current_md, md_idx, mds);
        } else if (current_md.type == MD::MISMATCH) {
            VariantFromSAM v;
            v.type = MISMATCH;
            v.length = 1;
            v.ref = current_md.ref;
            variants.push_back(v);
            current_cigar.length -= 1;
            set_next_item(current_md, md_idx, mds);
        } else if (current_cigar.type == CIGAR::INS && current_md.type == MD::MATCH) {
            variants.push_back(VariantFromSAM::from_cigar(current_cigar));
            set_next_item(current_cigar, cigar_idx, cigars);
        } else if (current_cigar.type == CIGAR::MATCH && current_md.type == MD::MATCH) {
            if (current_cigar.length > current_md.length) {
                VariantFromSAM v;
                v.type = MATCH;
                v.length = current_md.length;
                variants.push_back(v);
                current_cigar.length -= current_md.length;
                set_next_item(current_md, md_idx, mds);
                if (current_cigar.length == 0) {
                    set_next_item(current_cigar, cigar_idx, cigars);
                }
            } else {
                variants.push_back(VariantFromSAM::from_cigar(current_cigar));
                current_md.length -= current_cigar.length;
                set_next_item(current_cigar, cigar_idx, cigars);
                if (current_md.length == 0) {
                    set_next_item(current_md, md_idx, mds);
                }
            }
        } else {
            throw std::string("something_strange_during_VariantFromSAM::from_cigar_and_md");
        }
    }
    
    //
    // setting left & right pos for variants
    //
    int k = left_most_pos;
    for (int i = 0;i < variants.size();i++) {
        VariantFromSAM &v = variants[i];
        v.left = k;
        if (v.type == INS) {
            v.right = v.left;
        } else {
            k += v.length;
            v.right = v.left + v.length;
        }
    }
    
    //
    // setting obs for variants
    //
    int b = 0;
    for (int i = 0;i < variants.size();i++) {
        VariantFromSAM &v = variants[i];
        if (v.type != DEL) {
            v.obs = seq.substr(b, v.length);
            b += v.length;
        }
    }

    return variants;
}