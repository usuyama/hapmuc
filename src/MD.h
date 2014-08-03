#ifndef variant_parser_MD_h
#define variant_parser_MD_h

#include <cstdlib>
#include <iostream>
#include <string>
#include <cctype>

class MD {
    static int seek_non_digit_char(std::string &str, int current_idx) {
        int idx = current_idx;
        while (idx < str.length()) {
            if (std::isdigit(str[idx])) {
                idx++;
            } else {
                break;
            }
        }
        return idx;
    }
    
public:
    typedef enum { DEL, MISMATCH, MATCH } Type;
    Type type;
    int length;
    std::string ref;
    
    MD(Type _type, int _length, std::string _ref = "") {
        type = _type;
        length = _length;
        ref = _ref;
    }
    
    // String for mismatching positions.
    // Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    static std::vector<MD> from_md_string(std::string &md_str) {
        std::vector<MD> mds;
        int left = 1;
        // check first [0-9]+
        int right = seek_non_digit_char(md_str, left);
        mds.push_back(MD(MATCH, atoi(md_str.substr(left, right - left).c_str())));
        left = right;
        // next
        while (left < md_str.length()) {
            if (md_str[left] == '^') {
                // deletion
                right = left + 1;
                // find the last base for this deletion
                while (right < md_str.length()) {
                    if(std::isdigit(md_str[right])) {
                        break;
                    }
                    right++;
                }
                left++; // for removing ^
                mds.push_back(MD(DEL, right - left, md_str.substr(left, right - left)));
            } else {
                // mimatch
                mds.push_back(MD(MISMATCH, 1, std::string(1, md_str[left])));
                right = left + 1;
            }
            left = right;
            // check next [0-9]+
            right = seek_non_digit_char(md_str, left);
            mds.push_back(MD(MATCH, atoi(md_str.substr(left, right - left).c_str())));
            left = right;
        }
        return mds;
    }
    
    friend std::ostream& operator<<(std::ostream &stream, const MD &md) {
        std::string type_strings[] = { "DEL", "MISMATCH", "MATCH" };
        stream << type_strings[md.type] << " " << md.length << " " << md.ref;
        return stream;
    }
};
        
#endif
