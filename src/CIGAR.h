#ifndef variant_parser_CIGAR_h
#define variant_parser_CIGAR_h

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

class CIGAR {
public:
    typedef enum { INS, DEL, MATCH, SOFTCLIP } Type;
    Type type;
    int length;
    
    CIGAR() {};
    
    CIGAR(CIGAR::Type _type, int _length) {
        type = _type;
        length = _length;
    }
    
    static std::vector<CIGAR> from_cigar_string(std::string cigar_str) {
        std::vector<CIGAR> cigars;
        std::string tmp_len = "";
        for (int i = 0;i < cigar_str.length();i++) {
            char &c = cigar_str[i];
            if (std::isdigit(static_cast<unsigned char>(c))) {
                tmp_len += c;
            } else {
                // make new cigar instance
                CIGAR cigar;
                switch (c)
                {
                    case 'M':
                        cigar.type = CIGAR::MATCH;
                        break;
                    case 'I':
                        cigar.type = CIGAR::INS;
                        break;
                    case 'D':
                        cigar.type = CIGAR::DEL;
                        break;
                    case 'S':
                        cigar.type = CIGAR::SOFTCLIP;
                        break;
                    case 'H':
                        throw std::string("hardclip_in_cigar");
                        break;
                    case 'N':
                        throw std::string("N_in_cigar");
                        break;
                    case 'P':
                        throw std::string("P_in_cigar");
                        break;
                    default:
                        throw std::string("unsupported_cigar");
                }
                cigar.length = atoi(tmp_len.c_str());
                cigars.push_back(cigar);
                tmp_len = "";
            }
        }
        return cigars;
    }
    
    friend std::ostream& operator<<(std::ostream &os, const CIGAR &cigar) {
        std::string type_strings[] = { "INS", "DEL", "MATCH", "SOFTCLIP" };
        os << type_strings[cigar.type] << " " << cigar.length;
        return os;
   }
        
};

#endif
