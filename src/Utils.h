#ifndef ReadTest_Utils_h
#define ReadTest_Utils_h

#include <vector>
#include <string>
#include <cmath>

class Utils {
public:
    static inline std::vector<std::string> split(std::string s, std::string delimiter) {
        std::vector<std::string> vec;
        size_t pos = 0;
        std::string token;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            vec.push_back(token);
            s.erase(0, pos + delimiter.length());
        }
        if (s.size() != 0) {
            vec.push_back(s);
        }
        return vec;
    }
    
    static inline bool nearlyEqual(double a, double b) {
        return std::abs(a - b) < 0.01;
    }
    
#include "float.h"
#include "math.h"
    // naively counts the number of supporting reads
    static inline std::vector<double> votes(const std::vector<std::vector<double> > &liks) {
        if (liks.size() == 0) {
            throw std::string("liks.size() == 0");
        }
        std::vector<double> vs = std::vector<double>(liks[0].size());
        for (int i = 0; i < liks.size(); i++) {
            double norm = 0.0;
            for (int j = 0; j < liks[i].size(); j++) {
                norm += exp(liks[i][j]);
            }
            for (int j = 0; j < liks[i].size(); j++) {
                vs[j] += exp(liks[i][j] - log(norm));
            }
        }
        return vs;
    }
};

#endif
