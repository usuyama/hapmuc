/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include "Variant.h"

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

void Variant::initFromString(const std::string &str) {
    int ok = 1;
    if (str.size() > 1) {
        if (str[0] == '-') {
            length = int(str.size()) - 1;
            seq = str.substr(1, length);
            ref = seq;
            obs = "-";
            type = Variant::DEL;
        } else if (str[0] == '+') {
            length = (int(str.size()) - 1);
            seq = str.substr(1, length);
            ref = "-";
            obs = seq;
            type = Variant::INS;
        } else if (str.size() == 4 && str[1] == '=' && str[2] == '>') {
            type = Variant::SNP;
            seq = str;
            ref = str.substr(0, 1);
            obs = str.substr(3, 1);
            length = 1;
        } else {
            ok = 0;
        }
    } else {
        ok = 0;
    }
    if (!ok) {
        std::cout << "input std::string: " << str << std::endl;
        throw std::string("Unrecognized variant");
    }
    this->originalString = str;
}

std::string Variant::getSymbol() const {
    std::stringstream ss;
    ss << startInGenome << ":";
    if (type == INS) {
        ss << "+" << obs;
    } else if (type == DEL) {
        ss << "-" << ref;
    } else {
        ss << ref << "=>" << obs;
    }
    return ss.str();
}

bool Variant::isIndel() const {
    if (type == INS || type == DEL) {
        return true;
    } else {
        return false;
    }
}

bool Variant::isSNP() const {
    return type == Variant::SNP;
}

Variant::Variant(const std::string &_str, int _startInGenome) {
    initFromString(_str);
    startInGenome = _startInGenome;
    if (type != INS) {
        endInGenome = _startInGenome + length;
    } else {
        endInGenome = _startInGenome;
    }
}

Variant::Variant() {
    startInGenome = -1;
    endInGenome = -1;
    originalString = "NA";
    seq = "NA";
}

// distance:
// if 0:
// ---A--
// --A---
// => false
// if 1: true
bool Variant::isOverlap(const Variant &var, int margin) const {
    int s1 = startInGenome, e1 = endInGenome;
    int s2 = var.startInGenome, e2 = var.endInGenome;
    if (s1 > s2) {
        std::swap(s1, s2);
        std::swap(e1, e2);
    }
    // s1<--->e1
    //    s2<--->e2
    return e1 > (s2 - margin);
}