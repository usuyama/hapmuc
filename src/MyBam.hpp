/*
 * Dindel
 * http://www.sanger.ac.uk/resources/software/dindel/
 *
 * Copyright 2010, Kees Albers
 * Released under the GPL Version 3.
 */
/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * MyBam.hpp
 *
 *  Created on: Aug 27, 2009
 *      Author: caa
 */

#ifndef MYBAM_HPP
#define MYBAM_HPP

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include "faidx.h"
#include "bam.h"

class MyBam {
public:
    MyBam() {
        initialized = false;
    };
    MyBam(const std::string &bamFile) {
        // load bam file
        initialized = false;
        init(bamFile);
    }

    MyBam(const MyBam &myBam) {
        initialized = false;
        init(myBam.fileName);
    }
    int getTID(const std::string &str) const {
        std::map<std::string, int>::const_iterator it = strToTID.find(str);
        if (it == strToTID.end()) {
            throw std::string("Cannot find ID!");
        } else {
            return it->second;
        }
    }

    void destroy() {
        if (initialized) {
            bam_close(bf);
            bam_header_destroy(bh);
            bam_index_destroy(idx);
        }
        initialized = false;
    }

    ~MyBam() {
        destroy();
    }

    bamFile bf;
    bam_header_t *bh;
    bam_index_t *idx;
    std::string fileName;
private:
    void init(const std::string &_fileName) {
        destroy();
        fileName = _fileName;
        bf = bam_open(fileName.c_str(), "r");
        if (!bf) {
            throw std::string("Cannot open BAM file.");
        }
        bh = bam_header_read(bf);
        for (int nt = 0; nt < bh->n_targets; nt++) {
            strToTID[std::string(bh->target_name[nt])] = nt;
        }

        idx = bam_index_load(fileName.c_str());
        initialized = true;
    }
    bool initialized;
    std::map<std::string, int> strToTID;
};

#endif
