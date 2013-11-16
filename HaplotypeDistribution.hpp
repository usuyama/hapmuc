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
#ifndef HAPLOTYPEDISTRIBUTION_HPP_
#define HAPLOTYPEDISTRIBUTION_HPP_
#include <string>
#include <assert.h>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <map>
#include <cmath>
#include "bam.h"
#include "Haplotype.hpp"
#include "HapBlock.hpp"
#include "foreach.hpp"
#include "VariantFile.hpp"
using namespace std;

class HaplotypeDistribution
{
// methods
public:
	HaplotypeDistribution(uint32_t _midPos, const string & refSeq, uint32_t refSeqStart);

	void insertRead(const bam1_t *b);
	static int fetchFuncInsertRead(const bam1_t *b, void *data);
	pair<vector<Haplotype>, vector<double> > enumerateHaplotypes(double th);
	friend ostream &operator<<(ostream &stream, const HaplotypeDistribution &hb);
	void insertSeq(Haplotype & seq, uint32_t seqStart);
	void check();
	size_t getNumberOfHaplotypes(uint32_t start, uint32_t end) const;
	size_t getNumberOfHaplotypes(uint32_t start, uint32_t end, double minFreq) const;
	void setFrequencies();
	~HaplotypeDistribution();
	set<Haplotype> getIndelsAtMidPos() const { return indelsAtMidPos; };
	vector<Variant> getIndelVariantsAtMidPos();
protected:
	void updateBlock(HapBlock *hb, const Haplotype & seq, uint32_t seqStart);
	void newBlock(HapBlock *hb);
	void deleteBlock(int idx);
	void splitBlock(int idx, const Haplotype & seq, uint32_t seqStart);
	int getFirstOverlappingBlock(uint32_t seqStart, uint32_t seqEnd) const;

	uint32_t len;
	uint32_t pos0, pos1, midPos;
	vector<HapBlock*> hapBlocks;
	map<int, HapBlock*> insertions;

	set<Haplotype> indelsAtMidPos;

};


class HDHapBlock
{
public:
	HDHapBlock() { };
	vector<Haplotype> haps;
	uint32_t start, end;
	int type;
};

#endif
