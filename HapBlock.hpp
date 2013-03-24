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
#ifndef HAPBLOCK_HPP_
#define HAPBLOCK_HPP_
#include <stdint.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include "Haplotype.hpp"
using namespace std;

class HapBlock
{
public:
	map<Haplotype, int> haplotypes;
	bool operator<(const HapBlock & hb) const { return pos0<hb.pos0; };


	HapBlock(const Haplotype & seq, uint32_t start);
	HapBlock(const HapBlock & hb, uint32_t _start, uint32_t _len);
	bool hasHaplotype(const Haplotype & seq, uint32_t seqStart);
	uint32_t start() const { return pos0; };
	uint32_t end() const { return pos1; };
	uint32_t length() const { return end()-start()+1; };
	size_t size() const { return haplotypes.size(); };
	void insert(const Haplotype & seq);// { haplotypes[seq]++; }
	vector<pair<Haplotype,int> > getHaplotypes();
	void setFrequencies();
	friend ostream &operator<<(ostream &stream, const HapBlock &hb);
	static void showVector(ostream &stream,const vector<HapBlock*> & hapBlocks, uint32_t midPos);
	void setType(int _type) { type=_type; };
	int getType() const { return type; };

	static const int NORMAL=0;
	static const int INSERT=1;
private:
	uint32_t pos0, pos1;
	int type;
};

#endif /*HAPBLOCK_HPP_*/
