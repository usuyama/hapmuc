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
#include "HapBlock.hpp"
#include <sstream>
#include <iostream>
using namespace std;

HapBlock::HapBlock(const HapBlock & hb, uint32_t _start, uint32_t _len)
{
	assert(hb.end()>=_start+_len-1);

	pos0=_start;
	pos1=_start+_len-1;
	if (pos1<pos0) {
		cout << "SMALLER" << endl;
	}
	type = HapBlock::NORMAL;
	assert(pos1>=pos0);
	haplotypes.clear();
	bool found=false;
	for (map<Haplotype, int>::const_iterator it=hb.haplotypes.begin();it!=hb.haplotypes.end();it++) {
		Haplotype newHap=Haplotype(it->first, _start-hb.pos0, _len);
		if (newHap.type==Haplotype::Ref) found=true;
		map<Haplotype, int>::iterator hit=haplotypes.find(newHap);
		if (hit==haplotypes.end()) {
			haplotypes[newHap]=it->second;
		} else {
			if (newHap.type==Haplotype::Ref) hit->first.type=Haplotype::Ref;
			hit->second+=it->second;
		}
		//haplotypes[ Haplotype(it->first, _start-hb.pos0, _len) ]+=it->second;
		// += because subhaplotype may occur multiple times
	}
}

void HapBlock::insert(const Haplotype & seq)
{
	map<Haplotype, int>::iterator hit=haplotypes.find(seq);
	if (hit==haplotypes.end()) {
		haplotypes[seq]=1;
	} else {
		if (seq.type==Haplotype::Ref) hit->first.type=Haplotype::Ref;
		hit->second++;
	}
}

HapBlock::HapBlock(const Haplotype & h, uint32_t start)
{
	pos0=start;
	pos1=start+h.size()-1;
	if (pos1<pos0) {
		cout << pos0 << " " << pos1 << " " << endl;
		cout << "h: " << h << endl;
	}
	assert(pos1>=pos0);
	haplotypes[h]=1;
	type=HapBlock::NORMAL;
}

void HapBlock::setFrequencies()
{
	int sum=0;
	for (map<Haplotype, int>::iterator it=haplotypes.begin();it!=haplotypes.end();it++) {
		sum+=it->second;
	}
	for (map<Haplotype, int>::iterator it=haplotypes.begin();it!=haplotypes.end();it++) {
		(it->first).freq=double(it->second)/double(sum);
	}
}
ostream &operator<<(ostream &stream, const HapBlock &hb)
{
	// construct matrix
	vector<string> output(hb.length());
	vector<int> counts;
	vector<double> freqs;
	for (map<Haplotype, int>::const_iterator it=hb.haplotypes.begin();it!=hb.haplotypes.end();it++)
	{
		for (size_t y=0;y<hb.length();y++) {
			if ((it->first).size()>y) output[y]+=((it->first)[y]); else output[y]+='.';
			output[y]+=' ';
		}
		counts.push_back(it->second);
		freqs.push_back(it->first.freq);
	}

	stream << "start: " << hb.start() << " end: " << hb.end() << " numHap: " << hb.haplotypes.size() << endl;
	for (size_t y=0;y<output.size();y++) cout << output[y] << endl;
	for (size_t y=0;y<counts.size();y++) cout << freqs[y] << " "; cout << endl;
	for (size_t y=0;y<counts.size();y++) cout << counts[y] << " "; cout << endl;
	for (map<Haplotype, int>::const_iterator it=hb.haplotypes.begin();it!=hb.haplotypes.end();it++) cout << it->first.type << " ";
	return stream;
}

bool HapBlock::hasHaplotype(const Haplotype & seq, uint32_t seqStart)
{
	//cout << "hasHaplotype(" << seq << "," << seqStart << "): ";
	for (map<Haplotype, int>::iterator it=haplotypes.begin();it!=haplotypes.end();it++) {
		if (it->first.compare(seqStart-start(), seq.size(), seq)==0) { it->second++; /*cout << "true" << endl;*/ return true; };
	}
	//cout << "false" << endl;
	return false;
}

void HapBlock::showVector(ostream &stream,const vector<HapBlock*> & hapBlocks,uint32_t midPos)
{
	size_t nb=hapBlocks.size();
	vector<size_t> length(nb,0), num(nb,0), pos(nb,0);
	vector<HapBlock*> hbs(nb);
	size_t y=0,x=0,c=0;
	const size_t offset=20;
	size_t indelPos=0;
	for (x=0;x<nb;x++) if (hapBlocks[x]!=NULL){
		pos[c]=offset+y;
		if (midPos>=hapBlocks[x]->start() && midPos<=hapBlocks[x]->end()) indelPos=pos[c];
		length[c]=hapBlocks[x]->length();
		y+=length[c];
		hbs[c]=hapBlocks[x];
		num[c]=hbs[c]->size();
		c++;
	}

	/*
	for (map<int, HapBlock *>::const_iterator it=hb.insertions.begin();it!=hb.insertions.end();it++,x++) {
		pos[c]=y;
		length[c]=it->second->length();
		y+=length[c];
		hbs[c]=it->second;
		num[c]=hbs[c]->size();
		c++;
	}
	*/

	size_t maxLen=*max_element(num.begin(), num.end());
	vector<string> lines(maxLen*2+1,string(offset+y,' '));

	lines[1][1]='R'; lines[1][2]='E'; lines[1][3]='F';
	//for (size_t x=0;x<lines.size();x++) { lines[x][0]='\t'; };
	for (size_t i=0;i<pos.size();i++) {
			//cout << "o: " << o << " o.size() : " << o.size() << " pos[i]: " << pos[i] << endl;
		lines[0][pos[i]]='|';
		/*
		size_t j=1;
		for (map<Haplotype, int>::const_iterator it=hbs[i]->haplotypes.begin();it!=hbs[i]->haplotypes.end();it++) {
			string u=it->first.seq;
			//cout << "u: " << u << endl;
			for (size_t l=0;l<u.size();l++) lines[j][pos[i]+l]=u[l];
			j++;
		}
		j=maxLen+1;
		for (map<Haplotype, int>::const_iterator it=hbs[i]->haplotypes.begin();it!=hbs[i]->haplotypes.end();it++) {
			string o;
			ostringstream os(ostringstream::out);
			os << int(round(-log(it->first.freq))); o=os.str();
			for (size_t l=0;l<o.size();l++) lines[j][pos[i]+l]=o[l];
			j++;
		}
		*/
		// order haplotypes such that reference sequence is top, then sorted based on frequency
		vector<Haplotype> haps; Haplotype refHap;
		for (map<Haplotype, int>::const_iterator it=hbs[i]->haplotypes.begin();it!=hbs[i]->haplotypes.end();it++) if (it->first.type!=Haplotype::Ref) haps.push_back(it->first); else refHap=it->first;
		class SortFunc
		{
		public:
			static bool sortFunc(const Haplotype & h1, const Haplotype & h2)  { return h1.freq<h2.freq; };
		};
		sort(haps.begin(),haps.end(), SortFunc::sortFunc);
		haps.push_back(refHap);

		size_t j=1;
		for (int k=int(haps.size())-1;k>=0;k--) {
			string u=haps[k].seq;
			//cout << "u: " << u << endl;
			for (size_t l=0;l<u.size();l++) lines[j][pos[i]+l]=u[l];
			j++;
		}
		j=maxLen+1;

		for (int k=int(haps.size())-1;k>=0;k--) {
			string o;
			ostringstream os(ostringstream::out);
			os << int(round(-log(haps[k].freq))); o=os.str();
			for (size_t l=0;l<o.size();l++) lines[j][pos[i]+l]=o[l];
			j++;
		}


	}
	lines[0][indelPos]='X';
	for (size_t j=0;j<lines.size();j++) {
		stream << lines[j] << endl;
	}
}
