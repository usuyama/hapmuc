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
#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_
#include <stdint.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include "Variant.hpp"
#include "MLAlignment.hpp"
#include <boost/foreach.hpp>
#include <ext/hash_map>
#include <set>
//#include "Fast.hpp"
using namespace std;
using __gnu_cxx::hash;
namespace std { using namespace __gnu_cxx; }

const char NUCLEOTIDES[]={'A','T', 'G','C'};





class Haplotype //: public string
{
public:
	// ContainsInDel means the haplotype is contains a small non-zero length segment
	// that was identified as an InDel from Cigar by WH alignment
	//typedef enum { Normal, In, Del, HasIn, HasDel, HasInDel } Type;
	typedef int Type;
	static const int Ref=1;
	static const int Normal=Ref<<1;
	static const int In=Ref<<2;
	static const int Del=Ref<<3;



	mutable Type type;
	// nfreq is the product of the frequencies of haplotypes that are not indels
	mutable double freq, conf, nfreq;
	uint32_t pos;
	string seq;
	string indel; // if haplotype has indel for a given position
	string align; // annotates for each base in the sequence this haplotype was aligned to whether the equal to the reference=R, snp=S, deletion=D, insertions cannot be recorded this way
	//vector<pair<string, double> > haps;
	map<int, AlignedVariant  > indels, snps;
	MLAlignment ml;

	size_t size() const { return seq.size(); };
	char & operator[](size_t idx) {  return seq[idx]; };
	const char & operator[](size_t idx) const {  return seq[idx]; };
	Haplotype & operator+=(char c) { seq+=c; return *this; };
	bool operator<(const Haplotype & h) const { return seq<h.seq; };

	/*
	bool operator<(const Haplotype & h) const
	{
		if (seq!=h.seq) {
			if (type<h.type) return true;
			else if (type==h.type) return seq<h.seq;
		} else return seq<h.seq;
        };
	*/


	int compare ( size_t pos1, size_t n1, const Haplotype & h ) const { return seq.compare(pos1,n1,h.seq); };
	Haplotype & insert ( size_t pos1, size_t n, char c ) { seq.insert(pos1,n,c); return *this; };
	void reserve(size_t n) { seq.reserve(n); };
	Haplotype & append(const string & str) { seq.append(str); return *this; };

	Haplotype(const Haplotype & h, size_t pos0, size_t n)
	{
		seq=h.seq.substr(pos0, n);
		conf=h.conf;
		freq=h.freq;
		type=h.type;
		nfreq=h.nfreq;
		indel=h.indel;
		align=h.align;
		pos=h.pos;
		snps=h.snps;
		indels=h.indels;
		ml=h.ml;
		//haps=h.haps;
	};
	Haplotype()
	{
		type=Normal;
		conf=0.0;
		freq=0.0;
		nfreq=0.0;
		pos=0;

	};
	Haplotype(Type _type)
	{
		type=_type;
		conf=0.0;
		freq=0.0;
		nfreq=0.0;
		pos=0;
	}
	Haplotype(Type _type, const string & _seq )
	{
		seq=_seq;
		type=_type;
		conf=0.0;
		freq=0.0;
		nfreq=0.0;
		pos=0;
	}
	Haplotype(const Haplotype &h)
	{
		seq=h.seq;
		conf=h.conf;
		freq=h.freq;
		type=h.type;
		nfreq=h.nfreq;
		indel=h.indel;
		align=h.align;
		pos=h.pos;
		snps=h.snps;
		indels=h.indels;
		ml=h.ml;
		//haps=h.haps;
	}

	Haplotype & operator=(const Haplotype & h)
	{
		if (&h!=this) {

			seq=h.seq;
			conf=h.conf;
			freq=h.freq;
			type=h.type;
			nfreq=h.nfreq;
			indel=h.indel;
			pos=h.pos;
			snps=h.snps;
			indels=h.indels;
			align=h.align;
			ml=h.ml;
			//haps=h.haps;
		}
		return *this;
	}

	string getIndel(int relPos) const
	{
		map<int, AlignedVariant>::const_iterator it=indels.find(relPos);
		if (it==indels.end()) {
			char a=align[relPos];
			if (a=='R') return string("*REF"); else return string("RR=>")+=a;
		} else {
			const AlignedVariant & av=it->second;
			//if (av.getType()==Variant::SNP)  throw string("Haplotype::getIndel error");
			return av.getString();
		}
	}

	string getSNP(int relPos) const
	{
		map<int, AlignedVariant>::const_iterator it=snps.find(relPos);
		if (it==snps.end()) {
			char a=align[relPos];
			if (a=='R') return string("*REF"); else return string("RR=>")+=a;
		} else {
			return it->second.getString();
		}
	}
	Haplotype filtered() const
	{
		/*
		Haplotype hap=*this, newhap=*this;
		newhap.seq.clear();
		transform(hap.seq.begin(), hap.seq.end(), hap.seq.begin(), ::toupper);

		for (size_t x=0;x<hap.seq.size();x++) {
			if (hap.seq[x]!='_' && hap.seq[x]!='#') newhap+=hap.seq[x];
		}
		*/
		return *this;
	}

	void addRefVariant(int rp)
	{
		map<int, AlignedVariant>::const_iterator it;

		// first do indels
		int offset=0;
		// get base position in haplotype of rp (relative position in reference sequence)

		bool addVariant=true;
		it=indels.begin();
		while (it!=indels.end() && it->first<=rp) {
			if (it->second.getType()==AlignedVariant::DEL){
				if (it->first+it->second.size()<=rp) {
					offset-=it->second.size();
				} else {
					// deletion deleted rp from reference
					//addVariant=false;
					break;
				}
			}
			if (it->second.getType()==AlignedVariant::INS) offset+=it->second.size();
			it++;
		}

		if (addVariant) {
			int readStart=rp+offset;
			int readEnd=rp+offset;
			int hapStart=rp;
			int hapEnd=rp;

			if (indels.find(rp)==indels.end()) {
				// no indel at position relPos
				string gt;
				char a=align[rp];
             //   cout << "no indel" << a << endl;
				if (a!='R') {
					gt=string("R=>"); gt+=a;
				} else gt=string("*REF");
				indels[rp]=AlignedVariant(gt, hapStart, hapEnd, readStart, readEnd);
			}

			if (snps.find(rp)==snps.end()) {
				// no snp at position relPos
				string gt;
				char a=align[rp];
                //cout << "no snp" << a << endl;

				if (a!='R') {
					gt=string("R=>"); gt+=a;
				} else gt=string("*REF");
				snps[rp]=AlignedVariant(gt, hapStart, hapEnd, readStart, readEnd);
			}
		}
	}


	int countIndels() const
	{
		int num = 0;
		for (map<int, AlignedVariant>::const_iterator it = indels.begin();it!=indels.end();it++) {
			if (it->second.getType() == Variant::INS || it->second.getType() == Variant::DEL) num++;
		}
		return num;
	}

	int countSNPs() const
	{
		int num = 0;
		for (map<int, AlignedVariant>::const_iterator it = snps.begin();it!=snps.end();it++) {
			if (it->second.getType() == Variant::SNP && !it->second.isRef()) num++;
		}
		return num;
	}
	/*
	int getRefPos(int pos) const
	{
		// returns position of base in haplotype with respect to reference it was aligned to
		if (!indels.size()) return pos; else {
			int offset=0;
			map<int, AlignedVariant>::const_iterator it=indels.begin();
			while (it!=indels.end() && pos>it->second.getPos()) {
				offset-=it->second.length;
				it++;
			}
			return pos+offset;
		}
	}
	*/

	friend ostream &operator<<(ostream &stream, const Haplotype &h)
	{
		stream << "type: " << h.type << " seq: " << h.seq << " len: " << h.size() << " nfreq: " << h.nfreq << " freq: " << h.freq << " indel: " << h.indel;
		return stream;
	}

	/*
	void printHaps() const
	{
		cout << "freq: " << nfreq << " length: " << seq.size() << endl;
		for (size_t i=0;i<haps.size();i++) cout << "[" << i << " |" << haps[i].first << "|," << haps[i].second << "]";
		cout << endl;

		for (size_t i=0;i<haps.size();i++) cout << haps[i].first; cout << endl;
		for (size_t i=0;i<haps.size();i++) {
			cout << int(round(-log(haps[i].second)));
			if (haps[i].first.size()>1) cout << string(haps[i].first.size()-1,' ');
		}
		cout << endl;
		cout << seq << endl;

		cout << endl;
	}
	*/

};


class HapHash
{
public:
	HapHash(unsigned int _kmer, const Haplotype & hap )
	{
		kmer=_kmer;
		mask=( 1<< (2*kmer) )-1;
		makeHash(hap);
	}
	unsigned int getKmer() const { return kmer; };
	unsigned int getMask() const { return mask; };

	const set<int> & lookup(const string & seq, int pos) const {
		int v=convert(seq,pos);
		Hash::const_iterator it=hash.find(v);
		if (it==hash.end()) return emptySet; else return it->second;
	}

	inline const set<int> & lookup(unsigned int key) const
	{
		Hash::const_iterator it=hash.find(key);
		if (it==hash.end()) return emptySet; else return it->second;

	}
	inline unsigned int convert(const string & seq, int pos) const
	{
		if (pos+kmer>seq.size()) throw string("HapHash string too short");
		int v=0;
		for (int x=pos, y=0;x<int(pos+kmer);x++,y++) {
			v |= (map_char(seq[x]) << (2*y) );
		}
		return v;
	}
	inline unsigned int pushBack(const unsigned int & key, const char & c) const
	{
		return (key >> 2) | (map_char(c) << (2*(kmer-1)));
	}
	void print() {
		for (Hash::const_iterator it=hash.begin();it!=hash.end();it++) {
			//cout << "hash[" << it->first << "]: ";
			BOOST_FOREACH(int i, it->second) {
			//	cout << " " << i;
			}
			//cout << endl;

		}


	}
	inline int map_char(const char & c) const
	{
		// TODO do something with N's in reads?
		if (c=='A') return 0; else if (c=='C') return 1; else if (c=='G') return 2; else if (c=='T') return 3; else return 0; //throw string("Haplotype/Read in hash has N's");
	}
	typedef hash_map<unsigned int, set<int> > Hash;

protected:
	unsigned int kmer;
	unsigned int mask;
	const Haplotype *hap_ptr;
	Hash hash;
	set<int> emptySet;



	void makeHash(const Haplotype & hap)
	{
		for (int x=0;x<int(hap.size())-int(kmer);x++) hash[convert(hap.seq,x)].insert(x);
	}
};



#endif /*HAPLOTYPE_HPP_*/

