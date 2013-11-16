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
#ifndef READ_HPP_
#define READ_HPP_
#include <cmath>
#include "Haplotype.hpp"
#include "bam.h"
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <ext/hash_map>
#include "MyBam.hpp"

using namespace std;

class Read
{

public:
	class FetchReadData
	{
	public:
		FetchReadData(vector<Read> * _reads,int _poolID, vector<MyBam *> * _myBams, int _numReads = 0, int _maxNumReads = 100000)
		{
			reads=_reads;
			poolID=_poolID;
			myBams = _myBams;
			numReads = _numReads;
			maxNumReads = _maxNumReads;
		}
		vector<Read> * reads;
		vector<MyBam *> * myBams;
		int poolID;
		int numReads;
		int maxNumReads;
	};

	Read()
	{
		init(0);
	}
	Read(uint32_t _pos)
	{
		init(pos);
	}
	void init(uint32_t _pos)
	{
		pos=_pos;
		posStat.first=double(_pos);
		posStat.second=1.0;
		initBam=false;
		onReverseStrand=false;
		poolID=-1;
		mateLen = -1;
		matePos = -1;
	}
	Read(const Read & r)
	{
		initBam=false;
		copy(r,0);
	}
	Read(const Read & r, int poolID)
	{
		initBam=false;
		copy(r,0);
	}
	void copy(const Read &r, int x)
	{
        seq_name=r.seq_name;
	 	 seq=r.seq;
		 pos=r.pos;
		 qual=r.qual;
		 posStat=r.posStat;
		 mapQual=r.mapQual;
		 onReverseStrand=r.onReverseStrand;
		 poolID=r.poolID;
		 matePos = r.matePos;
		 mateLen = r.mateLen;
		 bamHeader = r.bamHeader;
		if (initBam) {
			delete[] bam->data;
			delete bam;
			initBam=false;
		}

		 if (r.initBam) {
 			bam=new bam1_t;
			*bam=*r.bam;
			bam->data=new uint8_t[r.bam->m_data];
			bam->m_data=r.bam->m_data;
			for (int m=0;m<r.bam->m_data;m++) bam->data[m]=r.bam->data[m];
			initBam=true;
		}
	}
	Read & operator=(const Read & r)
	{
		if (&r!=this) {
			copy(r,1);
		}
		return *this;
	}
	Read(const bam1_t *b, int _poolID, bam_header_t * _bamHeader, const string & overrideLibName = string("") )
	{
		const bam1_core_t *c=&b->core;
		uint32_t len=c->l_qseq;
		double mapPhred=(double) c->qual;
		mapQual=(1.0-pow(10.0, -mapPhred/10.0));
		if (mapQual<0.0 || mapQual>1.0 || isnan(mapQual) || isinf(mapQual)) throw string("Phred error.");
		if (mapQual<1e-16) mapQual=1e-16;
		if (mapQual>1-1e-16) mapQual=1-1e-16;
       // cout << "init read" << endl;
	//	cout << "mapPhred: " << mapPhred << "mapqual: " << mapQual << endl;

		pos=c->pos;

        seq_name = (string)reinterpret_cast<char*>(bam1_qname(b));
		seq.reserve(len);
		qual.reserve(len);
		for(size_t x=0;x<len;x++) {
			seq+=( bam_nt16_rev_table[ bam1_seqi(bam1_seq(b), x) ] );

			// convert phred to probability
			double basePhred=(double)  ( ( (uint8_t*) bam1_qual(b))[x] );

			double q=(1.0-pow(10.0, -basePhred/10.0));
			if (q<0.0 || q>1.0 ||isnan(q) || isinf(q)) throw string("Phred error.");
			if (q<1e-16) q=1e-16;
			if (q>1.0-1e-16) q=1.0-1e-16;
			qual.push_back( q ); // base quality is on log10 scale
            //cout << basePhred << " " << q << " ";
		}
       // cout << endl;
		posStat=computePositionStatistics(b);

		bam = new bam1_t;
		*bam=*b;
		bam->data=new uint8_t[b->m_data];
		bam->m_data=b->m_data;
		for (int m=0;m<b->m_data;m++) bam->data[m]=b->data[m];
		initBam=true;

		if (bam->core.flag & BAM_FREVERSE)  onReverseStrand=true; else onReverseStrand=false;
		poolID=_poolID;
		matePos = bam->core.mpos;
		mateLen = -1;

		this->bamHeader = _bamHeader;
        /*
		LibraryCollection::const_iterator it;
		if (overrideLibName.empty()) {
			it = libraries.find( this->getLibraryName() );
		} else {
			it = libraries.find( overrideLibName );
		}

		if (it == libraries.end()) {
			deleteBam();
			initBam = false;
			throw string("Cannot find library: ").append(this->getLibraryName());
		} else {
			library = (const Library *) & (it->second);
		}*/
	}
	uint32_t getEndPos() const
	{
		return bam->core.n_cigar? bam_calend(&bam->core, bam1_cigar(bam)) : bam->core.pos + 1;
	}
    /*
	string getLibraryName() const
	{
		if (this->isPaired()) {
			const char *p = bam_get_library((bam_header_t *) this->bamHeader, this->bam);
			if (p) {
				return string(p);
			} else {
				return string("dindel_default");
			}
		} else {
			return string("single_end");
		}
	}*/

	int32_t getBAMMatePos() const { return bam->core.mpos; }
	bool isUnmapped() const { return (bam->core.flag & BAM_FUNMAP) != 0 ; }
	bool mateIsUnmapped() const { return (bam->core.flag & BAM_FMUNMAP) != 0; }
	bool mateIsReverse() const { return (bam->core.flag & BAM_FMREVERSE) != 0; }
	bool isReverse() const { return (bam->core.flag & BAM_FREVERSE) != 0; }
	bool isPaired() const { return (bam->core.flag & BAM_FPAIRED) != 0; }
	void complement()
	{
		for (size_t s=0;s<this->seq.seq.size();s++) {
			char & nuc = this->seq.seq[s];
			if (nuc == 'A') nuc = 'T';
			else if (nuc == 'T') nuc = 'A';
			else if (nuc == 'C') nuc = 'G';
			else if (nuc == 'G') nuc = 'C';
		}
	}
	void reverse()
	{
		string newseq = this->seq.seq;
		size_t len = newseq.size();
		for (size_t x=0;x<newseq.size();x++) newseq[len-x-1]=this->seq.seq[x];
		this->seq.seq = newseq;
	}

	string getAuxData() const
	{
		stringstream os;
		uint8_t *s = bam1_aux(bam);

		while (s < bam->data + bam->data_len) {
			uint8_t type, key[2];
			key[0] = s[0]; key[1] = s[1];
			s += 2; type = *s; ++s;
			//printf("\t%c%c:", key[0], key[1]);
			os << "\t" << key[0] << key[1];
			/*
			if (type == 'A') { printf("A:%c", *s); ++s; }
			else if (type == 'C') { printf("i:%u", *s); ++s; }
			else if (type == 'c') { printf("i:%d", *s); ++s; }
			else if (type == 'S') { printf("i:%u", *(uint16_t*)s); s += 2; }
			else if (type == 's') { printf("i:%d", *(int16_t*)s); s += 2; }
			else if (type == 'I') { printf("i:%u", *(uint32_t*)s); s += 4; }
			else if (type == 'i') { printf("i:%d", *(int32_t*)s); s += 4; }
			else if (type == 'f') { printf("f:%g", *(float*)s); s += 4; }
			else if (type == 'Z' || type == 'H') { printf("%c:", type); while (*s) putchar(*s++); ++s; }
			*/
			if (type == 'A') { os << "A:"<<(char)*s; ++s; }
			else if (type == 'C') { os << "i:" << (unsigned int) *s; ++s; }
			else if (type == 'c') { os << "i:" << (int) *s; ++s; }
			else if (type == 'S') { os << "i:" << *(uint16_t*)s; s += 2; }
			else if (type == 's') { os << "i:" << *(int16_t*)s; s += 2; }
			else if (type == 'I') { os << "i:" << *(uint32_t*)s; s += 4; }
			else if (type == 'i') { os << "i:" <<  *(int32_t*)s; s += 4; }
			else if (type == 'f') { os << "f:" << *(float*)s; s += 4; }
			else if (type == 'Z' || type == 'H') { os << type <<":"; while (*s) os << (char) (*s++); ++s; }
		}
		return os.str();
	}

	// compute mean and standard deviation of first base position
	static pair<double, double> computePositionStatistics(const bam1_t *b)
	{
		const bam1_core_t *c=&b->core;
		uint32_t *cigar=bam1_cigar(b);
		uint32_t k;
		int32_t pos=0, mean=0,totLen=0;

		uint32_t refPos = c->pos;
		double var=0.0;
		if (c->n_cigar==0) {
			return pair<double,double>(-1.0,-1.0);
		}

		for (k = 0; k < c->n_cigar; ++k) {
			int op = cigar[k] & BAM_CIGAR_MASK;
			int32_t len=cigar[k] >> BAM_CIGAR_SHIFT;

			if (op==BAM_CMATCH) {
				mean+=len*(pos-totLen);
				totLen+=len;
			}
			// update position for the next cigar
			if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CSOFT_CLIP || op ==BAM_CHARD_CLIP) {
				pos+=len;
			}
		}
		double dmean=double(mean)/double(totLen);

		pos=0;
		totLen=0;
		for (k = 0; k < c->n_cigar; ++k) {
			int op = cigar[k] & BAM_CIGAR_MASK;
			int32_t len=cigar[k] >> BAM_CIGAR_SHIFT;

			if (op==BAM_CMATCH) {
				var+=double(len)*(double(pos-totLen)-dmean)*(double(pos-totLen)-dmean);
				totLen+=len;
			}
			// update position for the next cigar
			if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CSOFT_CLIP || op ==BAM_CHARD_CLIP) {
				pos+=len;
			}
		}
		var=var/double(totLen);
		return pair<double,double>(dmean+double(refPos), var);
	}

	static void filterReads(vector<Read> & reads, int exclude, const string & match)
	{
		vector<Read> filteredReads;
		for (size_t r=0;r<reads.size();r++) {
			string str=reads[r].getAuxData();
			size_t found=str.find(match);
			if (exclude) {
				if (found==string::npos) {
					filteredReads.push_back(reads[r]);
				}
			} else {
				// include if match
				if (found!=string::npos) filteredReads.push_back(reads[r]);
			}
		}
		reads.swap(filteredReads);
	}

	static int fetchFuncVectorPooled(const bam1_t *b, void *data)
	{
		FetchReadData *ptr=(FetchReadData *) data;

		if (!( (b->core.flag & BAM_FDUP) || (b->core.flag & BAM_FQCFAIL) )) {
			try {
				ptr->reads->push_back(Read(b, ptr->poolID, (*(ptr->myBams))[ptr->poolID]->bh));
				ptr->numReads++;
			} catch (string s) {
                /*
				if (s.find("Cannot find library")!=string::npos) {
					string lib = s.substr(22, s.size()-22);
					string_hash<int>::iterator _it = ptr->unknownLib.find(lib);
					if (_it == ptr->unknownLib.end()) {
						 ptr->unknownLib[lib] = 0;
					} else _it->second++;
				}
				ptr->numUnknownLib++;
				ptr->reads->push_back(Read(b, *(ptr->libraries), ptr->poolID, (*(ptr->myBams))[ptr->poolID]->bh, string("single_end")));
				ptr->numReads++;
                 */
			}
		}
		if (ptr->numReads > ptr->maxNumReads) {
			throw string("Too many reads in region");
		}
		if (ptr->numReads % 10000 == 9999) cout << "numreads: " << ptr->numReads << endl;
		return 0;
	}

	friend ostream &operator<<(ostream &stream, const Read & read)
	{
		//cout << "pos: " << read.pos << " 1-mapping quality: " << 1.0-read.mapQual << " ";
		//for (size_t b=0;b<read.seq.size();b++) stream << read.seq[b];
		//for (size_t b=0;b<read.qual.size();b++) stream << " " << read.qual[b];
		return stream;
	};
	bam1_t * getBam() const { return bam; };

	size_t size() const { return seq.size(); };
	void setAllQual(double v) { qual.clear(); qual.reserve(seq.size()); for (size_t x=0;x<seq.size();x++) qual.push_back(v); };

    string seq_name;
	Haplotype seq;
	vector<double> qual;
	pair<double,double> posStat;
	// offset of read with respect to some reference position

	int32_t pos, matePos, mateLen;

	double mapQual;
	bool initBam;
	bool onReverseStrand;
	int poolID;

	bam_header_t * bamHeader;

	bam1_t *bam;

	void deleteBam()
	{
		delete[] bam->data; delete bam;
	}

	~Read()
	{
		if (initBam) {
			deleteBam();
		}
	}
};

#endif /*READ_HPP_*/
