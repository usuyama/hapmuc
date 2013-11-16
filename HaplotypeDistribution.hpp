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
friend class HDIterator2;
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

class HDIterator2
{
public:
	HDIterator2(const HaplotypeDistribution &hd, size_t maxHap, uint32_t pos, uint32_t left, uint32_t right, int _noIndelWindow=-1)
	{
		// noIndelWindow ignores indels around pos
		noIndelWindow=_noIndelWindow;
		// variants will be added at position pos

		hdPtr=&hd;
		midPos=pos;
		setupBlocks(hd, pos, left, right);
		setThresholds(maxHap);
		init();
	};

	void init()
	{
		for (size_t x=0;x<iter.size();x++) {
			iter[x]=0;
		}
		hap.seq.clear();
		_last=false;
	};
	void operator++()
	{
		size_t x;
		for (x=0;x<iter.size() && (++iter[x])==max[x];++x) {
			iter[x]=0;
			if (x==iter.size()-1) _last=true;
		}
	};
	bool last() const { return _last; };
	uint32_t start() const { return (*hapBlocks.begin())->start(); };
	uint32_t end() const { return (*hapBlocks.rbegin())->end(); };
	Haplotype getMaxFreqHap() const
	{
		Haplotype maxh;
		maxh.seq.clear();
		maxh.freq=1.0;
		maxh.nfreq=1.0;
		for (size_t x=0;x<hbs.size();x++) {
			double mf=0.0;
			size_t idx;
			for (size_t y=0;y<hbs[x].haps.size();y++) if (hbs[x].haps[y].freq>mf) { idx=y; mf=hbs[x].haps[y].freq; };
			maxh.freq*=mf;
			if (!hasIndel[x]) maxh.nfreq*=mf;
			maxh.append(hbs[x].haps[idx].seq);
		}
		return maxh;
	}

	operator Haplotype()
	{
		throw string("REIMPLEMENT");
		hap.seq.clear();
		hap.freq=1.0;
		hap.nfreq=1.0;
		hap.type=Haplotype::Normal;
		//hap.haps.clear();
		for (size_t x=0;x<iter.size();x++) {
			const Haplotype & h=hbs[x].haps[iter[x]];
			// do not append deletions as they are codes by '#'
			if (h.type==Haplotype::In || h.type==Haplotype::Normal) {
				hap.append(h.seq);
			}
			if (h.seq.size()>0) {
				hap.type|=h.type;
			}

			hap.freq*=h.freq;
			if (hasIndel[x]==0) {
				hap.nfreq*=h.freq;
			}
			//hap.haps.push_back(pair<string, double>(h.seq, h.freq));
		}
		return hap;
	}


	double getLogNumHaps() const { return logNumHap; };
	friend ostream &operator<<(ostream &stream, const HDIterator2 & hdi)
	{
		vector<HapBlock *> hb; hb.reserve(hdi.hapBlocks.size());
		for (list<HapBlock*>::const_iterator lit=hdi.hapBlocks.begin();lit!=hdi.hapBlocks.end();lit++) hb.push_back(*lit);
		HapBlock::showVector(stream, hb, hdi.midPos);
		return stream;
	}

	void generateHapsWithAlignedVariants(vector<Haplotype> & haps, const AlignedCandidates & variants, int print=0)
	{
		haps.clear();
	//	if (print) {
			cout << "Variants: ";
			BOOST_FOREACH(Variant var, variants.variants) {
				cout << "[" << var.size() << " " << var.getSeq() << "]";
			}
			cout << endl;
	//	}

		//map <Haplotype, Haplotype> pRef, pInd;
		set <Haplotype> setHap;
		vector <Haplotype> vecHap;
		vector <vector<int> > vecRefPos;

		size_t minLen=100000;
		init();
		while (!last()) {
			hap.seq.clear();
			hap.freq=1.0;
			hap.nfreq=1.0;
			hap.type=Haplotype::Normal;
			//hap.haps.clear();

			vector<int> refPos;
			for (size_t x=0;x<iter.size();x++) {
				const Haplotype & h=hbs[x].haps[iter[x]];
				int len = hbs[x].end-hbs[x].start+1;
				if (hbs[x].type == HapBlock::NORMAL) {
					int p = hbs[x].start;
					bool hasDel = false;
					for (size_t y=0;y<h.seq.size();y++) {
						int c=int(h.seq[y]);
						if (c>=35 && c<65) { hasDel = true; }
						refPos.push_back(p);
						p++;
					}
					if (hasDel == false && int(h.seq.size())!=len) throw string("What's going on here?");
				} else if (hbs[x].type == HapBlock::INSERT) {
					for (size_t y=0;y<h.seq.size();y++) {
						refPos.push_back(-1);
					}
				}
				hap.append(h.seq);
				hap.freq*=h.freq;
			}

			// effectuate deletions at positions outside midPos

			size_t y=0;
			while (y<hap.size()) {
				int c=int(hap[y]);
				if (c>=35 && c<65) {
					int len=c-int('#');
					if (len>int(hap.size()-y)) len=hap.size()-y;
					hap.seq.erase(y,len);
					refPos.erase(refPos.begin()+y,refPos.begin()+y+len);
				} else y++;
			}
			vecHap.push_back(hap);
			vecRefPos.push_back(refPos);
			++(*this);
		}

		// first add variants combinatorially, the add variants to the set of combinatorially generated haplotypes

		for (int ac = 1;ac>=0;ac--) {
            cout << "ac: " << ac << endl;
			size_t numHap = vecHap.size();

			bool addComb = false;
			if (ac==1) {
				addComb = true;
			} else numHap = vecHap.size();
            
			BOOST_FOREACH(AlignedVariant var, variants.variants) {
                cout << "addComb " << addComb << endl;
				if (addComb) {
					numHap = vecHap.size();
				}
                
				if (var.getAddComb()==addComb) {
					for (size_t h=0;h<numHap;h++) {
						Haplotype _hap=vecHap[h];
						bool changed=false;
                        
						cout << "******************************" << endl;
						cout << "var: " << var.getStartHap() << " " << var.getString() << endl;
						cout << " hap: " << vecHap[h].seq << endl;
						vector<int> refPos = vecRefPos[h];
                        BOOST_FOREACH(int x, refPos) {
                            cout << x << " ";
                        }
                        cout << endl;
						vector<int>::iterator it = find(refPos.begin(), refPos.end(), var.getStartHap());
						if (it!=refPos.end()) {

							int idx = distance(refPos.begin(), it);
							if (var.getType()==Variant::DEL) {
								// deletion
								_hap.seq.erase(idx, var.size());
								refPos.erase(refPos.begin()+idx, refPos.begin()+idx+var.size());
								changed=true;
							} else if (var.getType()==Variant::INS) {
								_hap.seq.insert(idx, var.getSeq());
                                refPos.insert(refPos.begin()+idx, (size_t) var.size(), -1);
								changed=true;
							} else if (var.getType()==Variant::SNP) {
								// snp
								const string & seq=var.getSeq();
								char nuc=seq[3];
								if (_hap.seq[idx]!=seq[3]) {
									_hap.seq[idx]=nuc;
									changed=true;
								}
							}
							if (changed) {
								cout << "_hap: " << _hap.seq << endl;
								vecHap.push_back(_hap);
								vecRefPos.push_back(refPos);
							}
						}
					}
				}
			}
		}
		for (size_t x=0;x<vecHap.size();x++) if (vecHap[x].size()<minLen) minLen=vecHap[x].size();

		BOOST_FOREACH(Haplotype hap, vecHap) {
			//setHap.insert(Haplotype(hap,0, minLen));
			setHap.insert(hap);
		}

		BOOST_FOREACH(Haplotype hap, setHap) {
			haps.push_back(hap);
		}
        BOOST_FOREACH(Haplotype hap, haps) {
			cout << hap.seq << " " << hap.indels.size() << " " << hap.snps.size() << endl;
		}
	}

protected:
	void setupBlocks(const HaplotypeDistribution &hd, uint32_t pos, uint32_t left, uint32_t right)
	{
     //  cout << "setupBlocks" << endl;
	//	cout << "_minFreq: " << _minFreq << endl;
		for (size_t x=0;x<hd.hapBlocks.size();x++) if (hd.hapBlocks[x]!=NULL) {
			if (x) {
				if (hd.hapBlocks[x-1]->end()>hd.hapBlocks[x]->start()) {
				//	cout << hd.hapBlocks[x-1]->end() << " " << hd.hapBlocks[x]->start() << endl;
				//	cout << "HD: " << endl << hd << endl;	
					throw string("Blocks are overlapping.");
				}
			}
			if (hd.hapBlocks[x]->start()>=left && hd.hapBlocks[x]->end()<=right) {
				if (hd.hapBlocks[x-1]->end()+1!=hd.hapBlocks[x]->start()) {
				//	cout << "NOT CONSECUTIVE" << endl;
				//	cout << hd.hapBlocks[x-1]->end() << " " << hd.hapBlocks[x]->start() << endl;
				//	cout << "HD: " << endl << hd << endl;	
	
					throw string("Blocks are not consecutive.");

				}

				hapBlocks.push_back(hd.hapBlocks[x]);
				//cout << *hd.hapBlocks[x] << endl;
			}
		}

		// insertions

		list<HapBlock*>::iterator lit=hapBlocks.begin();
		for (map<int, HapBlock*>::const_iterator it=hd.insertions.begin();it!=hd.insertions.end();it++) {
			if (it->second->start()>=left) {
				for (list<HapBlock*>::iterator lit2=lit;lit2!=hapBlocks.end();lit2++) {
					if (int((*lit2)->start())>=it->first)  {
						hapBlocks.insert(lit2, it->second);
						lit=lit2;
						break;
					}
				}
			}
		}

		// copy

		bool found=false;
		hbs.resize(hapBlocks.size());
		hasIndel.resize(hapBlocks.size());
		int x=0;
		for (lit=hapBlocks.begin();lit!=hapBlocks.end();lit++,x++) {
			uint32_t bs=(*lit)->start();
			uint32_t be=(*lit)->end();
			if (pos>=bs && pos<=be) {
				indelIdx=x;
				indelOffs=pos-bs;
				found=true;
			}
			hasIndel[x]=0;
			for (map<Haplotype,int>::const_iterator it=(*lit)->haplotypes.begin();it!=(*lit)->haplotypes.end();it++) {
				hbs[x].haps.push_back(it->first);
			}
			hbs[x].start=bs;
			hbs[x].end=be;
			hbs[x].type=(*lit)->getType();
			if (hbs[x].type==HapBlock::INSERT) hbs[x].end=hbs[x].start-1;
			// set frequency of reference haplotype in block
			// this makes sure that the reference haplotype is always included
			bool reffound=false;
			for (size_t y=0;y<hbs[x].haps.size();y++) {
				if (hbs[x].haps[y].type==Haplotype::Ref) {
					reffound=true;
					//hbs[x][y].freq=1.0;
				} else {
					//for (size_t z=0;z<hbs[x].haps[y].seq.size();z++) hbs[x].haps[y].seq[z] = tolower(hbs[x].haps[y].seq[z]);
				}
			}
			if (!reffound) {
				//cout << **lit << endl;

			}
			assert(reffound==true);
		}

		if (hbs.size() == 0) {
			throw string("Not enough blocks.");
		}

		//if (!found) throw string("Cannot find position of indel in haplotypedistribution.");
	//	cout << "maxFreqHap: " << getMaxFreqHap() << endl;


	};

	void setThresholds(size_t maxHap)
	{
		// hasIndel is currently set to zero for all blocks, because HaplotypeDistribution
		// does not includes indels at midPos
		// get lowest frequency
		vector<double> minFreq(hbs.size(),0.0);
		vector<int> elim(hbs.size(),1);
		size_t x=0;

		typedef vector<Haplotype>::iterator LHIt;
		LHIt it;

		double logMinHap=0.0;
		double logNH=0.0;

		for (x=0;x<hbs.size();x++) {
			logNH+=log(double(hbs[x].haps.size()));
		}




		double logMH=log(double(maxHap));
		if (logMH<logMinHap) logMH=logMinHap;


		// keep removing haplotypes until we have the desired number of haplotypes
		bool erased=true;
		while (logNH>logMH && erased) {
			erased=false;
			for (x=0;x<hbs.size();x++) {
				double mf=2.0;
				for (it=hbs[x].haps.begin();it!=hbs[x].haps.end();it++) {
					if (it->type!=Haplotype::Ref && it->freq<mf) mf=it->freq;
				}
				if (hbs[x].haps.size()<=1) { minFreq[x]=2.0; elim[x]=0; } else minFreq[x]=mf;
			}

			vector<double>::iterator mel=min_element(minFreq.begin(), minFreq.end());
			assert(mel!=minFreq.end());
			size_t y=distance(minFreq.begin(),mel);

			if (elim[y]==0) break;
			// erase the element
			for (it=hbs[y].haps.begin();it!=hbs[y].haps.end();it++) if (it->type!=Haplotype::Ref && it->freq<=*mel) {
				hbs[y].haps.erase(it);
				erased=true;
				break;
			}

			logNH=0.0;

			for (x=0;x<hbs.size();x++) {
				logNH+=log(double(hbs[x].haps.size()));

			}
			//cout << "logNH: " << logNH << " logMH: " << logMH << endl;
		}
		max.resize(hbs.size(),0);
		iter.resize(hbs.size(),0);
		for (x=0;x<hbs.size();x++) max[x]=hbs[x].haps.size();

		logNumHap=logNH;

		// check if we still have the reference sequence in every block
		for (size_t x=0;x<hbs.size();x++) {
			// cout << "hbs[" << x << "]: " << hbs[x].haps.size() << endl;
			bool reffound=false;
			for (size_t y=0;y<hbs[x].haps.size();y++) {
				if (hbs[x].haps[y].type==Haplotype::Ref) {
					reffound=true;
				}
			}
			if (!reffound) {
				cout << "x: " << x << endl;

			}
			if (!reffound) { throw string("Cannot find reference sequence."); };
		}
	}

	double logNumHap;
	bool _last;
	vector<int > iter, max;
	vector<size_t> hasIndel;
	list<HapBlock *> hapBlocks;
	Haplotype hap;
	vector<HDHapBlock > hbs;
	int indelIdx, indelOffs, noIndelWindow;
	uint32_t midPos;

	typedef list<HapBlock*>::iterator HBIt;
	const HaplotypeDistribution *hdPtr;
};

#endif
