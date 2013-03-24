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
#include <stdint.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "HaplotypeDistribution.hpp"

using namespace std;




HaplotypeDistribution::HaplotypeDistribution(uint32_t _midPos, const string & refSeq, uint32_t refSeqStart)
{
	pos0=0;
	pos1=0;
	midPos=_midPos;

	// add reference sequence
	uint32_t bs = 4;
	uint32_t rm = refSeq.size() % bs;
	int add = 1; if (rm==0) add=0;

	for (size_t x=0;x<(refSeq.size()/bs)+add;x++) {
		uint32_t start = refSeqStart+x*bs;
		Haplotype refHap(Haplotype::Ref, refSeq.substr(x*bs,bs));
		insertSeq(refHap, start);
	}
}


ostream &operator<<(ostream &stream, const HaplotypeDistribution &hb)
{
	size_t cnt=0;

	for (size_t x=0;x<hb.hapBlocks.size();x++) if (hb.hapBlocks[x]!=NULL){
		stream << "HAPLOTYPE BLOCK " << cnt++ << endl;
		stream << *hb.hapBlocks[x] << endl;
	}

	cnt=0;

	for (map<int, HapBlock *>::const_iterator it=hb.insertions.begin();it!=hb.insertions.end();it++) {
		stream << "INSERTION " << cnt++ << endl;
		stream << *it->second << endl;
	}
	return stream;
}



int HaplotypeDistribution::fetchFuncInsertRead(const bam1_t *b, void *data)
{
	( (HaplotypeDistribution*) data)->insertRead(b);
	return 0;
}

void HaplotypeDistribution::insertRead(const bam1_t* b)
{
	 if ((b->core.flag & BAM_FMUNMAP) != 0) return;
	/*
	for each cigar operation in read {

		get sequence corresponding to cigar n
		calc starting position in reference
		calc confidence of seq (product of mapping quality and base qualities)
		make haplotype-> seq
		insert_seq(seq, reference_pos)
	*/

	const bam1_core_t *c=&b->core;
	uint32_t *cigar=bam1_cigar(b);
	uint32_t k, l=0;
	uint32_t refPos = c->pos;
	int lastop=-1;
	uint32_t lastPos=refPos;
	//cout << "read: " << bam1_qname(b) << endl;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int32_t len=cigar[k] >> BAM_CIGAR_SHIFT;
	//	cout << "cigar" << k << endl;

		Haplotype seq; seq.reserve(len);
		seq.conf=(double) c->qual; // this scales the confidence from the individual base calls with the mapping confidence

		if (op==BAM_CINS || op==BAM_CMATCH || op==BAM_CSOFT_CLIP) {
			for(int32_t x=0;x<len;x++) {
				seq+=( bam_nt16_rev_table[ bam1_seqi(bam1_seq(b), l) ] );
				seq.conf+=(double) bam1_qual(b)[l]; // base quality is on log10 scale
				l++;
			}
		} else if (op==BAM_CDEL) {
			seq.insert(0, len, '#');
	//		cout << "INSDELETION" << endl;
		}


		if (op==BAM_CINS) seq.type=Haplotype::In; else if (op==BAM_CDEL) seq.type=Haplotype::Del; else seq.type=Haplotype::Normal;
		//cout << endl << " *** " << endl;
		//cout << "op: " << op  << " len: " << len << endl;
		//cout << "CIGARseq: " << seq << endl;
		//cout << "refPos " << refPos << endl;

		// now add it to the haplotype structure for this location
		if (seq.size()) {
			if (1) { //refPos!=midPos) {
				//cout << bam1_qname(b) << " refpos: " << refPos << " seq: " << seq << endl;
				if (seq[0]=='#') {
					// deletion, recode

					if (seq.size()>30) {
						cerr << string("Deletion is too long...") << endl;
						len = 30;

					}
					seq.seq.clear();
					seq.seq+=(char(int('#')+len));
				}
				insertSeq(seq, refPos);
			}
			else {
				if (op==BAM_CINS || op==BAM_CDEL) {
					indelsAtMidPos.insert(seq);
				} else insertSeq(seq, refPos);
			}
		}

		// check if previous seq was not an insert
		if (lastop!=-1 && lastop!=BAM_CINS) {
			if (lastPos==refPos && lastop!=BAM_CSOFT_CLIP && lastop!=BAM_CHARD_CLIP) throw string("Mag niet.");
			for (uint32_t pos=lastPos;pos<refPos;pos++) {
				map<int, HapBlock*>::iterator it=insertions.find(pos);
				if (it!=insertions.end()) {
					(it->second)->insert(Haplotype(Haplotype::In));
				}
			}
		}

		// update position for the next cigar
		lastPos=refPos;
		if (op == BAM_CMATCH || op == BAM_CDEL || op==BAM_CREF_SKIP) {
			refPos+=(uint32_t) len;
		} else if (op!=BAM_CINS && op != BAM_CSOFT_CLIP && op != BAM_CHARD_CLIP) throw string("I don't know how to smoke this CIGAR");
		lastop=op;
	}


}

vector<Variant> HaplotypeDistribution::getIndelVariantsAtMidPos()
{
	vector<Variant> variants;
	BOOST_FOREACH(Haplotype hap, indelsAtMidPos) {
		if (hap.type==Haplotype::In) {
			variants.push_back(Variant(string("+").append(hap.seq)));
		} else if (hap.type==Haplotype::Del) {
			variants.push_back(Variant(string("-").append(string(hap.seq.size(),'R'))));
		} else throw string("Unrecognized variant");
	}
	return variants;
}

void HaplotypeDistribution::setFrequencies()
{
	for (size_t x=0;x<hapBlocks.size();x++) {
		HapBlock * hb=hapBlocks[x];
		if (hb!=NULL) {
			hb->setFrequencies();
		}
	}
	// insertions

	for (map<int, HapBlock*>::iterator it=insertions.begin();it!=insertions.end();it++) {
		HapBlock * hb=it->second;
		if (hb!=NULL) {
			hb->setFrequencies();
		}
	}

}

void HaplotypeDistribution::updateBlock(HapBlock *hb, const Haplotype & seq, uint32_t seqStart)
{
	if (seq.size()!=hb->length() || seqStart!=hb->start()) throw string("updateBlock-seq mismatch.");
	hb->insert(seq);
}





bool sortFunc(const HapBlock *a, const HapBlock *b)
{
	if (a==NULL && b!=NULL) return false;
	else if (a!=NULL && b==NULL) return true;
	else if (a==NULL && b==NULL) return false;
	else {
		if (a->start()<b->start()) return true; else return false;
	}
}

void HaplotypeDistribution::check()
{
	sort(hapBlocks.begin(), hapBlocks.end(), sortFunc);
	vector<HapBlock*>::iterator it=find(hapBlocks.begin(), hapBlocks.end(),(HapBlock *) NULL );
	for (vector<HapBlock*>::iterator it2=it;it2!=hapBlocks.end();it2++) {
		if (*it2!=NULL) throw string("Error: NULLs not consecutive");
	}
	for (size_t x=0;x<hapBlocks.size();x++) {
		if (hapBlocks[x]->end()<hapBlocks[x]->start()) {
			cout << "CHECK SMALLER HD:" << endl;
			cout << *this << endl;

			throw string("Blocks are smaller!");
		}
	}
	for (size_t x=1;x<hapBlocks.size();x++) {
		if (hapBlocks[x-1]->end()+1!=hapBlocks[x]->start()) {
			cout << "CHECK CONSECUTIVE HD:" << endl;
			cout << *this << endl;

			throw string("Blocks are not consecutive!");
		}
	}

	for (size_t x=1;x<hapBlocks.size();x++) {
		if (hapBlocks[x-1]->end()>=hapBlocks[x]->start()) {
			cout << "CHECK HD:" << endl;
			cout << *this << endl;

			throw string("Blocks are overlapping!");
		}
	}
}


void HaplotypeDistribution::newBlock(HapBlock *hb)
{
	//checkBlock(hb);
	vector<HapBlock*>::iterator it=find(hapBlocks.begin(), hapBlocks.end(),(HapBlock *) NULL );
	if (it==hapBlocks.end()) {
		hapBlocks.push_back(hb);
	} else {
		*it=hb;
	}
	if (hb->start()<pos0) pos0=hb->start();
	if (hb->end()>pos1) pos1=hb->end();

	sort(hapBlocks.begin(), hapBlocks.end(), sortFunc);
}


/*
void HaplotypeDistribution::invalidateBlock(HapBlock *hb)
{
	for (uint32_t p=hb->start();p<=hb->end();p++) {
		posToBlock[p-pos0]=-1;
	}
}
*/
void HaplotypeDistribution::deleteBlock(int idx)
{
	delete hapBlocks[idx];
	hapBlocks[idx]=NULL;
	if (hapBlocks.size()>1 && idx!=(int)hapBlocks.size()-1) {
			hapBlocks[idx]=hapBlocks.back();
			hapBlocks.back()=NULL;
	}
}


void HaplotypeDistribution::splitBlock(int idx, const Haplotype & seq, uint32_t seqStart)
{
	// block **********
	// seq	    ***

	if (seq.size()==0) throw string("Empty haplotype!");

	uint32_t seqEnd=seqStart+seq.size()-1;
	HapBlock & block=*hapBlocks[idx];
	if (seqStart<block.start()||seqStart+seq.size()-1>block.end()) throw string("seq outside of block boundaries");

	uint32_t lenA=seqStart-block.start();
	uint32_t lenB=seq.size();
	if (lenB==0) throw string("Empty sequence!");
	uint32_t lenC=(block.end()==seqStart+seq.size()-1) ? 0 : block.end()-seqEnd;
	// cout << "block.start: " << block.start() << " block.end " << block.end() << " seqStart: " << seqStart << " seqEnd: " << seqStart+seq.size()-1 << " lenA: " << lenA << " lenB: " << lenB << " lenC: " << lenC << endl;
	if (1) { //!block.hasHaplotype(seq, seqStart)) {
		// split them
		// note that blocks are not overlapping, so we can append them to the vector
		HapBlock *hbA, *hbB, *hbC;
		if (lenA) hbA=new HapBlock(block, block.start(), lenA);
		if (lenB) hbB=new HapBlock(block, block.start()+lenA, lenB);



		if (lenC) hbC=new HapBlock(block, hbB->end()+1, lenC);

		deleteBlock(idx);

		if (lenB) newBlock(hbB);
		updateBlock(hbB, seq, seqStart);

		//if (lenA) cout << "hbA: " << *hbA << endl;
		//if (lenB) cout << "hbB: " << *hbB << endl;
		//if (lenC) cout << "hbC: " << *hbC << endl;

		if (lenA) newBlock(hbA);
		if (lenC) newBlock(hbC);
	}



}

int HaplotypeDistribution::getFirstOverlappingBlock(uint32_t seqStart, uint32_t seqEnd) const
{
	size_t x=0;
	while (x<hapBlocks.size()&&hapBlocks[x]!=NULL) {
		const HapBlock & hb=*hapBlocks[x];
		if ( (hb.end()>=seqStart) && hb.start()<=seqEnd ) return int(x);
		x++;
	}
	return -1;
}

void HaplotypeDistribution::insertSeq(Haplotype & seq, uint32_t seqStart)
{
	// this->check();
	if (seq.type==Haplotype::Normal || seq.type==Haplotype::Ref || seq.type==Haplotype::Del) {
		//cout << "insertSeq. seq: " << seq << " seqStart: " << seqStart << endl;
		uint32_t seqEnd = seqStart+(uint32_t) seq.size()-1;
		int idx=getFirstOverlappingBlock(seqStart, seqEnd);
		if (idx!=-1) {
			// cout << "found idx: " << idx << endl;
			HapBlock & block=*hapBlocks[idx];
			if (block.start()<seqStart) {
				assert(block.end()>=seqStart);
				if (seqEnd>block.end()) {
	//				cout << "typeA" << endl;
					uint32_t overlap=block.end()-seqStart+1;
					// block **********
					// seq          ********
					Haplotype seqA(seq, 0, overlap);
					splitBlock(idx, seqA, seqStart);
					Haplotype seqB(seq, overlap, seq.size());
					insertSeq(seqB, seqStart+overlap);
				} else {
	//				cout << "typeB" << endl;
					// block ***************
					// seq       *****
					splitBlock(idx, seq, seqStart);
				}
			} else
			{
				// block.start() >=seqStart
				if (block.end()>seqEnd) {
					// block.start() >=seqStart && block.end()>seqEnd
	//				cout << "typeC" << endl;
					// block      ***********
					// seq    *******
					uint32_t overlap=seqEnd-block.start()+1;
					assert(overlap>0 && overlap<=seq.size());

					Haplotype seqB(seq, seq.size()-overlap, overlap);

					splitBlock(idx, seqB, block.start());
	//				cout << "seqB: " << seqB << endl;

					if (overlap<seq.size()) {
						Haplotype seqA(seq, 0, seq.size()-overlap);
						// cout << "seqA: " << seqA << endl;
						newBlock(new HapBlock(seqA, seqStart));
					}

					// note that newBlock invalidates idx!

				} else
				{
					// block.start() >=seqStart && block.end()<=seqEnd
	//				cout << "typeD" << endl;
					// block      ********* ***
					// seq      ******************
					uint32_t lenA=block.start()-seqStart;
					uint32_t lenB=block.end()-block.start()+1;
					uint32_t lenC=seqEnd-block.end();
					if (lenA) {
						Haplotype seqA(seq, 0, lenA);
						newBlock(new HapBlock(seqA, seqStart));
					}

					assert(lenB>0);
					Haplotype seqB(seq, lenA, lenB);
					updateBlock(&block, seqB, seqStart+lenA);

					if (lenC) {
						Haplotype seqC(seq, lenA+lenB, lenC);
						insertSeq(seqC, seqStart+lenA+lenB);
					}
				}
			}

		} else {
			// cout << "newBLock" << endl;
			newBlock(new HapBlock(seq, seqStart));
		}
	}
	else if (seq.type==Haplotype::In) {
		map<int, HapBlock *>::iterator it=insertions.find(seqStart);
		if (it==insertions.end()) {
			insertions[seqStart]=new HapBlock(seq, seqStart);
			insertions[seqStart]->setType(HapBlock::INSERT);
			insertions[seqStart]->insert(Haplotype(Haplotype::Ref));//, string(seq.size(),'0')));
		} else
		{
			it->second->insert(seq);
		}
	}
	else throw string("Cannot handle this case.");
}

size_t HaplotypeDistribution::getNumberOfHaplotypes(uint32_t start, uint32_t end) const
{
	size_t n=1;
	size_t x=0;
	while (x<hapBlocks.size() && hapBlocks[x]!=NULL) {
		const HapBlock & hb=*hapBlocks[x];
		if (hb.end()>=start && hb.start()<=end) {
			n*=hapBlocks[x]->size();
		}
		x++;
	}
	for (map<int, HapBlock*>::const_iterator it=insertions.begin();it!=insertions.end();it++) {
		if (it->first>=(int) start && it->first<=(int)end) {
			n*=(it->second->size()); // we consider also haplotypes without the insertion
		}
	}
	return n;
}

size_t HaplotypeDistribution::getNumberOfHaplotypes(uint32_t start, uint32_t end, double minFreq) const
{
	size_t n=1;
	size_t x=0;
	while (x<hapBlocks.size() && hapBlocks[x]!=NULL) {
		const HapBlock & hb=*hapBlocks[x];
		if (hb.end()>=start && hb.start()<=end) {
			size_t nh=0;
			for (map<Haplotype, int>::const_iterator hit=hb.haplotypes.begin();hit!=hb.haplotypes.end();hit++) {
				const Haplotype & h=hit->first;
				if (h.type==Haplotype::In || h.type==Haplotype::Del) nh++;
					else { if (h.freq>minFreq) nh++; }
			}
			n*=nh;
		}
		x++;
	}
	for (map<int, HapBlock*>::const_iterator it=insertions.begin();it!=insertions.end();it++) {
		if (it->first>=(int) start && it->first<=(int)end) {
			n*=(it->second->size()); // we consider also haplotypes without the insertion
		}
	}
	return n;
}

HaplotypeDistribution::~HaplotypeDistribution()
{
	for (size_t x=0;x<hapBlocks.size();x++) if (hapBlocks[x]!=NULL) delete hapBlocks[x];
	for (map<int,HapBlock*>::iterator it=insertions.begin();it!=insertions.end();it++) delete it->second;
}
