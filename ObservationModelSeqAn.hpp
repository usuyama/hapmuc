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
 * ObservationModelSeqAn.hpp
 *
 *  Created on: Sep 4, 2009
 *      Author: caa
 */

#ifndef OBSERVATIONMODELSEQAN_HPP_
#define OBSERVATIONMODELSEQAN_HPP_

#include <seqan/align.h>
#include <seqan/align/align_base.h>
#include <seqan/graph_align.h>
#include "Haplotype.hpp"
#include "Read.hpp"
#include "ObservationModel.hpp"
#include "MLAlignment.hpp"
using namespace seqan;
using namespace std;

const int DEBUGCONVERTALIGNMENT = 0;

class Realign {
public:
	void getFlankingCoordinatesBetter(const Haplotype & hap, const Read & read, AlignedVariant & av)
	{
		int rightFlankHap, leftFlankHap,rightFlankRead, leftFlankRead;
		//cout << "variant: " << av.getString() << endl;
		//cout << "startHap: " << av.getStartHap() << endl;
		//cout << "startRead: " << av.getStartRead() << endl;


		if (av.getType()==Variant::DEL) {
			const string & seq = av.getSeq();
			int l = seq.size();

			string origSeq = hap.seq;

			int sh = av.getStartHap();
			origSeq.erase(sh, l);
			leftFlankHap = sh-1;
			rightFlankHap = sh+l;
			int newpos = sh;
			for (int x=sh-1;x>0;x--) {
				string newseq = hap.seq;
				newseq.erase(x, l);
				if (newseq == origSeq) {
					leftFlankHap = x-1;
					newpos = x;
				}
			}
			if (leftFlankHap<=0) leftFlankHap = 0;
			for (int x=sh+1;x<int(hap.seq.size()-l);x++) {
				string newseq = hap.seq;
				newseq.erase(x, l);
				if (newseq == origSeq) {
					rightFlankHap = x+l;
					newpos = x;
				}
			}

			leftFlankRead = av.getStartRead()- (sh-leftFlankHap)+1; if (leftFlankRead<0) leftFlankRead = 0;
			rightFlankRead = av.getStartRead()+1 + (rightFlankHap-sh-l); if (rightFlankRead>=int(read.seq.size())) leftFlankRead = read.seq.size()-1;



			//cout << "leftFlankHap: " << leftFlankHap << " rightFlankHap: " << rightFlankHap << endl;
			//cout << "leftFlankRead: " << leftFlankRead << " rightFlankRead: " << rightFlankRead << endl;
		} else if (av.getType()==Variant::INS) {
			const string & seq = av.getSeq();
			int l = seq.size();
			string newiseq;
			string origSeq = hap.seq;

			int sh = av.getStartHap();
			origSeq.insert(sh, seq);
			leftFlankHap = sh-1;
			rightFlankHap = sh;
			int newpos = sh;
			for (int x=sh-1;x>0;x--) {
				string newseq = hap.seq;
				string iseq = origSeq.substr(x, l);
				newseq.insert(x, iseq);
				int eq = 0;
				if (newseq == origSeq) {
					leftFlankHap = x-1;
					eq = 1;
					newpos = x;
					newiseq = iseq;
				}
				//cout << "x: " << x << " iseq: " << iseq << " eq: " << eq << endl;

			}
			if (leftFlankHap<=0) leftFlankHap = 0;
			for (int x=sh+1;x<int(hap.seq.size()-l);x++) {
				string newseq = hap.seq;
				string iseq = origSeq.substr(x, l);
				newseq.insert(x, iseq);
				int eq=0;
				if (newseq == origSeq) {
					rightFlankHap = x;
					eq=1;
					newpos = x;
				}
				//cout << "x: " << x << " iseq: " << iseq << " eq: " << eq << endl;

			}

			leftFlankRead = av.getStartRead()- (sh-leftFlankHap)+1; if (leftFlankRead<0) leftFlankRead = 0;
			rightFlankRead = av.getStartRead()+l + (rightFlankHap-sh); if (rightFlankRead>=int(read.seq.size())) leftFlankRead = read.seq.size()-1;



			//cout << "leftFlankHap: " << leftFlankHap << " rightFlankHap: " << rightFlankHap << endl;
			//cout << "leftFlankRead: " << leftFlankRead << " rightFlankRead: " << rightFlankRead << endl;
			//cout << "newiseq: " << newiseq << endl;
		} else {
			leftFlankRead = av.getStartRead()-1; if (leftFlankRead<0) leftFlankRead = 0;
			rightFlankRead = av.getStartRead()+1; if (rightFlankRead>=int(read.seq.size())) leftFlankRead = read.seq.size()-1;
			leftFlankHap = av.getStartHap()-1; if (leftFlankHap<0) leftFlankHap = 0;
			rightFlankHap = av.getStartHap()+1; if (rightFlankHap>=int(hap.seq.size())) leftFlankHap = hap.seq.size()-1;
		}
		av.setFlanking(leftFlankHap, rightFlankHap, leftFlankRead, rightFlankRead);
	}
};


template <typename TSource, typename TSpec>
inline void
convertAlignment(
          Align<TSource, TSpec> const & source, MLAlignment & ml, int hlen, int rlen, const Haplotype & hap, const Read & read)
{
        typedef Align<TSource, TSpec> const TAlign;
        typedef typename Row<TAlign>::Type TRow;
        typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
        typedef typename Position<TAlign>::Type TPosition;

        TPosition begin_ = beginPosition(cols(source));
        TPosition end_ = endPosition(cols(source));

	Realign realign;

	if (DEBUGCONVERTALIGNMENT) cout << "begin_ " << begin_ << "  end_ " << end_ << endl;


        ml.relPos = 0;
        bool fbfound=false;

        TRow& row_ = row(source, 0);
        typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
        TIter begin1_ = iter(row_, begin_);
        TIter end1_ = iter(row_, end_);


	ml.align=string(hlen,'R');
	ml.hpos=vector<int>(rlen,MLAlignment::LO);


	int b=0;
	int hs=0; // relative start of haplo

	int rb=0;

	while (isGap(row(source,0), begin_+b)) {
		ml.relPos--;
		if (!isGap(row(source,1),begin_+b)) {
			ml.hpos[rb]=MLAlignment::LO;
			rb++;
		}
		++b;
	}
	hs=b;

	if (DEBUGCONVERTALIGNMENT) cout << "relpos: " << ml.relPos << " hs: " << hs << endl;

	int hb=0; // number of haplotype bases

	while (begin_+b<end_ && rb<rlen) {
	        if (DEBUGCONVERTALIGNMENT)  cout << "b: " << b << " hb: " <<  hb << endl;
                if (isGap(row(source,0), begin_+b)) {
                        if (hb<hlen) {
				// insertion
				string seq("+");
				TIter it=iter(row(source,1), begin_+b);
				while (isGap(row(source,0),begin_+b) && begin_+b<end_) {
					seq+=convert<char>(*it);
					ml.hpos[rb]=MLAlignment::INS;
					++b;
					++it;
					++rb;
				}
				if (DEBUGCONVERTALIGNMENT) cout << "insertion: " << hb << " seq: " << seq << " readpos: " << rb-1 << " - " << rb-seq.size()+1 << endl;
				AlignedVariant av(seq, hb,hb, rb-seq.size()+1, rb-1);
				realign.getFlankingCoordinatesBetter(hap, read, av);
				ml.indels[hb]=av;
			} else {
				ml.hpos[rb]=MLAlignment::RO;
				++rb;
				++b;
			}

                } else {
                        if (!isGap( row(source,1), begin_+b)) {
                                if (!fbfound) {
                                        fbfound=true;
					ml.firstBase=hb;
                                }
                                if (row(source,1)[begin_+b]!=row(source,0)[begin_+b]) {
                                        // SNP
					string snp("X=>X");
					snp[0]=convert<char>(row(source,0)[begin_+b]);
					snp[3]=convert<char>(row(source,1)[begin_+b]);

					if (DEBUGCONVERTALIGNMENT) cout << "SNP: " << hb << " " << convert<char>(row(source,0)[begin_+b]) << "=>" << convert<char>(row(source,1)[begin_+b]) << endl;
					ml.snps[hb]=AlignedVariant(snp, hb, hb, rb,rb);
					realign.getFlankingCoordinatesBetter(hap, read, ml.snps[hb]);
					ml.align[hb]=snp[3];
                                }
				ml.hpos[rb]=hb;
				++rb;
				++b;
				++hb;
                        } else {
                                // deletion
				string seq("-");
				TIter it=iter(row(source,0), begin_+b);
				int len=0;
				while (isGap(row(source,1),begin_+b) && begin_+b<end_) {
					seq+=convert<char>(*it);
					ml.align[hb]='D';
					++b;
					++it;
					++hb;
					++len;
				}
				if (fbfound) {
					ml.indels[hb-len]=AlignedVariant(seq, hb-len,hb-1, rb-1,rb);
					realign.getFlankingCoordinatesBetter(hap, read, ml.indels[hb-len]);
					if (DEBUGCONVERTALIGNMENT) cout << "deletion: " << hb-len << " - " << hb-1 <<  " seq: " << seq << " readpos: " << rb-1 << " - " << rb << endl;

				}


                        }

        	}
	}
	ml.lastBase=hb;

	if (DEBUGCONVERTALIGNMENT) {
		cout << "mfb: " << ml.firstBase << " " << ml.lastBase << endl;
		for (int r=0;r<rlen;r++) cout << "[" << r << "," << ml.hpos[r] << "]"; cout << endl;
	}

}

/*
int main()
{



	// check
	// seqan::DnaString     _refSeq("ATGGCGTGACTGATCCTATCCCCGTT");
	// seqan::DnaString _hapSeq("TTATATGGCGTG");

	//seqan::DnaString _refSeq("ATGGCGTGACTGATCCTATCGTCGTT");
	//seqan::DnaString   _hapSeq("CCCGGTGACTCC");

	seqan::DnaString _refSeq("ATGGCGTGACTGATCCTATCGTCGTT");
	seqan::DnaString           _hapSeq("CTATCGTCTGTAGGTGTCCT");


	seqan::Score<int> score(-1, -460, -100,-960);

	seqan::Align<seqan::DnaString, seqan::ArrayGaps> align;
	seqan::resize(seqan::rows(align), 2);
	seqan::assignSource(seqan::row(align, 0), _refSeq);
	seqan::assignSource(seqan::row(align, 1), _hapSeq);
	cout << "Score = " << seqan::globalAlignment(align, score) << endl;
	cout << align << endl;

	MLAlignment ml;

	convertAlignment(align,ml, length(_refSeq),length(_hapSeq));

}
*/

class ObservationModelSeqAn
{
public:
	ObservationModelSeqAn(const Haplotype & _hap, const Read & r, uint32_t _hapStart, const ObservationModelParameters & _params, const seqan::Score<int> & _score)
	{
		score =_score;
		hap_ptr = &_hap;
		read_ptr = &r;
		hapStart=_hapStart;
		params=_params;
		aligned=false;

		// cout << "hap.seq: " << _hap.seq << endl;
		// cout << "read: " << r.seq.seq << endl;

	}

	void align()
	{
		if (aligned) return;
		seqan::DnaString _hapSeq(hap_ptr->seq);
		seqan::DnaString _readSeq(read_ptr->seq.seq);
		alignResult=MyAlign();

		seqan::resize(seqan::rows(alignResult), 2);
		seqan::assignSource(seqan::row(alignResult, 0), _hapSeq);
		seqan::assignSource(seqan::row(alignResult, 1), _readSeq);

		stringstream os;
		os << seqan::globalAlignment(alignResult, score) << endl;
		ml.ll=atof(os.str().c_str());


		if (DEBUGCONVERTALIGNMENT) {
			cout << alignResult << endl;
		}


		convertAlignment(alignResult, ml, length(_hapSeq), length(_readSeq), *hap_ptr, *read_ptr);

		reportVariants();
		aligned=true;
	}
	const MLAlignment &  getMLAlignment() { align(); return ml; }

protected:
	void reportVariants()
	{
		const Haplotype & hap = *hap_ptr;
		for (map<int,AlignedVariant>::const_iterator it=hap.indels.begin();it!=hap.indels.end();it++) {
			const AlignedVariant & av=it->second;
			if (av.isCovered(params.padCover, ml.firstBase, ml.lastBase)) ml.hapIndelCovered[it->first]=true; else ml.hapIndelCovered[it->first]=false;
		}
		for (map<int,AlignedVariant>::const_iterator it=hap.snps.begin();it!=hap.snps.end();it++) {
			const AlignedVariant & av=it->second;
			if (av.isCovered(params.padCover, ml.firstBase, ml.lastBase)) ml.hapSNPCovered[it->first]=true; else ml.hapSNPCovered[it->first]=false;
		}
	}

	typedef seqan::Align<seqan::DnaString, seqan::ArrayGaps> MyAlign;
	MyAlign alignResult;
	ObservationModelParameters params;
	const Haplotype *hap_ptr;
	const Read *read_ptr;
	uint32_t hapStart;


	MLAlignment ml;
	seqan::Score<int> score;
	bool aligned;

};

#endif /* OBSERVATIONMODELSEQAN_HPP_ */
