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
 * GetCandidates.cpp
 *
 *  Created on: Aug 27, 2009
 *      Author: caa
 */

#include <fstream>
#include <string>
#include "MyBam.hpp"
#include "foreach.hpp"
#include "ObservationModelSeqAn.hpp"
#include "GetCandidates.hpp"
#include "Variant.hpp"
#include "VariantFile.hpp"
#include "bam.h"
#include "sam.h"
#include "Fasta.hpp"
#include "StringHash.hpp"
#include <set>
#include "foreach.hpp"
using namespace std;

int GetCandidatesFromCIGAR::getIndelFromCIGARFetchFunc(const bam1_t *b, void *data)
{
	CFFData & dat = *( (CFFData *) data);
	vector<CIGARindel>  indels;
	HMap::iterator it;
	getIndelFromCIGAR(b, indels);
	BOOST_FOREACH(CIGARindel id, indels) {
		it = dat.hmap.find(id.refpos);
		if (it==dat.hmap.end()) dat.hmap[id.refpos][id]=1; else (it->second)[id]++;
	}
	return 0;
}

void GetCandidatesFromCIGAR::getIndelFromCIGARRegion(const vector<MyBam *> & myBams, const string & tid, int start, int end, const string & outputFileName, fasta::Fasta & fa)
{
	CFFData data;

	for (size_t b=0;b<myBams.size();b++) {
		bam_fetch(myBams[b]->bf, myBams[b]->idx, myBams[b]->getTID(tid), start, end, &data, &GetCandidatesFromCIGAR::getIndelFromCIGARFetchFunc);
	}

	ofstream ofile(outputFileName.c_str());
	if (!ofile.is_open()) throw string("Cannot open variants file ").append(outputFileName).append(" for writing.");
	outputIndels(tid, data.hmap,ofile,fa,1);
	ofile.close();
}

void GetCandidatesFromCIGAR::getIndelFromCIGAR(const bam1_t *b, vector<CIGARindel> & indels)
{
	const bam1_core_t *c=&b->core;
	uint32_t *cigar=bam1_cigar(b);
	uint32_t k, l=0;
	uint32_t refPos = c->pos;
	int lastop=-1;
	uint32_t lastPos=refPos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int32_t len=cigar[k] >> BAM_CIGAR_SHIFT;
		string seq;

		if (op==BAM_CINS || op==BAM_CMATCH || op==BAM_CSOFT_CLIP) {
			for(int32_t x=0;x<len;x++) {
				if (op==BAM_CINS) {
					seq+=( bam_nt16_rev_table[ bam1_seqi(bam1_seq(b), l) ] );
				}
                //cout << bam_nt16_rev_table[ bam1_seqi(bam1_seq(b), l)];
				l++;
			}
		} else if (op==BAM_CDEL) {
			seq.insert(0, len, 'D');
		}
        //cout << "  cigar #" << k << " " << op << " " << len << " " << refPos << " " << seq << endl;
		if (op==BAM_CINS || op==BAM_CDEL) {
			int ilen=len; if (op==BAM_CDEL) ilen=-ilen;
			indels.push_back(CIGARindel(refPos, ilen, seq));
		}

		// update position for the next cigar
		lastPos=refPos;
		if (op == BAM_CMATCH || op == BAM_CDEL || op==BAM_CREF_SKIP) {
			refPos+=(uint32_t) len;
		} else if (op!=BAM_CINS && op != BAM_CSOFT_CLIP && op != BAM_CHARD_CLIP) throw string("I don't know how to smoke this CIGAR");
		lastop=op;
	}
}

vector<AlignedVariant> GetCandidatesFromCIGAR::alignCIGAR(const string & tid, const CIGARindel & id, fasta::Fasta & fa)
{

	vector<AlignedVariant> variants;

	ObservationModelParameters alignParams("probabilistic");
	seqan::Score<int> score(-1, -460, -100,-960);

	Read rh1;
	rh1.pos=0;
	rh1.posStat.first=0;
	rh1.mapQual=1.0-1e-32;

	map<int, AlignedVariant> alIndel, alSNP;

	int width=params.alignWindow;

	if (abs(id.len)>width/3) width=abs(id.len)*3;

	int start=id.refpos-width;
	int end=id.refpos+width;

	string hap;
	try {
		hap=fa.getSequence(tid, start+1, end+1);
	} catch (string s) {
		cerr << "error: "<< s << endl;
		cerr << "start: " << start << " end: " << endl;
		return vector<AlignedVariant>();
	}
	//int startRef=start-params.refPad;
	//int refEnd=end+params.refPad;

	int startRef=start;
	int refEnd=end;

	string ref;
	try {
		ref=fa.getSequence(tid, startRef+1, refEnd+1);
	} catch (string s) {
		cerr << "error: "<< s << endl;
		cerr << "startRef: " << startRef << " refEnd: " << refEnd<< endl;
		return vector<AlignedVariant>();
	}

	Haplotype hRef; hRef.append(ref);



	// create haplotype with indel


	int pos=id.refpos-start;

	int testlen = (id.len>0)?0:-id.len;
	if (hap.size()<pos+testlen) {
		cerr << "Cannot align variant " << id.refpos << " " << id.len << " " << id.seq << endl;
		return variants;
	}

	if (id.len<0) {
		hap.erase(pos,-id.len);
	} else if (id.len>0) {
		hap.insert(pos, id.seq);
	}

	//cout << "hap: " << hap << endl;

	// align indel

	rh1.seq.seq=hap; // sequence with indel
	rh1.setAllQual(1.0-1e-16);

	try {
		ObservationModelSeqAn om(hRef, rh1, 0, alignParams, score);
		string align;
		om.align();
		const MLAlignment & ml=om.getMLAlignment();
		for(map<int, AlignedVariant>::const_iterator it=ml.indels.begin();it!=ml.indels.end();it++) if (it->second.getType()==AlignedVariant::INS  || it->second.getType()==AlignedVariant::DEL) {
			const AlignedVariant & aid = it->second;
			int pos=startRef+it->first;
			variants.push_back(AlignedVariant(aid.getString(),pos,pos,-1,-1));
		}
	} catch (const bad_alloc & ) {
		cout << "SeqAN Alloc error:  hRef.size(): " << hRef.size() << " rh1.size(): " << rh1.size() << endl;
		// cout << "hRef: " << hRef << endl;
		// cout << "rh1: " << rh1 << endl;
	}


	return variants;

}

void GetCandidatesFromCIGAR::outputIndels(const string & tid, const hash_map<int,map<CIGARindel, int> > & hmap, ofstream & ofile, fasta::Fasta & fa, int outputType=1)
{
   // cout << "outputIndels" << endl;
	hash_map<int,map<CIGARindel, int> >::const_iterator it=hmap.begin();
	hash_map<int,map<AlignedVariant, int> > realigned;


	/*
	ALWAYS realign indel
	if (fastaName.empty()) {
		BOOST_FOREACH(CIGARindel id, indels) {
			it = hmap.find(id.refpos);
			if (it==hmap.end()) hmap[id.refpos][id.seq]=1; else (it->second)[id.seq]++;
			stringstream os;
			os << tid << " " << id.refpos << " " << id.len << " " << id.seq;
			cout << os.str() << endl;
		}
	} else {
	*/
		// realign indel
	for  (it=hmap.begin();it!=hmap.end();it++) {
		for (map<CIGARindel, int>::const_iterator i2=it->second.begin();i2!=it->second.end();i2++) {
			const CIGARindel & id=i2->first;
			//
			//cout << "Here " << tid << " " << id.refpos << " " << id.len << " " << id.seq <<endl;
			vector<AlignedVariant> indels;
			indels=alignCIGAR(tid, id, fa);
			BOOST_FOREACH(AlignedVariant aid, indels)  if (aid.getType()==AlignedVariant::INS  || aid.getType()==AlignedVariant::DEL) {
				realigned[aid.getStartHap()][aid]=i2->second;
			}

		}
	}
		
	std::set<int> positions;
	for  (hash_map<int,map<AlignedVariant, int> >::const_iterator it=realigned.begin();it!=realigned.end();it++) {
		positions.insert(it->first);
	}


	//for  (hash_map<int,map<AlignedVariant, int> >::const_iterator it=realigned.begin();it!=realigned.end();it++) {
	for (std::set<int>::const_iterator posit = positions.begin(); posit != positions.end(); posit++) {
		const map<AlignedVariant, int>  & _variants = realigned[*posit];
		ostringstream ovar, ocnt;
		ovar << tid;
		ovar << " " << *posit;
		for (map<AlignedVariant, int>::const_iterator i2=_variants.begin();i2!=_variants.end();i2++) {
			const AlignedVariant & aid = i2->first;
			int len=aid.size();
			if (aid.getType()==AlignedVariant::DEL) len=-len;
			if (outputType==1) {
				ovar << " " << aid.getString();
				ocnt << " " << i2->second;
			} else if (outputType==2) {
				ovar << " " << len << " " << aid.getSeq();
				ocnt << " " << i2->second;
			} else throw string("Huh?");
		}
		ofile << ovar.str() << " #" << ocnt.str() << endl;

	}

}

void GetCandidatesFromCIGAR::realignCandidateFile(const string & _varFile, bool isOneBased, const string & outputFileName, const string & fastaName)
{
	hash_map<int,map<CIGARindel, int> > hmap;
	hash_map<int,map<CIGARindel, int> >::iterator it;

	fasta::Fasta fa(fastaName);

	VariantFile vf(_varFile);

	ofstream ofile(outputFileName.c_str());
	if (!ofile.is_open()) throw string("Cannot open ").append(outputFileName).append(" for writing CIGAR indels.");

	cout << "Realigning indels from variants file: " << _varFile << endl;

	string ctid="";
	while (!vf.eof()) {
			vector<Variant> variants;
			VariantFile::Candidates cand=vf.getLine(isOneBased);
			if (cand.variants.empty()) continue;

			if (cand.tid!=ctid) {
				if (hmap.size()) {
					outputIndels(ctid, hmap,ofile,fa);
					cout << "Wrote realigned candidate indel for target " << ctid << " to file " << outputFileName << endl;
				}
				hmap.clear();
				ctid=cand.tid;
			}

			BOOST_FOREACH(Variant var, cand.variants) if (var.isIndel()) {
				int len=var.size();
				if (var.getType()==Variant::DEL) len=-len;
				CIGARindel id(cand.pos,len, var.getSeq());
				it = hmap.find(id.refpos);
				if (it==hmap.end()) hmap[id.refpos][id]=1; else (it->second)[id]++;
			}
	}

	outputIndels(ctid, hmap,ofile,fa);
	cout << "Wrote realigned candidate indels for target " << ctid << " to file " << outputFileName << endl;


	ofile.close();
}

void GetCandidatesFromCIGAR::outputLibraries(LibInsertSize & libInsertSize, const string & outputFile)
{

	// open file
	ofstream ofile(outputFile.c_str());
	if (!ofile.is_open()) throw string("Cannot open ").append(outputFile).append(" for writing libraries.");

	for (LibIterator libit = libInsertSize.begin();libit!=libInsertSize.end();libit++) {
		string lib = string(libit->first);
		// compute mean and std
		InsertSizes & insertSizes = libit->second;
		InsIterator insit;

		long int tot = 0;
		double mean = 0.0, std = 0.0;

		std::set<int> isizes;

		for (insit = insertSizes.begin(); insit!=insertSizes.end();insit++) {
			tot += insit->second;
			isizes.insert(insit->first);
		}
		
		double cum = 0;
		int pct = int ( 0.9999 * double(tot));
		int median = tot/2;
		int max_isize = -1;
		int median_isize = -1;
		for (std::set<int>::const_iterator it = isizes.begin();it!=isizes.end();it++) {
			cum += insertSizes[*it];
			if (median_isize == -1 && cum>median) {
				median_isize = *it;
			}
		}
		isizes.clear();
		max_isize = median_isize * 10;
		//cout << "tot: " << tot << " pct: " << pct << " cum: " << cum << " max_isize: " << max_isize <<  " median: " << median <<  " median_isize: " << median_isize << endl;


		double dtot = double(tot);
		for (insit = insertSizes.begin(); insit!=insertSizes.end();insit++) if (insit->first<max_isize) {
//			cout << "isize: " << insit->first << " count: " << insit->second << endl;
			mean += double(insit->first)*double(insit->second)/dtot;
		}
		for (insit = insertSizes.begin(); insit!=insertSizes.end();insit++) if (insit->first<max_isize) {
			double dist = double(insit->first)-mean;
			std += double(insit->second)/dtot*dist*dist;
		}
		cout << "Library: " << lib << " mean: " << mean << " stddev: " << sqrt(std) << endl;
		// create histogram in vector
		int len = int(mean+5*sqrt(std));
		vector<long int> histo(len,2), inthisto(len,2);

		for (insit = insertSizes.begin(); insit!=insertSizes.end();insit++) {
			int isize = insit->first;
			if (isize<len) {
				histo[isize]=insit->second;
			}
		}

		// smooth histogram out a little
		int L = 5;
		for (int i=0;i<len;i++) {
			int min = i-L; if (min<0) min = 0;
			int max = i+L; if (max>len) max = len;
			int n = 0;
			long int sum = 0;
			for (int j=min;j<max;j++,n++) {
				sum += histo[j];
			}
			inthisto[i] = (sum+1)/(n+1);
		}


		// write histogram to file
		ofile << "#LIB " << lib << endl;
		for (int i=0;i<len;i++) {
			ofile << i << " " << inthisto[i] << endl;
		}
	}
	ofile.close();
}
void GetCandidatesFromCIGAR::get(const string & _bamFile, const string & outputFileName, const string & fastaName)
{
	// also get histogram
	LibInsertSize libInsertSize;

	fasta::Fasta fa(fastaName);

	samfile_t *bf;
	bf=samopen(_bamFile.c_str(), "rb", 0);
	bam1_t *b=bam_init1();

	string outputFileVariants = outputFileName;
	string outputFileLibraries = outputFileName;
	outputFileVariants.append(".variants.txt");
	outputFileLibraries.append(".libraries.txt");




	ofstream ofile(outputFileVariants.c_str());
	if (!ofile.is_open()) throw string("Cannot open ").append(outputFileName).append(" for writing CIGAR indels.");

	hash_map<int,map<CIGARindel, int> > hmap;
	hash_map<int,map<CIGARindel, int> >::iterator it;

	int oldtid=-1;
	string _oldtid="";
	cout << "Parsing indels from CIGAR strings..." << endl;
	long int numread = 0;


	string defaultlib("dindel_default");
	while (samread(bf,b)>=0) {

		int btid = b->core.tid;
		if (btid<0) continue; // unmapped read
		const char *tidptr = bf->header->target_name[(b->core).tid];
		if (!tidptr) continue;
		string tid=string(tidptr);
		if ((b->core).tid!=oldtid) {
			if (oldtid!=-1) {
				outputIndels(_oldtid, hmap,ofile,fa);
				cout << "Wrote indels in CIGARS for target " << _oldtid << " to file " << outputFileName << endl;
			}
			oldtid=(b->core).tid;
			_oldtid=tid;
			hmap.clear();

		}

		vector<CIGARindel> indels;
		getIndelFromCIGAR(b, indels);
        // count indels
		BOOST_FOREACH(CIGARindel id, indels) {
			it = hmap.find(id.refpos);
			if (it==hmap.end()) hmap[id.refpos][id]=1; else (it->second)[id]++;
		}

        
        
		// get insertsize
		//cout << int(b->core.flag & BAM_FPAIRED) << " " << int (b->core.flag & BAM_FPROPER_PAIR) << " tid: " << int(b->core.tid) << " mtid: " <<  b->core.mtid << " fdup: " << int(b->core.flag & BAM_FDUP) << " fqcfail: " <<  int( b->core.flag & qBAM_FQCFAIL) << endl;

		if ((b->core.flag & BAM_FPAIRED) && (b->core.flag & BAM_FPROPER_PAIR) && (b->core.tid == b->core.mtid) && !( (b->core.flag & BAM_FDUP) || (b->core.flag & BAM_FQCFAIL) )) {
			const char *p = bam_get_library((bam_header_t *) bf->header, b);

			string & lib = defaultlib;
			if (p) lib = string(p);

			//cout << "lib: " << lib << endl;

			int isize = abs(b->core.isize);

			LibIterator lit = libInsertSize.find(lib);
			if (lit == libInsertSize.end()) {
				libInsertSize[lib]=hash_map<int, long int>();
				lit = libInsertSize.find(lib);
			}

			InsIterator iit = lit->second.find(isize);
			if (iit == lit->second.end()) {
				(lit->second)[isize]=1;
			} else {
				iit->second += 1;
			}

		}
		numread++;
		if (numread % 1000000==999999) {
			cout << "Number of reads read: " << numread+1 << endl;
		}
	}
	outputIndels(_oldtid, hmap,ofile,fa);
	outputLibraries(libInsertSize, outputFileLibraries);

	cout << "Wrote indels in CIGARS for target " << _oldtid << " to file " << outputFileName << endl;
	cout << "Wrote library insert sizes to " << outputFileLibraries << endl;
	cout << "done!" << endl;

	bam_destroy1(b);

	ofile.close();
}

GetCandidatesFromCIGAR::GetCandidatesFromCIGAR()
{

}


GetCandidatesFromCIGAR::~GetCandidatesFromCIGAR()
{

}
