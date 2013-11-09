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
 * Fast.cpp
 *
 *  Created on: Feb 25, 2009
 *      Author: caa
 */

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
#include <sstream>
#include <algorithm>
#include "bam.h"
#include "Haplotype.hpp"
#include "Faster.hpp"
#include "Utils.hpp"
#include "foreach.hpp"
using namespace std;
const int DEBUGS=0;

ObservationModelS::ObservationModelS(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & _params) : params(_params)
{

	hap_ptr = &_hap;
	read_ptr = &r;
	if (params.maxLengthIndel>(int) hap_ptr->size()) throw string("hapSize error.");
	hlen=(int) hap_ptr->seq.size();
	rlen=(int) read_ptr->size();

	this->hapStart=hapStart;

	likelihoodComputed=false;
	bMidError=true;
	computeBMid();
	setupReadLikelihoods();

}

void ObservationModelS::computeBMid()
{
	const Read & read = *read_ptr;
	const Haplotype & hap = *hap_ptr;
	uint32_t hapEnd=hapStart+hap.size();
	uint32_t mReadStart=uint32_t(read.posStat.first);
	uint32_t readEnd=mReadStart+uint32_t(read.size())-1;
	uint32_t olStart, olEnd;
	int mid;

	bMidError=true;
	if (mReadStart>hapEnd) {
		bMid=0;
	} else if (readEnd<hapStart) {
		bMid=int(read.size())-1;
	} else {
		olStart=(hapStart>mReadStart)?hapStart:mReadStart;
		olEnd=(hapEnd>readEnd)?readEnd:hapEnd;
		mid=(int(olEnd)-int(olStart))/2+int(olStart);
		bMid=mid-int(mReadStart);
		bMidError=false;
	}

	if (bMid<0) { bMid=0; };
	if (bMid>=int(read.size())) { bMid=int(read.size())-1; };

	if (DEBUGS) cout << "bMid: " << bMid << endl;

}


void ObservationModelS::setupReadLikelihoods()
{
	const Read & read = *read_ptr;

	logMatch.resize(read.size());
	logMismatch.resize(read.size());
	cumLogMatch.resize(read.size());
	// initialize with prior
	llMatch=0.0;
	if (params.modelType=="probabilistic") {
		for (size_t r=0;r<read.size();r++) {
			double rq=read.qual[r];
			double pr=rq*(1.0-params.pMut);
			double eq=log(.25+.75*pr);
			double uq=log(.75+1e-10-.75*pr);
			logMatch[r]=eq;
			logMismatch[r]=uq;
			llMatch+=eq;
			cumLogMatch[r]=llMatch;
		}
	} else {
		throw string("Model not implemented.");
	}


	double mq=1.0-read.mapQual;
	if (-10.0*log10(mq)>params.capMapQualFast) {
		mq=pow(10.0,-params.capMapQualFast/10.0);
	}

	pOffFirst=mq;
	pOffFirstHMQ=1e-10;

	llOff=log(pOffFirst)+llMatch+double(rlen)*log(1.0-params.pError);
	llOffHMQ=log(pOffFirstHMQ)+llMatch+double(rlen)*log(1.0-params.pError);

}



void ObservationModelS::AlignHash(const HapHash & hash)
{

	const Read & read = *read_ptr;
	hash_map<int,int> hposFreq; // will keep track of frequencies of relative positions of read wrt haplotype
	hash_map<int,int>::iterator it;

	unsigned int kmer = hash.getKmer();

	size_t x=0, xl=read.size()-kmer;

	unsigned int key=hash.convert(read.seq.seq,x);
	for (;x<xl+1;x++) {
		//const set<int> & hpSet=hash.lookup(read.seq.seq,x);
		const set<int> & hpSet = hash.lookup(key);
		if (DEBUGS) cout << "hash: " << x << " :";
		BOOST_FOREACH(int hp, hpSet) {

			int rpfb=hp-x; // relative position of first base wrt haplotype
			if (DEBUGS) cout << " " << rpfb;
			it=hposFreq.find(rpfb);
			// todo weight according to bMid?
			if (it==hposFreq.end()) hposFreq[rpfb]=1; else it->second++;
		}
		if (DEBUGS) cout << endl;
		if (x!=xl) key = hash.pushBack(key, read.seq.seq[x+kmer]);
	}

	// sort according to frequency
	map<int, set<int> > freqToPos;
	for (it=hposFreq.begin();it!=hposFreq.end();it++) {

		if (DEBUGS) cout << "il : " << it->first << " " << it->second << endl;
		freqToPos[it->second].insert(it->first);
	}
	// do alignment with top 15 frequency hash lookups

	const int maxRelPos=15;

	vector<int> relPos; relPos.reserve(maxRelPos);

	int tot=0;

	for (map<int,set<int> >::reverse_iterator rit=freqToPos.rbegin(); rit!=freqToPos.rend() ;rit++) {
		BOOST_FOREACH(int rp, rit->second) {
			if (tot<maxRelPos) {
				relPos.push_back(rp);
				if (DEBUGS) cout << "rp: " << rp << " freq: " << rit->first << endl;
				tot++;
			} else goto _end;
		}
	}
	_end:

	if (DEBUGS) cout << "done"<<endl;
	// run HMM with sparse set of positions
	SStateHMM(relPos);

}

MLAlignment ObservationModelS::align(const HapHash & hash)
{
	AlignHash(hash);
	likelihoodComputed=true;
	reportVariants();
	return ml;
}
/*
inline void ObservationModelS::doTransition(int cr, int nr, const vector<int> & state, vector<double> & alpha, vector<double> & bt, const vector<double> & tr, const int & S)
{
	int r=cr;
	if (state[cr]==-1) {
		// current readbase not fixed
		if (state[nr]==-1) {
			// next base is not fixed
			for (int cs=0;cs<S;cs++) {
				for (int ns=0;ns<S;ns++) {
					double nv=obs[r*S+cs]+alpha[r*S+cs]+tr[cs*S+ns];
					if (nv>alpha[cr*S+ns]+EPS) { alpha[cr*S+ns]=nv; bt[cr*S+ns]=ns; }
				}
			}

		} else {
			// next base is fixed
			for (int cs=0;cs<S;cs++) {
				ns=state[nr];
				double nv=obs[r*S+cs]+alpha[r*S+cs]+tr[cs*S+ns];
				if (nv>alpha[cr*S+ns]+EPS) { alpha[cr*S+ns]=nv; bt[cr*S+ns]=ns; }
			}
		}
	} else {
		// current readbase is fixed
		if (state[nr]==-1) {
			// next base is not fixed
			int cs=state[r];
			for (int ns=0;ns<S;ns++) {
				double nv=obs[r*S+cs]+alpha[r*S+cs]+tr[cs*S+ns];
				if (nv>alpha[cr*S+ns]+EPS) { alpha[cr*S+ns]=nv; bt[cr*S+ns]=ns; }
			}
		} else {
			// next base is fixed
			int cs=state[r];
			int ns=state[nr];
			double nv=obs[r*S+cs]+alpha[r*S+cs]+tr[cs*S+ns];
			if (nv>alpha[cr*S+ns]+EPS) { alpha[cr*S+ns]=nv; bt[cr*S+ns]=ns; }
		}
	}
}

inline void ObservationModelS::doTransitionNF(int cr, int nr, const vector<int> & state, vector<double> & alpha, vector<double> & bt, const vector<double> & tr, const int & S)
{
	int r=cr;
	// next base is not fixed
	for (int cs=0;cs<S;cs++) {
		for (int ns=0;ns<S;ns++) {
			double nv=obs[r*S+cs]+alpha[r*S+cs]+tr[cs*S+ns];
			if (nv>alpha[cr*S+ns]+EPS) { alpha[cr*S+ns]=nv; bt[cr*S+ns]=ns; }
		}
	}
}
*/

void ObservationModelS::SStateHMM(vector<int> & relPos)
{
	// note that this HMM does not keep track of the last base before the insertion, so after the insertion it may transition not to the next haplotype base
	// also, the length of the insertion must be present as the difference between one of the positions in relPos vector.

	// int p1 and p2 are relative positions of first readbase with respect to the haplotype
	if (DEBUGS) cout << "hlen: " << hlen << " rlen: " << rlen << endl;
	const double EPS=1e-7;
	int readLen=read_ptr->size();


	relPos.push_back(-readLen);
	std::sort(relPos.begin(), relPos.end());

	mapState=vector<int>(readLen,0);

	//if (DEBUGS){ cout << "relPos: "; for (int x=0;x<relPos.size();x++) cout << " " << relPos[x]; cout << endl; }


	int S=relPos.size();
	int T=2*S;            // total number of states per slice

	// note that obs will encode observation potentials only for the non-inserted states
	vector<double> tr(S*S, -1000.0), trI(S*S, -1000.0), alpha(readLen*T,-1000.0), obs(readLen*S,0);

	// NOTE alpha is defined as the message that readbase r sends to its neighbour, where neighbour depends on the readbase and bmid

	vector<int> bt(readLen*T,0); // backtracking matrix for Viterbi

	// setup state array
	// initialize to all undetermined
	vector<int> state(readLen,-1);

	// initialize obs_lik (log-emission-probabilities) for every read-base

	for (int r=0;r<readLen;r++) {
		for (int s=0;s<S;s++) {
			int p1=relPos[s];
			if (p1+r>=0 && p1+r<hlen) {
				obs[r*S+s]=(read_ptr->seq.seq[r]==hap_ptr->seq[p1+r])?logMatch[r]:logMismatch[r];
			} else {
				// this corresponds to LO/RO in ObservationModelFB
				obs[r*S+s]=logMatch[r];
			}
		}

		// obs[r*S+S-1]=logMatch[r]; // assume match if insertion
		if (DEBUGS) { cout << "obs: "; for (int s=0;s<S;s++) cout << " " << -int(round(obs[r*S+s])); cout << endl; }
	}




	// todo : add code to fix state to OffHaplotype if to the left or right of a fixed base?


	// setup transition-matrix

	vector<double> prior(T, -1000.0), priorHMQ(T, -1000.0);

 	// p1 <- p1
	// p1 <- p2
	// p1 <- I

	// p2 <- p1
	// p2 <- p2
	// p2 <- I

	// I <- p1
	// I <- p2
	// I <- I


	// setup prior distribution for bMid
	for (int ins=0;ins<2;ins++) {
		double pins=(ins==0)?log(1.0-params.pError):log(params.pError);
		for (int y=0;y<S;y++) {
			int x=y+ins*S;
			int hp=relPos[y]+bMid;
			if (hp>=0 && hp<hlen) {
				prior[x]=log(1.0-pOffFirst)+pins;
				priorHMQ[x]=log(1.0-pOffFirstHMQ)+pins;
			} else {
				prior[x]=log(pOffFirst)+pins;
				priorHMQ[x]=log(pOffFirstHMQ)+pins;
			}
			if (DEBUGS) cout << "prior[" << x << "]: " << prior[x] << " " << priorHMQ[x] << endl;
		}
	}

	double logpInsgNoIns = log(params.pError);
	double logpInsgIns = -0.25;
	double logpNoInsgIns = log(1-exp(logpInsgIns));
	//double logpNoInsgNoIns = log(1.0-params.pError);



	// transitions between relPos
	for (int s1=0;s1<S;s1++) for (int s2=0;s2<S;s2++) {
		double ll=-1000.0;
		// relpos to relpos
		// for non-inserted states only deletions are allowed.
		// you can only transition to a lower relPos from an insertion-state (ie x>=S)
		if (s1!=s2) {
			double d=fabs(double(relPos[s1]-relPos[s2]));
			ll=(d-1.0)*logpInsgIns+log(params.pError);
			trI[s1*S+s2]=(d-1.0)*logpInsgIns;
		} else if (s1==s2) {
			ll=log(1.0-params.pError);
		}

		// Pr[s1 | s2 ]
		tr[s1*S+s2]=ll;
	}

	if (DEBUGS) for (int s1=0;s1<S;s1++) {
		cout << "tr["<< s1 << "]: "; for (int s2=0;s2<S;s2++) cout << " " << tr[s1*S+s2]; cout << endl;
	}
	// from left to bMid

	for (int r=0;r<bMid;r++) {
		int cr=r;
		//doTransition(cr, nr, state, alpha, bt, tr);

		for (int cs=0;cs<S;cs++) {
			double pv=obs[r*S+cs]; if (r) pv+=alpha[(r-1)*T+cs];

			// transition to non-inserted from non-inserted
			for (int ns=cs;ns<S;ns++) {
				double nv=pv+tr[cs*S+ns];
				if (nv>alpha[cr*T+ns]+EPS) { alpha[cr*T+ns]=nv; bt[cr*T+ns]=cs; }
			}

			// r          <---  r+1
			// to non-ins from ins
			int ns=cs+S;
			double nv=pv+logpNoInsgIns;
			if (nv>alpha[cr*T+ns]+EPS) { alpha[cr*T+ns]=nv; bt[cr*T+ns]=cs; }


			// insertion states

			// r          <---  r+1
			// ins        <---  ins

			int ics=cs+S;
			ns=ics;
			nv=logMatch[r]+logpInsgIns; if (r) nv += alpha[(r-1)*T+ics];
			if (nv>alpha[cr*T+ns]+EPS) { alpha[cr*T+ns]=nv; bt[cr*T+ns]=ics; }


			// ins       <---   noins
			ics=cs+S;		  //  must transition to a lower relPos in case of insertion and going from left to right
			for (int ns=0;ns<cs;ns++) if (relPos[cs]-r>=relPos[ns]) {
				nv=logMatch[r]+trI[cs*S+ns]+logpInsgNoIns; if (r) nv += alpha[(r-1)*T+ics];
				if (nv>alpha[cr*T+ns]+EPS) { alpha[cr*T+ns]=nv; bt[cr*T+ns]=ics; }
			}


		}

		if (DEBUGS) { cout << "alpha_fw: "; for (int x=0;x<T;x++) cout << " " << alpha[r*S+x]; cout << endl; }
	}

	if (DEBUGS) cout << endl;

	// from right to bMid

	for (int r=readLen-1;r>bMid;r--) {
		int cr=r;
		//doTransition(cr, nr, state, alpha, bt, tr);

		for (int cs=0;cs<S;cs++) {
			double pv=obs[r*S+cs]; if (r<readLen-1) pv+=alpha[(r+1)*T+cs];

			// transition to non-inserted from non-inserted
			for (int ns=0;ns<=cs;ns++) {
				double nv=pv+tr[cs*S+ns];
				if (nv>alpha[cr*T+ns]+EPS) { alpha[cr*T+ns]=nv; bt[cr*T+ns]=cs; }
			}

			// r          <---  r-1
			// to ins     from no-ins
			double nv=logMatch[r]+logpInsgNoIns; if (r<readLen-1) nv += alpha[(r+1)*T+cs+S];

			if (nv>alpha[cr*T+cs]+EPS) { alpha[cr*T+cs]=nv; bt[cr*T+cs]=cs+S; }

			int ns;

			// insertion states

			// r          <---  r-1
			// ins        <---  ins

			int ics=cs+S;
			ns=ics;
			nv=logMatch[r]+logpInsgIns; if (r<readLen-1) nv+= alpha[(r+1)*T+ics];
			if (nv>alpha[cr*T+ns]+EPS) { alpha[cr*T+ns]=nv; bt[cr*T+ns]=ics; }

			// r	       <---   r-1
			// noins       <---   ins
			ics=cs+S;		  //  must transition to a lower relPos in case of insertion and going from left to right
			for (int ns=cs+1;ns<S;ns++) if (relPos[cs]>relPos[ns]-r) {
				nv=obs[r*S+cs]+logpNoInsgIns+trI[cs*S+ns]; if (r<readLen-1) nv +=  alpha[(r+1)*T+cs];
				if (nv>alpha[cr*T+ns+S]+EPS) { alpha[cr*T+ns+S]=nv; bt[cr*T+ns+S]=cs; }
			}


		}
			// r               r-1
			// ins     <-----  noins

		if (DEBUGS) { cout << "alpha_bw: "; for (int x=0;x<T;x++) cout << " " << alpha[r*T+x]; cout << endl; }
	}


	double max=-HUGE_VAL;
	int xmax=0;

	for (int ins=0;ins<2;ins++)
	for (int y=0;y<S;y++) {
		int x=ins*S+y;
		double obsv=(ins==0)?obs[bMid*S+y]:logMatch[bMid];
		alpha[bMid*T+x]=obsv+prior[x];
		if (bMid<readLen-1) alpha[bMid*T+x]+=alpha[(bMid+1)*T+x];
		if (bMid>0) alpha[bMid*T+x]+=alpha[(bMid-1)*T+x];

		if (alpha[bMid*T+x]>max) {
			max=alpha[bMid*T+x];
			xmax=x;
		}
	}

	if (DEBUGS) { cout << "alpha_bmid: "; for (int x=0;x<T;x++) cout << " " << alpha[bMid*T+x]; cout << endl; }

	// check position of bMid on haplotype

	int hp=relPos[xmax%S]+bMid;
	if (hp>=0 || hp < hlen) {
		// bMid is an insertion

		ml.offHap=false;

	} else {
		// not an insertion
		// is it on or off the haplotype?

		ml.offHap=true;


	}


	ml.ll=max;

	max=-HUGE_VAL;
	xmax=0;

	if (DEBUGS) cout << "alpha_bmid_HMQ: ";
	for (int ins=0;ins<2;ins++)
	for (int y=0;y<S;y++) {
		int x=ins*S+y;
		double obsv=(ins==0)?obs[bMid*S+x]:logMatch[bMid];
		double v=obsv+priorHMQ[x];
		if (bMid<readLen-1) v+=alpha[(bMid+1)*T+x];
		if (bMid>0) v+=alpha[(bMid-1)*T+x];

		if (v>max) {
			max=v;
			xmax=x;
		}
		if (DEBUGS) cout << " " << v;
	}
	if (DEBUGS) cout << endl;

	hp=relPos[xmax%S]+bMid;
	if (hp>=0 || hp < hlen) {
		// bMid is an insertion
		ml.offHapHMQ=false;

	} else {
		// not an insertion
		// is it on or off the haplotype?
		ml.offHapHMQ=true;
	}

	state[bMid]=xmax;

	// backtrack to get the map state

	for (int b=bMid; b>0;b--) {
		state[b-1]=bt[(b-1)*T+state[b]];
	}

	for (int b=bMid;b<readLen-1;b++) {
		state[b+1]=bt[(b+1)*T+state[b]];
	}

	if (DEBUGS){ cout << "state: "; for (int r=0;r<readLen;r++) cout << "[" << r << " " << read_ptr->seq.seq[r] << " " << state[r] << "]"; cout << endl;}


	// convert relative positions to absolute positions, using LO, RO, x convention

	int lhp=1;
	for (int r=0; r<readLen; r++) {
		if (state[r]==-1) throw string("error in mapstate fast");
		if (state[r]<S) {
			int hp=relPos[state[r]]+r;
			if (hp>=0 && hp<hlen) {
				mapState[r]=hp+1;
				lhp=hp+1;
			} else if (hp<0) mapState[r]=0; else mapState[r]=hlen; // LO and RO
			if (DEBUGS) cout << "ms: " << r << " " << state[r] << " " << hp << endl;

		} else {
			// insertion
			mapState[r]=hlen+2+lhp;

		}
		if (DEBUGS) cout << "ms: " << r << " " << state[r] << " mapstate " << mapState[r] << endl;
	}


}


void ObservationModelS::reportVariants()
{
	int hapSize=hlen;
	int readSize=rlen;
	int numS=hapSize+2;

	const Read & read = *read_ptr;
	const Haplotype & hap = *hap_ptr;


	ml.align=string(hapSize, 'R');
	ml.indels.clear();
	ml.snps.clear();

	ml.firstBase=-1;
	ml.lastBase=-1;
	ml.hapIndelCovered.clear();
	ml.hapSNPCovered.clear();
	ml.hpos.clear();
	ml.hpos.resize(readSize);


	int b=0;
	while (b<readSize) {
		// only report variants for bases that are on the haplotype
		int s=mapState[b];
		if ( (s%numS)>0 && (s%numS)<=hapSize ) {
			if (s>=numS) { // insertion
				int pos=(s%numS)-1+1; // position of insertion wrt haplotype MAINTAIN CONVENTION OF INSERTION BEFORE BASE X
				int len=0; // length of insertion
				int rpos=b; // start base of insertion in read
				while (b<readSize && mapState[b]>=numS) {
					ml.hpos[b]=MLAlignment::INS;
					b++;
					len++;
				}
				int readStart=rpos;
				int readEnd=b-1;
				int hapStart=pos;
				int hapEnd=pos;
				string seq=read.seq.seq.substr(rpos,len);
				ml.indels[pos]=AlignedVariant(string("+").append(seq),  hapStart, hapEnd, readStart, readEnd);
				b--;
			} else {
				ml.hpos[b]=s-1;
				// update firstBase and lastBase
				if (ml.firstBase==-1) ml.firstBase=s-1; else if (s-1<ml.firstBase) ml.firstBase=s-1;
				if (ml.lastBase==-1) ml.lastBase=s-1; else if (s-1>ml.lastBase) ml.lastBase=s-1;


				// check for SNP
				if (read.seq[b]!=hap.seq[s-1]) {
					string snp;
					snp+=hap.seq[s-1];
					snp.append("=>");
					snp+=read.seq[b];
					int readStart=b;
					int readEnd=b;
					int hapStart=s-1;
					int hapEnd=s-1;


					ml.snps[s-1]=AlignedVariant(snp,hapStart, hapEnd, readStart, readEnd);
					ml.align[s-1]=read.seq[b];
				}
				// check for deletion
				if (b<readSize-1) {
					int ns=mapState[b+1];
					if (ns<numS && ns-s>1) { // make sure next state is not an insertion..
						int pos=s+1-1;
						int len=-(ns-s-1);
						//indels[pos]=ReportVariant(len, hap.seq.substr(pos, -len), b);

						for (int y=pos;y<-len+pos;y++) ml.align[y]='D';
						int readStart=b;
						int readEnd=b+1;
						int hapStart=pos;
						int hapEnd=pos-len-1;
						string seq=hap.seq.substr(pos,-len);
						ml.indels[pos]=AlignedVariant(string("-").append(seq), hapStart, hapEnd, readStart, readEnd);
					}
				}

			}

		} else {// on haplotype
			if (s%numS==0) ml.hpos[b]=MLAlignment::LO; else ml.hpos[b]=MLAlignment::RO;

		}
		b++;
	}

	for (map<int,AlignedVariant>::const_iterator it=hap.indels.begin();it!=hap.indels.end();it++) {
		const AlignedVariant & av=it->second;
		if (av.isCovered(params.padCover, ml.firstBase, ml.lastBase)) ml.hapIndelCovered[it->first]=true; else ml.hapIndelCovered[it->first]=false;
	}
	for (map<int,AlignedVariant>::const_iterator it=hap.snps.begin();it!=hap.snps.end();it++) {
		const AlignedVariant & av=it->second;
		if (av.isCovered(params.padCover, ml.firstBase, ml.lastBase)) ml.hapSNPCovered[it->first]=true; else ml.hapSNPCovered[it->first]=false;
	}


}

void ObservationModelS::printAlignment(size_t hapScrPos)
{
	// count how many bases in the read are left of the haplotype
	if (!likelihoodComputed) throw string("Must align() first!");
	int hapSize=hlen;
	int readSize=rlen;
	int numS=hapSize+2;

	const Read & read = *read_ptr;
	const Haplotype & hap = *hap_ptr;


	string leftHap, rightHap;
	string rhap(hap.size(),' ');
	string ins;

	bool insact=false;
	int b=0;
	while (b<readSize) {
		// only report variants for bases that are on the haplotype
		int s=mapState[b];
		char nuc=read.seq.seq[b];
		if (s%numS==0) {
			//
			leftHap+=nuc;
		} else if ( (s%numS)>0 && (s%numS)<=hapSize ) {
			if (s>=numS) { // insertion
				if (!insact) {
					insact=true;
					ins+='[';
					stringstream os; os << (s%numS);
					ins.append(os.str());
					ins+=' ';
				}

				ins+=nuc;

			} else {
				if (insact) ins+=']';
				insact=false;
				rhap[s-1]=nuc;

				if (b<readSize-1) {
					int ns=mapState[b+1];
					if (ns<numS && ns-s>1) {
						int len=ns-s-1;
						rhap.replace(s, len, string(len,'_'));
					}


				}


			}

		} else {
			rightHap+=nuc;
		}
		b++;
	}
	if (insact) ins+=']';

	stringstream os;
	os << readSize << " " << ml.offHap << " " << ml.indels.size() << " " << ml.firstBase << " " << ml.lastBase << " " << ml.ll << " ";
	for (map<int,AlignedVariant>::const_iterator it=hap.indels.begin();it!=hap.indels.end();it++) {
		if  (ml.hapIndelCovered[it->first]) os << "1 "; else os << "0 ";
	}
	string prefix=os.str();

	int leftHapSpace=int(hapScrPos)-int(prefix.size());
	if (leftHapSpace<0) leftHapSpace=0;

	string prLeftHap=string(leftHapSpace,' ');

	if (int(leftHap.size())>leftHapSpace) {
		prLeftHap=leftHap.substr(leftHap.size()-leftHapSpace, leftHapSpace);
	} else if (leftHap.size()>0) {
		prLeftHap.replace(leftHapSpace-leftHap.size(), leftHap.size(), leftHap);
	}

	cout << prefix<<prLeftHap<<rhap<<rightHap << " " << ins << " read: " << read.seq.seq << endl;


	for (map<int,AlignedVariant>::const_iterator it=hap.indels.begin();it!=hap.indels.end();it++) {
		cout << " " << it->first;
	}
	cout << endl;

	cout << endl;

	for (int x=0;x<readSize;x++) {
		cout << "[" << x << ":" << ml.hpos[x] << "]";
	}
	cout << endl;
}

ObservationModelS::~ObservationModelS()
{


}


