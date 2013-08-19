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
#include <vector>
#include <cmath>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "ObservationModelFB.hpp"
#include "Haplotype.hpp"
#include "Read.hpp"
#include "ReadIndelErrorModel.hpp"
using namespace std;
#define DEBUGHMM

ObservationModelFB::ObservationModelFB(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & _params) : read(r), params(_params)
{

	Init(_hap, hapStart);

}

void ObservationModelFB::Init(const Haplotype & _hap, uint32_t hapStart)
{
	hap=_hap.filtered();
	memAllocated=false;
	HMMInitialized=false;
	HMMConsistent=false;
	obsInitialized=false;
	forwardDone=false;
	backwardDone=false;
	likelihoodComputed=false;
	makeObsVector=false;
	marginalsComputed=false;
	if (params.maxLengthDel>(int) hap.size()) throw string("hapSize error.");


	//bMid=_bMid; compute bMid position
	uint32_t hapEnd=hapStart+hap.size();
	uint32_t mReadStart=uint32_t(read.posStat.first);
	uint32_t readEnd=mReadStart+uint32_t(read.size())-1;
	uint32_t olStart, olEnd;
	int mid;
	if (read.isUnmapped()) {
		bMid = int ( read.size()/2 );
        } else {
		if (mReadStart>hapEnd) {
			bMid=int(read.size()/2);
			/*
			cout << "hapStart: " << hapStart << " hapEnd: " << hapEnd << endl;
			cout << "mReadStart: " << mReadStart << " readEnd: " << readEnd << endl;
			cout << "read.posStat.first: " << read.posStat.first << " read.posStat.second: " << read.posStat.second << endl;
			cout << "Read: " << read << endl;
			cout << "BMID error: read is not on haplotype. Changing ObservationModelParameters." << endl;
			*/
			//cerr << "BMIDE" << endl;
			params.baseQualThreshold=0.0;
		} else if (readEnd<hapStart) {
			bMid=int(read.size()/2);
			params.baseQualThreshold=0.0;
			/*
			cout << "hapStart: " << hapStart << " hapEnd: " << hapEnd << endl;
			cout << "mReadStart: " << mReadStart << " readEnd: " << readEnd << endl;
			cout << "read.posStat.first: " << read.posStat.first << " read.posStat.second: " << read.posStat.second << endl;
			cout << "Read: " << read << endl;
			cout << "BMID error: read is not on haplotype. Changing ObservationModelParameters" << endl;
			*/
			//cerr << "BMIDE" << endl;
		} else {
			olStart=(hapStart>mReadStart)?hapStart:mReadStart;
			olEnd=(hapEnd>readEnd)?readEnd:hapEnd;
			mid=(int(olEnd)-int(olStart))/2+int(olStart);
			bMid=mid-int(mReadStart);
		}
	}

	/*
	if (bMid<0 || bMid>=int(read.size())) {
		cout << "hapStart: " << hapStart << " readStart: " << mReadStart << " readEnd: " << readEnd << " olStart: " << olStart << " olEnd: " << olEnd << " bMid: " << bMid << " mid: " << mid << endl;
		throw string("error");
	}
	*/

	if (params.bMid!=-1) bMid=params.bMid;

	if (bMid<0) { cout << "BMIDERROR" << endl; bMid=0; };
	if (bMid>=int(read.size())) { cout << "BMIDERROR" << endl; bMid=int(read.size())-1; };

	this->hapStart=hapStart;
}


void ObservationModelFB::changeHaplotype(const Haplotype & newHap)
{
	if (hap.size()!=newHap.filtered().size()) {
		cout << "hap: " << hap << " newHap: " << newHap << endl;
		throw string("New haplotype must have same length as old haplotype.");
	}

	hap=newHap.filtered();

	obsInitialized=false;
	HMMConsistent=false;
	forwardDone=false;
	backwardDone=false;
	likelihoodComputed=false;

}

void ObservationModelFB::calcLikelihoodFromLastSlice()
{
	throw string("CHANGE ME! PRIOR NOT CALCULATED IN RIGHT PLACE");
	if (likelihoodComputed) return;
	double *alpha_l=alpha[bMid];
	double *beta_l=beta[bMid];
	double *obs_l=obs[bMid];
	logLikelihood=0.0;
	logLikelihoodNoPrior=0.0;
	likOffHap.resize(2);
	likOffHap[0]=0.0;
	likOffHap[1]=0.0;

	int y=0;
	double max=0.0;
	int maxidx=0;
	for (int x=0;x<2*numS;x++, y++) {

		double v=alpha_l[y]+obs_l[y]+beta_l[y];
		double mar=exp(v);
		if (mar>max) {
			max=mar;
			maxidx=x;
		}
		logLikelihood+=mar;
		v=alpha_l[y]+obs_l[y]+beta_l[y]-prior[x];
		double marnp=exp(v);
		logLikelihoodNoPrior+=marnp;
		if ((x%numS)==0) likOffHap[0]+=marnp; else if ((x%numS)!=ROState) likOffHap[1]+=marnp;
	}
	logLikelihood=log(logLikelihood);
	logLikelihoodNoPrior=log(logLikelihoodNoPrior);
	likOffHap[0]=log(likOffHap[0]);
	likOffHap[1]=log(likOffHap[1]);

	ml.ll=logLikelihood;
	ml.llOff=likOffHap[0];
	ml.llOn=likOffHap[1];
	if ((maxidx%numS)==0 || (maxidx%numS)==ROState) {
		ml.offHapHMQ=true;
		ml.offHap=true;
	}
#ifdef DEBUGHMM
	//cout << "calcLikelihoodFromLastSlice(): " << logLikelihood << endl;
	//cout << "here: " << scientific << setprecision(10) <<  log(likOffHap[0]*exp(prior[0])+likOffHap[1]*exp(prior[1])) << " " << logLikelihood << endl;
#endif
	likelihoodComputed=true;
}

MLAlignment ObservationModelFB::calcLikelihood()
{

	initHMM();
	setupReadObservationPotentials();
	computeForwardMessages();
	calcLikelihoodFromLastSlice();

	return ml;
}


void ObservationModelFB::setupTransitionProbs()
{
	logpLOgLO=log(1.0-params.pFirstgLO);
	logpFirstgLO=log(params.pFirstgLO);

	numT=params.maxLengthDel+2;
	logPTrans.resize(numT);
	// maxT is the transition which corresponds to a normal-operation base extension
	logPTrans[1]=log(1.0-params.pError);
	double norm=0.0;
	for (int x=1;x<numT;x++) if (x!=1) {
		double p=-fabs(1.0-double(x));
		logPTrans[x]=p;
		norm+=exp(p);
	}
	norm=log(norm/params.pError);
	for (int x=1;x<numT;x++) if (x!=1) logPTrans[x]-=norm;

	// check norm
	norm=0.0;
	for (int x=1;x<numT;x++) norm+=exp(logPTrans[x]);
	assert(fabs(norm-1.0)<1e-15);

	logpInsgIns=-1.0;
	logpNoInsgIns=log(1.0-exp(logpInsgIns));
	logpInsgNoIns=log(params.pError);
	logpNoInsgNoIns=log(1-params.pError);
	/*
	cout << "logpInsgIns: " << logpInsgIns << endl;
	cout << "logpNoInsgIns: " << logpNoInsgIns << endl;
	cout << "logpInsgNoIns: " << logpInsgNoIns << endl;
	cout << "logpNoInsgNoIns: " << logpNoInsgNoIns << endl;
	*/

}


void ObservationModelFB::setupReadObservationPotentials()
{
	if (obsInitialized) return;
	int b;
	if (params.modelType=="probabilistic") {
		for (b=0;b<readSize;b++) {
			double rq=read.qual[b];
			char nuc=read.seq[b];
			double *obs_b=obs[b];

			double *obs_b_ins=&obs_b[numS];
			double *obs_b_noins=obs_b;
			double pr=rq*(1.0-params.pMut);
			double eq=log(.25+.75*pr);
			double uq=log(.75+1e-10-.75*pr);

/* TODO: re-think how to deal with unmapped reads
			obs_b_ins[0]=eq; // left of haplotype
			obs_b_ins[hapSize+1]=eq; // right of haplotype

			obs_b_noins[0]=eq; // left of haplotype
			obs_b_noins[hapSize+1]=eq; // right of haplotype
*/
            double penalty = -0.1;
            obs_b_ins[0]=eq + penalty; // left of haplotype
			obs_b_ins[hapSize+1]=eq + penalty; // right of haplotype

			obs_b_noins[0]=eq + penalty; // left of haplotype
			obs_b_noins[hapSize+1]=eq + penalty; // right of haplotype


			for (int y=0;y<hapSize;y++) {
				// given an insertion in the read assume match to prevent favoring of the insertion based on low base qualities
				obs_b_ins[y+1]=eq;
				if (hap[y] == 'N' || hap[y]==nuc) {
					obs_b_noins[y+1]=eq;
				} else {
					obs_b_noins[y+1]=uq;
				}
			}
		}
	} else throw string("Unsupported observation model.");





	//throw string("Check priors!");
	//cout << "prior[0]: " << prior[0] << endl;
	if (params.forceReadOnHaplotype) {
		forceOnHap();
	}

	obsInitialized=true;
}

void ObservationModelFB::computeBMidPrior(vector<double> & _prior, double mapQual)
{
	double pOffFirst;
	double mq=1.0-mapQual;
	if (-10.0*log10(mq)>params.mapQualThreshold) {
		mq=pow(10.0,-params.mapQualThreshold/10.0);
	}
	pOffFirst=mq;

	_prior=vector<double>(2*numS,0.0);
	vector<double> pinsert = vector<double>(numS,0.0);
	if (params.mapUnmappedReads && read.isPaired()) {
		//cout << "read: " << bam1_qname(read.getBam()) << " read.pos: " << read.pos << " read.matePos: " << read.matePos << " read.mateLen: " << read.mateLen <<  " read.isUnmapped: " << read.isUnmapped() << " mateUnmapped: " << read.mateIsUnmapped() << " bmid: " << bMid << " hapStart: " << hapStart << " read.tid: " << read.getBam()->core.tid << " read.mtid: " << read.getBam()->core.mtid << endl;
		//

		if (!read.mateIsUnmapped() && read.mateLen != -1 && read.getBam()->core.tid == read.getBam()->core.mtid) {
			if (read.mateIsReverse()) {
				for ( int x=1;x<hapSize+1;x++) pinsert[x] = log(read.getLibrary().getProb(abs(hapStart+x-bMid-int(read.matePos+read.mateLen))));
			} else {
				for ( int x=1;x<hapSize+1;x++) pinsert[x] = log(read.getLibrary().getProb(abs(hapStart+x+readSize-bMid-int(read.matePos))));
			}
			pinsert[0] = log(read.getLibrary().getNinetyFifthPctProb());
			// for (int x=0;x<hapSize+1;x++) cout << " " << pinsert[x];
			// cout << endl;
		}

	}

	for (size_t i=0;i<2;i++) {
		double logpIns=(i==1)?(logpInsgNoIns):log(1.0-exp(logpInsgNoIns));
		_prior[i*numS+0]=log(pOffFirst)+logpIns+pinsert[0];
		_prior[i*numS+ROState]=-100.0;
		for (int x=1;x<hapSize+1;x++) {
			_prior[i*numS+x]=pinsert[x]+log((1.0-pOffFirst))+logpIns;
		}
	}

}

void ObservationModelFB::forceOnHap()
{
	for (int b=0;b<readSize;b++) {
		obs[b][0]=-1000.0;
		obs[b][ROState]=-1000.0;
		obs[b][numS]=-1000.0;
		obs[b][ROState+numS]=-1000.0;
	}

}

void ObservationModelFB::getObsVector(int b, double *vec) const
{
	throw string("Not implemented");
//	for (int y=0;y<4;y++) vec[y]=obsVector[(b<<2)+y];
}

void ObservationModelFB::initHMM()
{
	if (HMMInitialized) return;
	hapSize=hap.size();
	numS=hapSize+2;
	readSize=read.seq.size();
	ROState=hapSize+1;


	allocateMemory();

	for (int x=0;x<2*numS;x++) {
			alpha[0][x]=0.0;
			beta[readSize-1][x]=0.0;
	}

	HMMInitialized=true;
	HMMConsistent=false;
	forwardDone=false;
	backwardDone=false;
	likelihoodComputed=false;
	marginalsComputed=false;

	setupTransitionProbs();

}

void ObservationModelFB::allocateMemory()
{
	if (memAllocated) return; //throw string("Memory already allocated.");
	mar.reserve(readSize);
	obs.reserve(readSize);
	alpha.reserve(readSize);
	beta.reserve(readSize);

	for (int b=0;b<readSize;b++) {
		mar.push_back(new double[numS*2]);
		obs.push_back(new double[numS*2]);
		alpha.push_back(new double[numS*2]);
		beta.push_back(new double[numS*2]);
		xmar.push_back(new double[numS]);
	}
	if (makeObsVector) { obsVector=new double[4*readSize]; };
	memAllocated=true;
}

void ObservationModelFB::deleteMemory()
{
	if (memAllocated) {
		for (int b=0;b<readSize;b++) {
			delete[] mar[b];
			delete[] obs[b];
			delete[] alpha[b];
			delete[] beta[b];
			delete[] xmar[b];
		}
		if (makeObsVector) delete[] obsVector;
		memAllocated=false;
	}

}



void ObservationModelFB::passMessageOneInc(double *alpha_l, const double *alpha_l_1, const double *obs_l_1)
{
	// P(x^l=x^{l-1}+d|x^{l-1})=(10^{logpSkip*d})
	for (int x=0;x<2*numS;x++) alpha_l[x]=0.0;

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	alpha_l[0]+=( exp(obs_l_1[0]+alpha_l_1[0]+logpLOgLO+logpNoInsgNoIns) );
	alpha_l[1]+=( exp(obs_l_1[0]+alpha_l_1[0]+logpFirstgLO+logpNoInsgNoIns ) );
	for (int x=1;x<=hapSize;x++ ) {
		double tmp=obs_l_1[x]+alpha_l_1[x]+logpNoInsgNoIns;
		for (int y=1;y<numT;y++) {
			int newx=x+y;
			if (newx>hapSize) newx=ROState;
			alpha_l[newx]+=exp(tmp+logPTrans[y]);
		}
	}
	// RO -> RO pROgRO=1.0;
	alpha_l[ROState]+=exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpNoInsgNoIns);

	// x^l, i^l=0 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		alpha_l[numS+x]+=exp(obs_l_1[x]+alpha_l_1[x]+logpInsgNoIns);
	}

	// x^l, i^l=1 = > x^{l+1}, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		alpha_l[numS+x]+=exp(obs_l_1[numS+x]+alpha_l_1[numS+x]+logpInsgIns);
	}

	// x^l, i^l=1 = > x^{l+1}=x^l+1, i^{l+1}=0
	alpha_l[0]+=exp(obs_l_1[numS+0]+alpha_l_1[numS+0]+logpNoInsgIns); // cannot go from insertion on to the haplotype
	for (int x=1;x<=hapSize+1;x++ ) {
		int newx=x+1; if (newx>ROState) newx=ROState;
		alpha_l[newx]+=exp(obs_l_1[numS+x]+alpha_l_1[numS+x]+logpNoInsgIns);
	}


	// convert back to log

	for (int x=0;x<2*numS;x++) alpha_l[x]=log(alpha_l[x]);
}

// doing a pass for a P(X_{l+1}|X_l) potential where the next state is lower than current state

void ObservationModelFB::passMessageOneDec(double *alpha_l, const double *alpha_l_1,const double *obs_l_1)
{
	// P(x^l=x^{l-1}+d|x^{l-1})=(10^{logpSkip*d})
	for (int x=0;x<2*numS;x++) alpha_l[x]=0.0;

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	alpha_l[ROState]+=( exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpLOgLO+logpNoInsgNoIns) );
	alpha_l[hapSize]+=( exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpFirstgLO+logpNoInsgNoIns ) );
	for (int x=1;x<=hapSize;x++ ) {
		double tmp=obs_l_1[x]+alpha_l_1[x]+logpNoInsgNoIns;
		for (int y=1;y<numT;y++) {
			int newx=x-y;
			if (newx<0) newx=0;
			alpha_l[newx]+=exp(tmp+logPTrans[y]);
		}
	}
	// RO -> RO pROgRO=1.0;
	alpha_l[0]+=exp(obs_l_1[0]+alpha_l_1[0]+logpNoInsgNoIns);

	// x^l, i^l=0 = > x^{l+1}=x^l-1, i^{l+1}=1
	alpha_l[numS+ROState]+=( exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpLOgLO+logpInsgNoIns ) );
	alpha_l[numS+hapSize]+=( exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpFirstgLO+logpInsgNoIns ) );
	for (int x=0;x<=hapSize;x++ ) {
		int newx=x-1; if (newx<0) newx=0;
		alpha_l[numS+newx]+=exp(obs_l_1[x]+alpha_l_1[x]+logpInsgNoIns);
	}


	/*
	for (int x=0;x<=hapSize+1;x++ ) {
		alpha_l[numS+x]+=exp(obs_l_1[x]+alpha_l_1[x]+logpInsgNoIns);
	}
	*/

	// x^l, i^l=1 = > x^{l+1}, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		alpha_l[numS+x]+=exp(obs_l_1[numS+x]+alpha_l_1[numS+x]+logpInsgIns);
	}



	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=0
	for (int x=0;x<=hapSize+1;x++ ) {
		alpha_l[x]+=exp(obs_l_1[numS+x]+alpha_l_1[numS+x]+logpNoInsgIns);
	}

	/*
	alpha_l[ROState]+=exp(obs_l_1[numS+ROState]+alpha_l_1[numS+ROState]+logpNoInsgIns); // cannot go from insertion on to the haplotype
	for (int x=0;x<=hapSize;x++ ) {
		int newx=x-1; if (newx<0) newx=0;
		alpha_l[newx]+=exp(obs_l_1[numS+x]+alpha_l_1[numS+x]+logpNoInsgIns);
	}
	*/


	// convert back to log
	for (int x=0;x<2*numS;x++) alpha_l[x]=log(alpha_l[x]);
}

void ObservationModelFB::passMessageTwoInc(double *beta_l, const double *beta_l_1,const double *obs_l_1)
{

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	beta_l[0]=( exp(obs_l_1[0]+beta_l_1[0]+logpLOgLO+logpNoInsgNoIns) ) + ( exp(obs_l_1[1]+beta_l_1[1]+logpFirstgLO+logpNoInsgNoIns ) );
	for (int x=1;x<=hapSize;x++ ) {
		// double tmp=beta_l_1[x]+logpNoInsgNoIns;
		beta_l[x]=0.0;
		for (int y=1;y<numT;y++) {
			int newx=x+y;
			if (newx>hapSize) newx=ROState;
			beta_l[x]+=exp(logPTrans[y]+logpNoInsgNoIns+beta_l_1[newx]+obs_l_1[newx]);
		}
	}

	// RO -> RO pROgRO=1.0;
	beta_l[ROState]=exp(obs_l_1[ROState]+beta_l_1[ROState]+logpNoInsgNoIns);

	//
	// x^l, i^l=0 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		beta_l[x]+=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns);
	}

	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
	}

	// x^l, i^l=1 = > x^{l+1}=x^l+1, i^{l+1}=0
	beta_l[0+numS]+=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgIns); // cannot go from insertion on to the haplotype
	for (int x=1;x<=hapSize+1;x++ ) {
		int newx=x+1; if (newx>ROState) newx=ROState;
		beta_l[numS+x]+=exp(obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns);
	}


	// convert back to log

	for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}



void ObservationModelFB::passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1)
{
	// P(x^l=x^{l-1}+d|x^{l-1})=(10^{logpSkip*d})
	//for (int x=0;x<2*numS;x++) beta_l[x]=0.0;

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	beta_l[ROState]=( exp(obs_l_1[ROState]+beta_l_1[ROState]+logpLOgLO+logpNoInsgNoIns) )+( exp(obs_l_1[hapSize]+beta_l_1[hapSize]+logpFirstgLO+logpNoInsgNoIns ) );
	for (int x=1;x<=hapSize;x++ ) {
		beta_l[x]=0.0;
		for (int y=1;y<numT;y++) {
			int newx=x-y;
			if (newx<0) newx=0;
			beta_l[x]+=exp(obs_l_1[newx]+logPTrans[y]+beta_l_1[newx]+logpNoInsgNoIns);
		}
	}
	// RO -> RO pROgRO=1.0;
	beta_l[0]=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns);

	// x^l, i^l=0 = > x^{l+1}=x^l-1, i^{l+1}=1
	beta_l[ROState]+=(exp(obs_l_1[numS+ROState]+beta_l_1[numS+ROState]+logpLOgLO+logpInsgNoIns)+exp(obs_l_1[numS+hapSize]+beta_l_1[numS+hapSize]+logpFirstgLO+logpInsgNoIns)); // cannot go from insertion on to the haplotype
	for (int x=0;x<=hapSize;x++) {
		int newx=x-1; if (newx<0) newx=0;
		beta_l[x]+=exp(obs_l_1[numS+newx]+beta_l_1[numS+newx]+logpInsgNoIns);
	}

	/*
	 * from forward OneDec
	alpha_l[numS+ROState]+=( exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpLOgLO+logpInsgNoIns ) );
	alpha_l[numS+hapSize]+=( exp(obs_l_1[ROState]+alpha_l_1[ROState]+logpFirstgLO+logpInsgNoIns ) );
	for (int x=0;x<=hapSize;x++ ) {
		int newx=x-1; if (newx<0) newx=0;
		alpha_l[numS+newx]+=exp(obs_l_1[x]+alpha_l_1[x]+logpInsgNoIns);
	}
	*/
	/*
	for (int x=0;x<=hapSize+1;x++ ) {
		beta_l[x]+=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns);
	}
	*/

	// x^l, i^l=1 = > x^{l+1}, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
	}

	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=0
	for (int x=0;x<=hapSize+1;x++ ) {
		beta_l[numS+x]+=exp(obs_l_1[x]+beta_l_1[x]+logpNoInsgIns);
	}
	// convert back to log
	for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}


void ObservationModelFB::computeForwardMessages()
{
	if (forwardDone) return;
	/*
	for (int b=1;b<readSize;b++) {
			passMessageTwoDec(alpha[b], alpha[b-1], obs[b-1]);
		}
	*/

	for (int b=1;b<=bMid;b++) {
		passMessageTwoDec(alpha[b], alpha[b-1], obs[b-1]);
	}
	for (int b=readSize-1;b>bMid;b--) {
	     passMessageTwoInc(beta[b-1], beta[b], obs[b]);
	}


	forwardDone=true;
}

void ObservationModelFB::computeBackwardMessages()
{
	if (backwardDone) return;
	/*
	for (int b=readSize-1;b>0;b--) {
		passMessageOneDec(beta[b-1], beta[b], obs[b]);
	}
	*/


	for (int b=bMid;b>0;b--) {
		passMessageOneDec(beta[b-1], beta[b], obs[b]);
	}
	for (int b=bMid+1;b<readSize;b++) {
		passMessageOneInc(alpha[b], alpha[b-1], obs[b-1]);
	}


	backwardDone=true;
}

bool ObservationModelFB::_badValue(double v)
{
	if (isnan(v)||isinf(v)) return true; else return false;
}

bool ObservationModelFB::hasErrors()
{
	for (int l=0;l<readSize;l++) {
		for (int x=0;x<hapSize+2;x++) {
			if (_badValue(mar[l][x]) || _badValue(alpha[l][x]) || _badValue(beta[l][x]) ) {
				return true;
			}
		}
	}
	return false;
}


void ObservationModelFB::computeMarginals()
{
	// also normalizes
	if (marginalsComputed) return;
#ifdef DEBUGHMM
	vector<double> logL(readSize);
	cout << "log-likelihoods: ";
#endif
	if (!HMMInitialized) initHMM();
	if (!obsInitialized) setupReadObservationPotentials();
	if (!forwardDone) computeForwardMessages();
	// if (!likelihoodComputed) calcLikelihoodFromLastSlice();
	if (!backwardDone) computeBackwardMessages();

	for (int l=0;l<readSize;l++) {
		double sum=0.0;
		for (int x=0;x<2*numS;x++) {
			mar[l][x]=exp(alpha[l][x]+beta[l][x]+obs[l][x]);
			sum+=mar[l][x];
		}
		for (int x=0;x<2*numS;x++) mar[l][x]/=sum;
		//if (l==0) logLikelihood=log(sum);
#ifdef DEBUGHMM
		logL[l]=log(sum);
		cout << "[" << l << "," << scientific << setprecision(10) << logL[l] << "]";

#endif
	}
#ifdef DEBUGHMM
	cout << "end" << endl;
#endif
	HMMConsistent=true;
	likelihoodComputed=true;
	marginalsComputed=true;
#ifdef DEBUGHMM
	set<int> sl;
	for (int l=0;l<readSize;l++) sl.insert(int(logL[l]*1000.0));
	if (sl.size()!=1) {
		//printMarginals();
		cout << "inconsistent" << endl;
		//throw string("HMM inconsistent!");
	}
#endif
}

void ObservationModelFB::computeXMarginals()
{

	if (!marginalsComputed) computeMarginals();
	for (int b=0;b<readSize;b++) {
		double *_mar=mar[b];
		double *_xmar=xmar[b];
		for (int x=0;x<numS;x++) _xmar[x]=0.0;

		int y=0;
		for (int ins=0;ins<2;ins++) {
			for (int x=0;x<numS;x++) _xmar[x]+=_mar[y++];
		}
	}
	cout.precision(2);

}

void ObservationModelFB::printMarginalsInt( const vector<double*> & pot)
{
	cout.precision(2);

	for (size_t b=0;b<pot.size();b++) {
		cout << "base["  << b << "]: " << endl;
		int y=0;
		for (int ins=0;ins<2;ins++) {
			cout << " ins: " << ins <<  " ";
			for (int x=0;x<numS;x++) {
				stringstream os; os<<fixed <<  (pot[b][y++]);
				string s=os.str();
				cout << string(s,0,5) <<" " ;
			}

			cout << endl;
			}
	}

}

void ObservationModelFB::printMarginals()
{
	cout << "read: " << read.seq <<  " hap: " << hap.seq << endl;
	cout << "logLikelihood: " << logLikelihood << endl;

	/*
	cout << "obs: " << endl;
	printMarginalsInt(obs);
	cout << "alpha: " << endl;
	printMarginalsInt(alpha);
	cout << "beta: " << endl;
	printMarginalsInt(beta);
	*/

	cout << "obs: " << endl;
	printMarginalsInt(obs);
	cout << "alpha: " << endl;
		printMarginalsInt(alpha);
	cout << "mar: " << endl;
	printMarginalsInt(mar);
}


void ObservationModelFB::printStatistics()
{
	cout << "bMid: " << bMid << endl;
	cout << "pTrans: "; for (int x=0;x<numT;x++) cout << exp(logPTrans[x]) << " "; cout << endl;
}

void ObservationModelFB::printAlignment(size_t hapScrPos)
{
	if (!HMMConsistent) computeMarginals();
	computeXMarginals();
	// for every base determine most likely position
	vector<double> maxP(readSize,-HUGE_VAL), entropy(readSize), obsLik(readSize);
	vector<int> maxIdx(readSize, 0);

	int min=ROState+1;
	bool isIncreasing=true;
	for (int b=0;b<readSize;b++) {
		double max=-HUGE_VAL;
		int idx;
		double *m=xmar[b];
		entropy[b]=0.0;
		for (int s=0;s<hapSize+2;s++) {
			entropy[b]+=m[s]*exp(m[s]);
			if (m[s]>max) { max=m[s]; idx=s; };
		}
		maxP[b]=exp(max);
		maxIdx[b]=idx;
		if (b && maxIdx[b]!=0 && maxIdx[b]!=ROState) if (maxIdx[b]-maxIdx[b-1]!=1) isIncreasing=false;
		if (idx<min) min=idx;

		obsLik[b]=obs[b][idx];
		if (b==bMid) obsLik[b]-=prior[idx];

		/*
		if (idx!=0 && idx!=ROState) {
			char hn=hap.seq[idx-1];
			char rn=read.seq[b];

			//if (hn!=rn) cout << "b: " << b << " idx: " << idx << " " << hn << " " << rn << " ol: " << obsLik[b] << endl;
			//if (hn==rn && obsLik[b]<-.5) cout << "b: " << b << " idx: " << idx << " " << hn << " " << rn << " ol: " << obsLik[b] << endl;
		}
		*/
	}

	// number of bases left and right off the haplotype
	//printMarginalsInt(mar);
	size_t nLeft=0, nRight=0;
	for (int b=0;b<readSize;b++) {
		if (maxIdx[b]==0) nLeft++; else if (maxIdx[b]==ROState) nRight++;
	}


	size_t rskip=0;
	if (nLeft>hapScrPos) { rskip=1+nLeft-hapScrPos; nLeft=hapScrPos; };

	size_t offset=hapScrPos-nLeft;

	// print aligned read

	string readString(nLeft+nRight+hapSize+readSize,' ');
	size_t idxL=0, idxR=0;
	for (int b=rskip;b<readSize;b++) {
		char nuc=read.seq[b];
		if (read.qual[b]<params.baseQualThreshold) { nuc=::tolower(nuc); };
			if (maxIdx[b]==0) readString[idxL++]=nuc;
			else if (maxIdx[b]==ROState) readString[nLeft+hapSize+idxR++]=nuc;
			else readString[nLeft+maxIdx[b]-1]=nuc;

	}

	ostringstream os(ostringstream::out);
	os << isIncreasing << " " << logLikelihood << " " << (maxIdx[bMid]!=0 & maxIdx[bMid] !=ROState) << " ";
	size_t s=os.str().size();
	size_t rs,ds;
	if (offset>=s) {
		rs=0;
		ds=offset-s;
	} else {
		rs=s-offset;
		ds=0;
	}


	cout << os.str() << string(ds,' ') << string(readString, rs, readString.size()) << endl;
	/*
	cout << "obsLik:" << string(offset+min-1-7,' '); for (int x=0;x<readSize;x++) cout << int(round(-obsLik[x])); cout << endl;
	cout << "-logq: " << string(offset+min-1-7,' '); for (int x=0;x<readSize;x++) cout << int(round(-log(read.qual[x]))); cout << endl;
	cout << "leq:   " << string(offset+min-1-7,' ');
	for (int b=0;b<readSize;b++) {
		double pr=read.mapQual*read.qual[b];
		double eq=log(.25+.75*pr);
		cout << int(round(-eq));
	}
	cout << endl;
	cout << "luq:   " << string(offset+min-1-7,' ');
	for (int b=0;b<readSize;b++) {
			double pr=read.mapQual*read.qual[b];
			double uq=log(1e-16+1.0-pr)+log(.25);
			cout << int(round(-uq));
	}
	cout << endl;

	for (int b=0;b<readSize;b++) cout << maxIdx[b] << " ";
	cout << " isIncreasing: " << isIncreasing << " hapSize:" << hapSize << " ROState: " << ROState << endl;
	*/
}

ObservationModelFB::~ObservationModelFB()
{
	deleteMemory();
}

ObservationModelFBMax::ObservationModelFBMax(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & _params)
{

	read=r;
	params=_params;
	Init(_hap, hapStart);

}


inline void ObservationModelFBMax::updateMax(double & destValue, int & destIdx, const double newValue, const int  newIdx)
{
	if (newValue>destValue+EPS) {
		destValue=newValue;
		destIdx=newIdx;
	}
	else if (newValue>=destValue && newValue<=destValue+1e-5 && destIdx>newIdx) {
		destValue=newValue;
		destIdx=newIdx;
	}

}


// note that max-product works completely in the log-domain
void ObservationModelFBMax::passMessageTwoInc(double *beta_l, const double *beta_l_1,const double *obs_l_1, int *bt_l)
{

	// x^l, i^l=0 => x^{l+1}, i^{l+1}=0
	//beta_l[0]=( exp(obs_l_1[0]+beta_l_1[0]+logpLOgLO+logpNoInsgNoIns) ) + ( exp(obs_l_1[1]+beta_l_1[1]+logpFirstgLO+logpNoInsgNoIns ) );
	beta_l[0]=-HUGE_VAL;
	updateMax(beta_l[0], bt_l[0], obs_l_1[0]+beta_l_1[0]+logpLOgLO+logpNoInsgNoIns, 0);
	updateMax(beta_l[0], bt_l[0], obs_l_1[1]+beta_l_1[1]+logpFirstgLO+logpNoInsgNoIns, 1);

	for (int x=1;x<=hapSize;x++ ) {
		// double tmp=beta_l_1[x]+logpNoInsgNoIns;
		beta_l[x]=-HUGE_VAL;
		for (int y=1;y<numT;y++) {
			int newx=x+y;
			if (newx>hapSize) newx=ROState;
			//beta_l[x]+=exp(logPTrans[y]+logpNoInsgNoIns+beta_l_1[newx]+obs_l_1[newx]);
			updateMax(beta_l[x], bt_l[x], logPTrans[y]+logpNoInsgNoIns+beta_l_1[newx]+obs_l_1[newx], newx);
		}
	}

	// RO -> RO pROgRO=1.0;
	//beta_l[ROState]=exp(obs_l_1[ROState]+beta_l_1[ROState]+logpNoInsgNoIns);
	beta_l[ROState]=-HUGE_VAL;
	updateMax(beta_l[ROState], bt_l[ROState],obs_l_1[ROState]+beta_l_1[ROState]+logpNoInsgNoIns, ROState);

	//
	// x^l, i^l=0 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		//beta_l[x]+=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns);
		updateMax(beta_l[x], bt_l[x], obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns, numS+x);
	}

	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		//beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
		beta_l[numS+x]=obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns;
		bt_l[numS+x]=numS+x;
	}

	// x^l, i^l=1 = > x^{l+1}=x^l+1, i^{l+1}=0
	//beta_l[0+numS]+=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgIns); // cannot go from insertion on to the haplotype
	updateMax(beta_l[0+numS], bt_l[0+numS], obs_l_1[0]+beta_l_1[0]+logpNoInsgIns, 0);
	for (int x=1;x<=hapSize+1;x++ ) {
		int newx=x+1; if (newx>ROState) newx=ROState;
		//beta_l[numS+x]+=exp(obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns);
		updateMax(beta_l[numS+x], bt_l[numS+x], obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns, newx);
	}


	// convert back to log

	// for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}




/*
void ObservationModelFBMax::passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1, int *bt_l)
{
	// P(x^l=x^{l-1}+d|x^{l-1})=(10^{logpSkip*d})
	//for (int x=0;x<2*numS;x++) beta_l[x]=0.0;

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	//beta_l[ROState]=( exp(obs_l_1[ROState]+beta_l_1[ROState]+logpLOgLO+logpNoInsgNoIns) )+( exp(obs_l_1[hapSize]+beta_l_1[hapSize]+logpFirstgLO+logpNoInsgNoIns ) );

	beta_l[ROState]=-HUGE_VAL;
	updateMax(beta_l[ROState], bt_l[ROState], obs_l_1[ROState]+beta_l_1[ROState]+logpLOgLO+logpNoInsgNoIns, ROState);
	updateMax(beta_l[ROState], bt_l[ROState], obs_l_1[hapSize]+beta_l_1[hapSize]+logpFirstgLO+logpNoInsgNoIns, hapSize);


	for (int x=1;x<=hapSize;x++ ) {
		beta_l[x]=-HUGE_VAL;
		for (int y=1;y<numT;y++) {
			int newx=x-y;
			if (newx<0) newx=0;
			//beta_l[x]+=exp(obs_l_1[newx]+logPTrans[y]+beta_l_1[newx]+logpNoInsgNoIns);
			updateMax(beta_l[x], bt_l[x], obs_l_1[newx]+logPTrans[y]+beta_l_1[newx]+logpNoInsgNoIns, newx);
		}
	}
	// RO -> RO pROgRO=1.0;
	//beta_l[0]=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns);
	beta_l[0]=obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns;
	bt_l[0]=0;

	// x^l, i^l=0 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		//beta_l[x]+=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns);
		updateMax(beta_l[x], bt_l[x], obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns, numS+x);
	}

	// x^l, i^l=1 = > x^{l+1}, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		//beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
		beta_l[numS+x]=obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns;
		bt_l[numS+x]=numS+x;
	}

	// x^l, i^l=1 = > x^{l+1}=x^l+1, i^{l+1}=0
	//beta_l[numS+ROState]+=exp(obs_l_1[ROState]+beta_l_1[ROState]+logpNoInsgIns); // cannot go from insertion on to the haplotype
	updateMax(beta_l[numS+ROState], bt_l[numS+ROState], obs_l_1[ROState]+beta_l_1[ROState]+logpNoInsgIns, ROState);

	for (int x=0;x<=hapSize;x++ ) {
		int newx=x-1; if (newx<0) newx=0;
		//beta_l[numS+x]+=exp(obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns);
		updateMax(beta_l[numS+x], bt_l[numS+x],obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns, newx );
	}
	// convert back to log
	// for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}
*/
void ObservationModelFBMax::passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1, int *bt_l)
{
	// P(x^l=x^{l-1}+d|x^{l-1})=(10^{logpSkip*d})

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	beta_l[ROState]=-HUGE_VAL;
	updateMax(beta_l[ROState], bt_l[ROState], obs_l_1[ROState]+beta_l_1[ROState]+logpLOgLO+logpNoInsgNoIns, ROState);
	updateMax(beta_l[ROState], bt_l[ROState], obs_l_1[hapSize]+beta_l_1[hapSize]+logpFirstgLO+logpNoInsgNoIns, hapSize);

	for (int x=1;x<=hapSize;x++ ) {
		beta_l[x]=-HUGE_VAL;
		for (int y=1;y<numT;y++) { //numT = maxLengthDel
			int newx=x-y;
			if (newx<0) newx=0;
			//beta_l[x]+=exp(obs_l_1[newx]+logPTrans[y]+beta_l_1[newx]+logpNoInsgNoIns);
			updateMax(beta_l[x], bt_l[x], obs_l_1[newx]+logPTrans[y]+beta_l_1[newx]+logpNoInsgNoIns, newx);
		}
	}
	// RO -> RO pROgRO=1.0;
	//beta_l[0]=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns);
	beta_l[0]=obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns;
	bt_l[0]=0;


	// x^l, i^l=0 = > x^{l+1}=x^l-1, i^{l+1}=1
	//beta_l[ROState]+=(exp(obs_l_1[numS+ROState]+beta_l_1[numS+ROState]+logpLOgLO+logpInsgNoIns)+exp(obs_l_1[numS+hapSize]+beta_l_1[numS+hapSize]+logpFirstgLO+logpInsgNoIns)); // cannot go from insertion on to the haplotype
	updateMax(beta_l[ROState], bt_l[ROState],obs_l_1[numS+ROState]+beta_l_1[numS+ROState]+logpLOgLO+logpInsgNoIns, numS+ROState);
	updateMax(beta_l[ROState], bt_l[ROState],obs_l_1[numS+hapSize]+beta_l_1[numS+hapSize]+logpFirstgLO+logpInsgNoIns, numS+hapSize);

	for (int x=0;x<=hapSize;x++) {
		int newx=x-1; if (newx<0) newx=0;
		//beta_l[x]+=exp(obs_l_1[numS+newx]+beta_l_1[numS+newx]+logpInsgNoIns);
		updateMax(beta_l[x], bt_l[x],obs_l_1[numS+newx]+beta_l_1[numS+newx]+logpInsgNoIns, numS+newx);
	}



	// x^l, i^l=1 = > x^{l+1}, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
	//	beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
		beta_l[numS+x]=obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns;
		bt_l[numS+x]=numS+x;
	}

	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=0
	for (int x=0;x<=hapSize+1;x++ ) {
	//	beta_l[numS+x]+=exp(obs_l_1[x]+beta_l_1[x]+logpNoInsgIns);
		updateMax(beta_l[numS+x], bt_l[numS+x],obs_l_1[x]+beta_l_1[x]+logpNoInsgIns, x);
	}
	// convert back to log
	//for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}


void ObservationModelFBMax::runHMM()
{
	if (HMMConsistent) return;
	initHMM();
	setupReadObservationPotentials();
	computeForwardMessages();
	calcLikelihoodFromLastSlice();
	computeMAPState();
	HMMConsistent=true;
}

MLAlignment ObservationModelFBMax::calcLikelihood()
{
	runHMM();
	reportVariants();
  //printMarginals();
	return ml;
}

void ObservationModelFBMax::calcLikelihoodFromLastSlice()
{
	if (likelihoodComputed) return;
	double *alpha_l=alpha[bMid];
	double *beta_l=beta[bMid];
	double *obs_l=obs[bMid];
	logLikelihood=-HUGE_VAL;
	logLikelihoodNoPrior=0.0;
	likOffHap.resize(2);
	likOffHap[0]=-HUGE_VAL;
	likOffHap[1]=-HUGE_VAL;

	int mapStateRMQ=0;

	double llHMQ = -HUGE_VAL;

	vector<double> priorRMQ, priorHMQ;
	computeBMidPrior(priorRMQ, read.mapQual);
	computeBMidPrior(priorHMQ, 1.0-1e-10);
#ifdef NDEBUG
    cout << "bMid, numS = " << bMid << ", " << numS << ", " << endl;
#endif
	int y=0;
	for (int x=0;x<2*numS;x++, y++) {
		double v=alpha_l[y]+obs_l[y]+beta_l[y]+priorRMQ[y];
#ifdef NDEBUG
     cout << "[ " << x << " " << y << "] " << alpha_l[y] << ", " << obs_l[y] << ", " << beta_l[y] << ", " << priorRMQ[y] << ", " << priorHMQ[y];
     cout << ", " << v << endl;
#endif
		if (v>logLikelihood+EPS) {
			logLikelihood=v;
			mapStateRMQ=x;
		}

		if ((x%numS)==0) {
			if (v>likOffHap[0]) likOffHap[0]=v;
		} else if ((x%numS)!=ROState) {
			if (v>likOffHap[1]) likOffHap[1]=v;
		}

		v=alpha_l[y]+obs_l[y]+beta_l[y]+priorHMQ[y];
		if (v>llHMQ+EPS) {
			llHMQ=v;
			mapState[bMid]=x;
		}


	}
#ifdef NDEBUG
	cout << "read: " << bam1_qname(this->read.getBam()) << " read.pos: " << read.pos << " matePos: " << this->read.matePos << " lib: " <<  " prior[" << bMid << ":" << mapState[bMid] << "]: " << priorRMQ[mapStateRMQ] << endl;
	cout << "lib: " << this->read.getLibraryName() << endl;
#endif
	ml.ll=logLikelihood;
	if ((mapState[bMid]%numS)==0 || (mapState[bMid]%numS)==ROState) {
		ml.offHapHMQ=true;
	}else {
		ml.offHapHMQ=false;
	}

	if ((mapStateRMQ%numS)==0 || (mapStateRMQ%numS)==ROState) {
		ml.offHap=true;
	}else {
		ml.offHap=false;
	}
#ifdef NDEBUG
    cout << "mapStateRMQ, numS, mapState[bMid] = " << mapStateRMQ << ", " << numS << ", " << mapState[bMid] << endl;
    cout << "offHapHMQ, offHap = "<< ml.offHapHMQ << "," << ml.offHap << endl;
#endif

	ml.llOff=likOffHap[0];
	ml.llOn=likOffHap[1];
	// now recompute mapState: we want to only show alignments to the haplotype

#ifdef DEBUGHMM
	//cout << "calcLikelihoodFromLastSlice(): " << logLikelihood << endl;
	//cout << "here: " << scientific << setprecision(10) <<  log(likOffHap[0]*exp(prior[0])+likOffHap[1]*exp(prior[1])) << " " << logLikelihood << endl;
#endif
	likelihoodComputed=true;
}



void ObservationModelFBMax::computeMAPState()
{
    //cout << "computeMAPState" << endl;
	// now backtrack
	for (int b=bMid; b>0;b--) {
		//cout << "mapState[" << b << "]: " << mapState[b] << " btf[]: " << btf[b][mapState[b]] << endl;
		mapState[b-1]=btf[b][mapState[b]];
	}

	for (int b=bMid;b<readSize-1;b++) {
		mapState[b+1]=btb[b][mapState[b]];
		//cout << "mapState[" << b << "]: " << mapState[b] << " btf[]: " << btb[b][mapState[b]] << endl;
	}

	//cout << "mapState: "; for (int b=0;b<readSize;b++) cout << " " << mapState[b]; cout << endl;

}

/*
void ObservationModelFBMax::reportVariants(map<int, ReportVariant > & indels, map<int, ReportVariant > & snps, string & align)
{
	runHMM();

	align=string(hapSize, 'R');
	indels.clear();
	snps.clear();

	ml.firstBase=-1;
	ml.lastBase=-1;

	int b=0;
	while (b<readSize) {
		// only report variants for bases that are on the haplotype
		int s=mapState[b];
		if ( (s%numS)>0 && (s%numS)<=hapSize ) {
			if (s>=numS) { // insertion
				int pos=(s%numS)-1; // position of insertion wrt haplotype
				int len=0; // length of insertion
				int rpos=b; // start base of insertion in read
				while (b<readSize && mapState[b]>=numS) {
					b++;
					len++;
				}
				indels[pos]=ReportVariant(len, read.seq.seq.substr(rpos, len), b-len);
				indels[pos]=ReportVariant()
			} else {
				// update firstBase and lastBase
				if (ml.firstBase==-1) ml.firstBase=s-1; else if (s-1<ml.firstBase) ml.firstBase=s-1;
				if (ml.lastBase==-1) ml.lastBase=s-1; else if (s-1>ml.lastBase) ml.lastBase=s-1;


				// check for SNP
				if (read.seq[b]!=hap.seq[s-1]) {
					string snp;
					snp+=hap.seq[s-1];
					snp.append("=>");
					snp+=read.seq[b];
					snps[s-1]=(ReportVariant(0,snp, b));
					align[s-1]=read.seq[b];
				}
				// check for deletion
				if (b<readSize-1) {
					int ns=mapState[b+1];
					if (ns<numS && ns-s>1) { // make sure next state is not an insertion..
						int pos=s+1-1;
						int len=-(ns-s-1);
						indels[pos]=ReportVariant(len, hap.seq.substr(pos, -len), b);
						for (int y=pos;y<-len+pos;y++) align[y]='D';
					}
				}

			}

		}

		b++;
	}


}
*/
/*
void getUniqueCoordinates(const Haplotype & hap, const Read & read, AlignedVariant & av)
{
	int rightFlankHap, leftFlankHap,rightFlankRead, leftFlankRead;

	if (av.getType()==Variant::INS || av.getType() == Variant::DEL) {
		const string & seq = av.getSeq();
		int l = seq.size();

		if (1) {


			int p = av.getStartHap();
			cout << "startHap: " << p << endl;
			cout << "startRead: " << av.getStartRead() << endl;

			while (p+l<=hap.seq.size()) {
				string ss = hap.seq.substr(p,l);
				if (ss!=seq) {

					break;
				}
				p+=l;
			}

			rightFlankHap = p;

			p = av.getStartHap()-l;

			while (p>=0) {
				string ss = hap.seq.substr(p,l);
				if (ss!=seq) {
					p+=l;
					break;
				}
				p-=l;
			}

			leftFlankHap = p-1;
			if (leftFlankHap<0) leftFlankHap = 0;

			cout << "leftFlankHap: " << leftFlankHap << " rightFlankHap: " << rightFlankHap << endl;
		}

		// do read flanks

		if (av.getType() == Variant::INS) {
			int p = av.getStartRead();
			cout << "startRead: " << p << endl;
			while (p+l<=read.seq.size()) {
				string ss = read.seq.seq.substr(p,l);
				if (ss!=seq) {

					break;
				}
				p+=l;
			}

			rightFlankRead = p;

			p = av.getStartRead()-l;

			while (p>=0) {
				string ss = read.seq.seq.substr(p,l);
				if (ss!=seq) {
					p+=l;
					break;
				}
				p-=l;
			}

			int leftFlankRead = p-1;
			if (leftFlankRead<0) leftFlankRead = 0;

			cout << "leftFlankRead: " << leftFlankRead << " rightFlankRead: " << rightFlankRead << endl;



		} else if (av.getType() == Variant::DEL) {
					int p = av.getStartRead()+1;
					cout << "startRead: " << p << endl;
					while (p+l<=read.seq.size()) {
						string ss = read.seq.seq.substr(p,l);
						if (ss!=seq) {

							break;
						}
						p+=l;
					}

					rightFlankRead = p;

					p = av.getStartRead()+1-l;

					while (p>=0) {
						string ss = read.seq.seq.substr(p,l);
						if (ss!=seq) {
							p+=l;
							break;
						}
						p-=l;
					}

					int leftFlankRead = p-1;
					if (leftFlankRead<0) leftFlankRead = 0;

					cout << "leftFlankRead: " << leftFlankRead << " rightFlankRead: " << rightFlankRead << endl;



				}


	}



}
*/


void ObservationModelFBMax::reportVariants()
{
	runHMM();

	ml.align=string(hapSize, 'R');
	ml.indels.clear();
	ml.snps.clear();

	ml.firstBase=-1;
	ml.lastBase=-1;
	ml.hapIndelCovered.clear();
	ml.hapSNPCovered.clear();
	ml.hpos.clear();
	ml.hpos.resize(readSize);

	ml.nBQT=0;
	ml.nmmBQT=0;
	ml.mLogBQ=0.0;
	ml.nMMRight=0;
	ml.nMMLeft=0;
	ml.numIndels = 0;
	ml.numMismatch = 0;

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
				ml.numIndels++;
				b--;
				//getFlankingCoordinatesBetter(this->hap, this->read, ml.indels[pos]);

			} else {
				ml.hpos[b]=s-1;
				// update firstBase and lastBase
				if (ml.firstBase==-1) ml.firstBase=s-1; else if (s-1<ml.firstBase) ml.firstBase=s-1;
				if (ml.lastBase==-1) ml.lastBase=s-1; else if (s-1>ml.lastBase) ml.lastBase=s-1;

				if (read.qual[b]>params.checkBaseQualThreshold){
					ml.nBQT++;
					ml.mLogBQ+=log10(1.0-read.qual[b]);
				}

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

					if (read.qual[b]>params.checkBaseQualThreshold) {
						ml.nmmBQT++;

					}

					if (b<6) ml.nMMLeft++;
					if (b>readSize-6) ml.nMMRight++;


					if (read.qual[b]>0.95) ml.numMismatch++;



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
						//getFlankingCoordinatesBetter(this->hap, this->read, ml.indels[pos]);
						ml.numIndels++;
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

void ObservationModelFBMax::printAlignment(size_t hapScrPos)
{
	// count how many bases in the read are left of the haplotype
	calcLikelihood();




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
	os << readSize << " " << ml.offHap << " " << ml.indels.size() << " " << ml.firstBase << " " << ml.lastBase << " " << logLikelihood << " ";
	for (map<int,AlignedVariant>::const_iterator it=hap.indels.begin();it!=hap.indels.end();it++) {
        os << it->second.getString();
		if  (ml.hapIndelCovered[it->first]) os << "1 "; else os << "0 ";
	}
	string prefix=os.str();

	int leftHapSpace=int(hapScrPos);
	if (leftHapSpace<0) leftHapSpace=0;

	string prLeftHap=string(leftHapSpace,' ');

	if (int(leftHap.size())>leftHapSpace) {
		prLeftHap=leftHap.substr(leftHap.size()-leftHapSpace, leftHapSpace);
	} else if (leftHap.size()>0) {
		prLeftHap.replace(leftHapSpace-leftHap.size(), leftHap.size(), leftHap);
	}

	cout <<prLeftHap<<rhap<<rightHap << " " << endl << ins << " read: " << read.seq.seq << endl;
    cout << prefix << endl;
    /*
    for(int i=0;i < read.qual.size();i++) {
        cout << read.qual[i] << " ";
    }
    cout << endl;

	for (map<int,AlignedVariant>::const_iterator it=hap.indels.begin();it!=hap.indels.end();it++) {
		cout << " " << it->first;
	}
	cout << endl;
    for (int b=0;b<readSize;b++) {
		int s=mapState[b];
		cout << "[" << b << " " << read.seq.seq[b] << " " << hap.seq[(s%numS)-1] << " " << s << "]";
	}
     */
}



void ObservationModelFBMax::computeForwardMessages()
{
	if (forwardDone) return;

	for (int b=1;b<=bMid;b++) {
		passMessageTwoDec(alpha[b], alpha[b-1], obs[b-1], btf[b]);
	}
	for (int b=readSize-1;b>bMid;b--) {
	     passMessageTwoInc(beta[b-1], beta[b], obs[b], btb[b-1]);
	}

	forwardDone=true;
}

void ObservationModelFBMax::computeBackwardMessages()
{
	// no backward messages for this model
	backwardDone=true;
}

void ObservationModelFBMax::allocateMemory()
{
	if (memAllocated) return; //throw string("Memory already allocated.");
	mapState.resize(readSize, 0);


	obs.reserve(readSize);
	alpha.reserve(readSize);
	beta.reserve(readSize);
	btf.reserve(readSize);
	btb.reserve(readSize);

	for (int b=0;b<readSize;b++) {
		obs.push_back(new double[numS*2]);
		alpha.push_back(new double[numS*2]);
		beta.push_back(new double[numS*2]);
		if (b<=bMid) {
			btf.push_back(new int[numS*2]);
		} else btf.push_back(NULL);
		if (b>=bMid) {
			btb.push_back(new int[numS*2]);
		} else btb.push_back(NULL);
	}
	if (makeObsVector) { obsVector=new double[4*readSize]; };
	memAllocated=true;
}

void ObservationModelFBMax::deleteMemory()
{
	if (memAllocated) {
		for (int b=0;b<readSize;b++) {
			delete[] obs[b];
			delete[] alpha[b];
			delete[] beta[b];
			if (btf[b]!=NULL) delete[] btf[b];
			if (btb[b]!=NULL) delete[] btb[b];
		}
		if (makeObsVector) delete[] obsVector;
		memAllocated=false;
	}
}

ObservationModelFBMaxErr::ObservationModelFBMaxErr(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & _params)
{

	read=r;
	params=_params;
	Init(_hap, hapStart);
}



void ObservationModelFBMaxErr::setupTransitionProbs()
{
	logpLOgLO=log(1.0-params.pFirstgLO);
	logpFirstgLO=log(params.pFirstgLO);

	numT=params.maxLengthDel+2;
	logPTrans.resize(numT);
	// maxT is the transition which corresponds to a normal-operation base extension
	logPTrans[1]=log(1.0-params.pError);
	double norm=0.0;
	for (int x=1;x<numT;x++) if (x!=1) {
		double p=-fabs(1.0-double(x));
		logPTrans[x]=p;
		norm+=exp(p);
	}
	norm=log(norm/params.pError);
	for (int x=1;x<numT;x++) if (x!=1) logPTrans[x]-=norm;

	// check norm
	norm=0.0;
	for (int x=1;x<numT;x++) norm+=exp(logPTrans[x]);
	assert(fabs(norm-1.0)<1e-15);

	logpInsgIns=-.5 ;
	logpNoInsgIns=log(1.0-exp(logpInsgIns));
	logpInsgNoIns=log(params.pError);
	logpNoInsgNoIns=log(1-params.pError);
	/*
	cout << "logpInsgIns: " << logpInsgIns << endl;
	cout << "logpNoInsgIns: " << logpNoInsgIns << endl;
	cout << "logpInsgNoIns: " << logpInsgNoIns << endl;
	cout << "logpNoInsgNoIns: " << logpNoInsgNoIns << endl;
	*/

	// determine base-specific error probabilities
	ReadIndelErrorModel riem;

	logProbError = vector<double>(hapSize+2,log(1e-5));
	logProbNoError = vector<double>(hapSize+2,log(1-1e-5));


	int len=1;
	double perr=riem.getViterbiHPError(1);
	logProbError[1]=log(perr);
	logProbNoError[1]=log(1.0-perr);

	// NOTE X = ( LO, 0, 1, 2, 3, .. )
	for (int b=1;b<hapSize;b++) {
		if (hap.seq[b]==hap.seq[b-1]) {
			len++;
		} else {
			perr=riem.getViterbiHPError(len);
//			cout << "len: " << len << " perr: " << perr << endl;
			logProbError[b]=log(perr);
			logProbNoError[b]=log(1.0-perr);
			len=1;
		}
		//cout << "hap[" << b << "]: " << len << " " << logProbError[b+1]  << endl;
	}
	perr=riem.getViterbiHPError(len);
	// cout << "len: " << len << " perr: " << perr << endl;
	logProbError[hapSize-1]=log(perr);
	logProbNoError[hapSize-1]=log(1.0-perr);


	/*
	cout << "logProbError: " << endl;
	for (int x=0;x<=hapSize+1;x++) {
		cout << "x: " << x << " " << ((x>0&&x<=hapSize)?hap.seq[x-1]:'N') << " " << logProbError[x]  << " " << logProbNoError[x] <<endl;
	}
	*/

}

void ObservationModelFBMaxErr::passMessageTwoInc(double *beta_l, const double *beta_l_1,const double *obs_l_1, int *bt_l)
{           //                                       b-1                           b

	// x^l, i^l=0 => x^{l+1}, i^{l+1}=0
	//beta_l[0]=( exp(obs_l_1[0]+beta_l_1[0]+logpLOgLO+logpNoInsgNoIns) ) + ( exp(obs_l_1[1]+beta_l_1[1]+logpFirstgLO+logpNoInsgNoIns ) );
	beta_l[0]=-HUGE_VAL;
	updateMax(beta_l[0], bt_l[0], obs_l_1[0]+beta_l_1[0]+logpLOgLO+logpNoInsgNoIns, 0);
	updateMax(beta_l[0], bt_l[0], obs_l_1[1]+beta_l_1[1]+logpFirstgLO+logpNoInsgNoIns, 1);

	for (int x=1;x<=hapSize;x++ ) {
		// double tmp=beta_l_1[x]+logpNoInsgNoIns;
		beta_l[x]=-HUGE_VAL;
		for (int y=1;y<numT;y++) {
			int newx=x+y;
			if (newx>hapSize) newx=ROState;
			double lpn=logProbNoError[newx];
			double lpt=logProbError[newx];
			double lp=(y==1)?lpn:(lpt+double(y-1)*logpInsgIns);

			//beta_l[x]+=exp(logPTrans[y]+logpNoInsgNoIns+beta_l_1[newx]+obs_l_1[newx]);
			updateMax(beta_l[x], bt_l[x], lp+lpn+beta_l_1[newx]+obs_l_1[newx], newx);
		}
	}

	// RO -> RO pROgRO=1.0;
	//beta_l[ROState]=exp(obs_l_1[ROState]+beta_l_1[ROState]+logpNoInsgNoIns);
	beta_l[ROState]=-HUGE_VAL;
	updateMax(beta_l[ROState], bt_l[ROState],obs_l_1[ROState]+beta_l_1[ROState]+logProbNoError[ROState], ROState);

	//
	// x^l, i^l=0 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize;x++ ) {
		//beta_l[x]+=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgNoIns);
		updateMax(beta_l[x], bt_l[x], obs_l_1[numS+x]+beta_l_1[numS+x]+logProbError[x+1], numS+x);
	}
	int x=hapSize+1; updateMax(beta_l[x], bt_l[x], obs_l_1[numS+x]+beta_l_1[numS+x], numS+x);


	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
		//beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
		beta_l[numS+x]=obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns;
		bt_l[numS+x]=numS+x;
	}

	// x^l, i^l=1 = > x^{l+1}=x^l+1, i^{l+1}=0
	//beta_l[0+numS]+=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgIns); // cannot go from insertion on to the haplotype
	updateMax(beta_l[0+numS], bt_l[0+numS], obs_l_1[0]+beta_l_1[0]+logpNoInsgIns, 0);
	for (int x=1;x<=hapSize+1;x++ ) {
		int newx=x+1; if (newx>ROState) newx=ROState;
		//beta_l[numS+x]+=exp(obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns);
		updateMax(beta_l[numS+x], bt_l[numS+x], obs_l_1[newx]+beta_l_1[newx]+logpNoInsgIns, newx);
	}


	// convert back to log

	// for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}

void ObservationModelFBMaxErr::passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1, int *bt_l)
{ // 														b                 b-1
	// P(x^l=x^{l-1}+d|x^{l-1})=(10^{logpSkip*d})

	// x^l, i^l=0 = > x^{l+1}, i^{l+1}=0
	beta_l[ROState]=-HUGE_VAL;
	updateMax(beta_l[ROState], bt_l[ROState], obs_l_1[ROState]+beta_l_1[ROState]+logpLOgLO+logpNoInsgNoIns, ROState);
	updateMax(beta_l[ROState], bt_l[ROState], obs_l_1[hapSize]+beta_l_1[hapSize]+logpFirstgLO+logpNoInsgNoIns, hapSize);

	for (int x=1;x<=hapSize;x++ ) {
		beta_l[x]=-HUGE_VAL;
		double lpt = logProbError[x];
		double lpn = logProbNoError[x];
		for (int y=1;y<numT;y++) {
			int newx=x-y;
			if (newx<0) newx=0;
			double lp=(y==1)?lpn:(lpt+double(y-1)*logpInsgIns);
			//beta_l[x]+=exp(obs_l_1[newx]+logPTrans[y]+beta_l_1[newx]+logpNoInsgNoIns);
			updateMax(beta_l[x], bt_l[x], obs_l_1[newx]+lp+beta_l_1[newx]+lpn, newx);
		}
	}
	// RO -> RO pROgRO=1.0;
	//beta_l[0]=exp(obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns);
	beta_l[0]=obs_l_1[0]+beta_l_1[0]+logpNoInsgNoIns;
	bt_l[0]=0;


	// x^l, i^l=0 = > x^{l+1}=x^l-1, i^{l+1}=1
	//beta_l[ROState]+=(exp(obs_l_1[numS+ROState]+beta_l_1[numS+ROState]+logpLOgLO+logpInsgNoIns)+exp(obs_l_1[numS+hapSize]+beta_l_1[numS+hapSize]+logpFirstgLO+logpInsgNoIns)); // cannot go from insertion on to the haplotype
	updateMax(beta_l[ROState], bt_l[ROState],obs_l_1[numS+ROState]+beta_l_1[numS+ROState]+logpLOgLO+logProbError[ROState], numS+ROState);
	updateMax(beta_l[ROState], bt_l[ROState],obs_l_1[numS+hapSize]+beta_l_1[numS+hapSize]+logpFirstgLO+logProbError[hapSize], numS+hapSize);

	for (int x=1;x<=hapSize;x++) {
		int newx=x-1; if (newx<0) newx=0;
		//beta_l[x]+=exp(obs_l_1[numS+newx]+beta_l_1[numS+newx]+logpInsgNoIns);
		updateMax(beta_l[x], bt_l[x],obs_l_1[numS+newx]+beta_l_1[numS+newx]+logProbError[x], numS+newx);
	}



	// x^l, i^l=1 = > x^{l+1}, i^{l+1}=1
	for (int x=0;x<=hapSize+1;x++ ) {
	//	beta_l[numS+x]=exp(obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns);
		beta_l[numS+x]=obs_l_1[numS+x]+beta_l_1[numS+x]+logpInsgIns;
		bt_l[numS+x]=numS+x;
	}

	// x^l, i^l=1 = > x^{l+1}=x^l, i^{l+1}=0
	for (int x=1;x<=hapSize+1;x++ ) {
	//	beta_l[numS+x]+=exp(obs_l_1[x]+beta_l_1[x]+logpNoInsgIns);
		updateMax(beta_l[numS+x], bt_l[numS+x],obs_l_1[x]+beta_l_1[x]+logpNoInsgIns, x);
	}
	// convert back to log
	//for (int x=0;x<2*numS;x++) beta_l[x]=log(beta_l[x]);
}
