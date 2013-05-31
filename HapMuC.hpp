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
#ifndef HAPMUC_HPP_
#define HAPMUC_HPP_
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <ext/hash_map>

#include "MyBam.hpp"
#include "faidx.h"
#include "Haplotype.hpp"
#include "ObservationModel.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
//#include "Fast.hpp"
#include "MLAlignment.hpp"
#include "Read.hpp"
#include "StringHash.hpp"

#include "OutputData.hpp"
#include "Library.hpp"
#include "VariantFile.hpp"

const int SHIFTSTRAND = 1000000; // used to keep track of forward and reverse matches in ::filterHaplotypes

using namespace std;
using namespace boost;
using __gnu_cxx::hash;
typedef struct
{
	double pOff, pOn;
} HapReadLik;

struct lower_bound_t {	
	double lower_bound;
	double ln_p_x_given_z;
    double ln_p_z_given_pi;
    double ln_p_pi;
    double ln_q_z;
    double ln_q_pi;
    double ln_prior;
    vector<int> compatible;
public:
	lower_bound_t(){			
        lower_bound = 0.0;
        ln_p_x_given_z = 0.0;
        ln_p_z_given_pi = 0.0;
        ln_p_pi = 0.0;
        ln_q_z = 0.0;
        ln_q_pi = 0.0;
	}
    
    void printCOut() const {
        cout << "lb=" << lower_bound << ", ln_p_x_given_z = " << ln_p_x_given_z 
        << ", ln_p_z_given_pi = " << ln_p_z_given_pi << ", ln_p_pi = " << ln_p_pi << ", ln_q_z = " << ln_q_z
        << ", ln_q_pi = " << ln_q_pi << "ln_prior" << ln_prior << endl;
    }
};

#include <iostream>
#include <string>
#include <map>
using namespace std;

#include <sstream>
#include <fstream>
#include <utility>
#include <string.h>

template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}


class VariantCoverage {
public:
	VariantCoverage()
	{
		nf=0;
		nr=0;
	}
	VariantCoverage(int _nf, int _nr) 
	{
		nf = _nf;
		nr = _nr;
	}
	int nf, nr; // forward and reverse
};
class HapPairLik {
public:
    double ll;
    int h1, h2;
    int numIndFirst;
    int numIndSecond;
    int numFirst, numSecond; // number of reads mapped to first and second
    int numOffBoth; // number of reads that do not map to either haplotype
    double numOffBothError;
    map<int, VariantCoverage> hapIndelCoverage1, hapSNPCoverage1,hapIndelCoverage2, hapSNPCoverage2; // indels and snps in the _haplotype_ covered by the read
    
    operator double() const { return ll;};
};

class HapEstResult {
public:
    HapEstResult();
    VariantInfo info;
    HapEstResult(const AlignedVariant & _av, int _pos, double _prob, double _freq, int _nrf, int _nrr, VariantInfo vi, double _N1, double _N2, double _N3, double _N4) {
        av=_av;
        pos=_pos;
        prob=_prob;
        freq=_freq;
        nrf=_nrf;
        nrr=_nrr;
        info = vi;
        N1=_N1;
        N2=_N2;
        N3=_N3;
        N4=_N4;
    };
    AlignedVariant av;
    int pos;
    double prob;
    double freq;
    int nrf; // number of reads on reverse strand
    int nrr; // number of reads on forward strand
    double N1, N2, N3, N4;
    void printCOut() const {
        cout << "[" << av.getStartHap() << ", " << av.getString() << "] " 
        << pos << " " << prob << " " << freq << " " << nrf << " " << nrr << endl;
    }
};

class Parameters {
public:
    Parameters(const string & _tid, string _fileName, const string & modelType) : obsParams(modelType)
    {
        tid=_tid;
        fileName=_fileName;
        setDefaultValues();
    }
    void setDefaultValues()
    {
        bayesa0=0.001;
        width=30;
        maxHap=100;
        skipMaxHap=1000;
        maxReads=500;
        mapQualThreshold=0.995;
        glfNumHap=5;
        inferenceMethod="empirical";
        minReadOverlap=5;
        minCount=2;
        maxReadLength=40;
        numOutputTopHap=5;
        checkAllCIGARs=1;
        bayesType="all";
        
        fastWidth=4;
        analyzeLowFreq=false;
        analyzeLowFreqDiffThreshold=1.0;
        showHapDist=true;
        showHapAlignments=false;
        showCandHap=false;
        showReads=false;
        fastWidthOverlap=4;
        noIndelWindow=-1;
        mapUnmappedReads=false;
        priorIndel=1.0/10000;
        priorSNP=1.0/1000.0;
        filterReadAux=string("");
        quiet=true;
        computeML=false;
        computeMAP=false;
        doDiploid=false;
        slower=true;
        estimateHapFreqs=false;
        printCallsOnly=true;
        outputPooledLikelihoods=false;
        filterHaplotypes = false;
        vector<double> tmpa(3,1.0), tmpb(2,1.0), tmpc(2,1.0);
        a0 = tmpa;b0 = tmpb;c0 = tmpc;
        outputGLF=true;
        outputRealignedBAM=false;
        processRealignedBAM="no";
        changeINStoN = false;
        
        
        EMtol=1e-4;
    }
    
    OutputData makeOutputData(ostream & out)
    {
        OutputData oData(out);
        oData("msg")("index");
        oData("analysis_type");
        oData("tid")("lpos")("rpos")("center_position")("realigned_position");
        oData("ref_all")("num_reads")("num_hqreads");
        oData("post_prob_variant")("est_freq")("was_candidate_in_window");
        
        oData("num_mapped_to_first")("num_mapped_to_second");
        oData("num_off_hap")("loglik_hap_pair")("loglik_next_hap_pair");
        oData("first_var_cover_forward")("first_var_cover_reverse")("second_var_cover_forward")("second_var_cover_reverse");
        oData("first_called_all")("second_called_all")("loglik_called_genotype")("loglik_ref_ref")("alt_genotypes");
        return oData;
    }
    
    OutputData makeGLFOutputData(ostream & out)
    {
        OutputData oData(out);
        oData("msg")("index");
        oData("analysis_type");
        oData("tid")("lpos")("rpos")("center_position")("realigned_position")("was_candidate_in_window");
        oData("ref_all")("nref_all")("num_reads");
        oData("post_prob_variant")("qual")("est_freq")("logZ")("hapfreqs");
        
        oData("indidx")("msq")("numOffAll")("num_indel")("num_cover_forward")("num_cover_reverse")("num_unmapped_realigned");
        oData("var_coverage_forward")("var_coverage_reverse");
        oData("nBQT")("nmmBQT")("mLogBQ")("nMMLeft")("nMMRight");
        oData("glf");
        return oData;
    }
    
    OutputData makeMutationData(ostream & out) {
        OutputData oData(out);
        oData("chr")("start")("end");
        oData("ref")("obs");
        oData("ref_count_tumor")("obs_count_tumor")("ref_count_normal")("obs_count_normal");
        oData("missrate_tumor")("strandrate_tumor")("missrate_normal")("strandrate_normal");
        oData("TN1")("TN2")("TN3")("TN4");
        oData("NN1")("NN2")("NN3")("NN4");
        oData("fisher");
        oData("hap2_bf")("bf2");
        return oData;
    }
    
    OutputData makeGLFv2OutputData(ostream & out)
    {
        OutputData oData(out);
        oData("msg")("index");
        oData("analysis_type");
        oData("tid")("candidate_position")("realigned_position");
        oData("ref_all")("nref_all")("num_reads");
        oData("post_prob_variant")("est_freq");
        
        oData("indidx")("msq")("num_cover_forward")("num_cover_reverse");
        oData("glf");
        return oData;
    }   
    
    
    void print()
    {
        cout << "HapMuC parameters: " << endl;
        cout << "\ttid: " << tid << " width: " << width << " maxHap: " << maxHap << " maxReads: " << maxReads << " skipMaxHap: " << skipMaxHap << endl;
        cout << "\toutputFilename: " << fileName << endl;
        cout << "\tmapQualThreshold: " << mapQualThreshold << endl;
        //cout << "\tscaleError: " << scaleErr << endl;
        cout << "\tinferenceMethod: " << inferenceMethod << endl;
        //cout << "\tglfNumHap: " << glfNumHap << endl;
        cout << "\tanalyzeLowFreq: " << analyzeLowFreq << endl;
        cout << "\tanalyzeLowFreqDiffThreshold: " << analyzeLowFreqDiffThreshold << endl;
        cout << "\tshowHapDist: " << showHapDist << endl;
        cout << "\tminReadOverlap: " << minReadOverlap << endl;
        cout << "\tmaxReadLength: " << maxReadLength << endl;
        //cout << "\tminCount: " << minCount << endl;
        cout << "\tmaxHapReadProd: " << maxHapReadProd << endl;
        //cout << "\tfastWidth: " << fastWidth << endl;
        //cout << "\tfastWidthOverlap: " << fastWidthOverlap << endl;
        cout << "\tshowCandHap: " << showCandHap << endl;
        cout << "\tshowReads: " << showReads << endl;
        cout << "\tfilterHaplotypes: " << filterHaplotypes << endl;
        cout << "\tnoIndelWindow: " << noIndelWindow << endl;
        cout << "\tmapUnmappedReads: " << mapUnmappedReads << endl;
        
        cout << "\tnumOutputTopHap: " << numOutputTopHap << endl;
        
        cout << "\tcheckAllCIGARs: " << checkAllCIGARs << endl;
        cout << "\tchangeINStoN: " << changeINStoN << endl;
        
        cout << endl;
        cout << "\tquiet: " << quiet << endl;
        cout << "\tprintCallsOnly: " << printCallsOnly << endl;
        cout << "\tfaster: " << !slower << endl;
        cout << "\tdoDiploid: " << doDiploid << endl;
        cout << "\tdoEM: " << estimateHapFreqs << endl;
        
        cout << "\toutputPooledLikelihoods: " << outputPooledLikelihoods << endl;
        cout << "\toutputRealignedBAM: " << outputRealignedBAM << endl;
        cout << "\tprocessRealignedBAM: " << processRealignedBAM << endl;
        cout << "\tshowHapAlignments: " << showHapAlignments << endl;
        
        cout << "\tEM tol: " << EMtol << endl;
        cout << "\tbayesEM a0: " << bayesa0 << endl;
        cout << "\tbayesType: " << bayesType << endl;
        
        
        cout << "\tpriorIndel: " << priorIndel << endl;
        cout << "\tpriorSNP: " << priorSNP << endl;
        
        //cout << "\tmeanInsert: " << meanInsert << endl;
        //cout << "\tstdInsert: " << stdInsert << endl;
        
        cout << "\tfilterReadAux: " << filterReadAux << endl;
        
        cout << "Observation model parameters: " << endl;
        obsParams.print();
    }
    int noIndelWindow, numOutputTopHap, checkAllCIGARs, minReadOverlap, maxHapReadProd;
    uint32_t width, maxHap, maxReads, skipMaxHap, glfNumHap,  maxReadLength, minCount, fastWidth, fastWidthOverlap;
    double checkBaseQualThreshold;
    double mapQualThreshold, scaleErr, priorIndel, priorSNP, EMtol, bayesa0;
    string fileName, inferenceMethod, refFileName, tid, filterReadAux, bayesType, processRealignedBAM;
    bool analyzeLowFreq, showHapDist, showCandHap, showReads, showHapAlignments, alignAgainstReference, mapUnmappedReads, quiet, estimateHapFreqs, doDiploid, computeML, computeMAP, slower,printCallsOnly, outputPooledLikelihoods, filterHaplotypes;
    bool outputRealignedBAM, outputGLF, varFileIsOneBased, changeINStoN;
    double analyzeLowFreqDiffThreshold;
    double meanInsert, stdInsert;
    vector<double> a0, b0, c0;
    ObservationModelParameters obsParams;
};

class HapMuC
{
public:
	string getRefSeq(uint32_t lpos, uint32_t rpos);
    string get_var_symbol(string ref, string obs);
    AlignedCandidates getCandidateVariants(string line, vector<AlignedVariant> &close_somatic, vector<AlignedVariant> &close_germline);
    void mutationCall(const string & variantsFileName); 
	void addLibrary ( const string & name, const Library & lib) { libraries[name.c_str()]=lib; }
	void addLibrary ( const string & fileName) { libraries.addFromFile(fileName); }

protected:
	void outputHapsAndFreqs(ostream *output, const string & prefix, const vector<Haplotype> & haps, const vector<double> & freqs, uint32_t leftPos);
	void outputParams(int iter, int leftPos, const vector<int> & compatible, const vector<double> & resp, const vector<double> & ak, const vector<double> & ln_p_x_given_h, const vector<double> & lpi) ;
	void debug(const pair<Haplotype, Haplotype> & hp, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos);
	void debug(const pair<Haplotype, Haplotype> & hp1, const pair<Haplotype, Haplotype> & hp2, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos);
	void filterHaplotypes(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks,  vector<int> & filtered,  map<pair<int, AlignedVariant>, VariantCoverage> & varCoverage, bool doFilter);
	vector<MyBam *> myBams;
	vector<string> myBamsFileNames;
    void getReadsFromBams(vector<MyBam *> & Bams, uint32_t leftPos, uint32_t rightPos, vector<Read> & reads, uint32_t & oldLeftPos, uint32_t  & oldRightFetchReadPos, vector<Read *> & readBuffer, const bool reset);
    
    vector<MyBam *> tumorBams;
    vector<string> tumorBamsFileNames;
    vector<MyBam *> normalBams;
    vector<string> normalBamsFileNames;
    
    vector<AlignedVariant> parse_close_vars(string s);
    
	LibraryCollection libraries;

	class CIGAR : public vector<pair<int,int> >
	{
	public:
		typedef pair<int,int> CIGOp;
		int refPos;
	};
	CIGAR getCIGAR(const Haplotype & hap, const Read & read, const MLAlignment & ml, int refSeqStart);
	void writeRealignedBAMFile(const string & fileName, const vector<CIGAR> & cigars, const vector<Read> & reads, const vector<int> & onHap, const bam_header_t *bh);
	void writeUnalignedBAMFile(const string & fileName, const vector<Read> & reads, const vector<int> & onHap, const bam_header_t *bh);
	class InDel {
	public:
		InDel()
		{
			count[0]=0;
			count[1]=0;
		}
		typedef enum { In, Del} Type;
		Type type;
		size_t count[2];
	};

public:

	Parameters params;

	HapMuC(const string & bfName, const Parameters & _params, int multipleFiles);
    HapMuC(const string & normalBF, const string & tumorBF, const Parameters & _params);
	~HapMuC();


	map<uint32_t, InDel> indels;
	class ScanStats
	{
		public:
		ScanStats()
		{
			numUnmappedMate=0;
		}
		int numUnmappedMate;
	};
	ScanStats scanStats;
protected:
	faidx_t *fai;

};



class FFData
{
public:
	uint32_t start, end;
	HapMuC *det;
	map<string, int> unmappedMate;
	map<int, int > insHisto, delHisto;
};

#endif /*HapMuC_HPP_*/