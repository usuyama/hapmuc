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
#ifndef OBSERVATIONMODELFB_HPP_
#define OBSERVATIONMODELFB_HPP_
#include <vector>
#include "Haplotype.hpp"
#include "Read.hpp"
#include "MLAlignment.hpp"
#include "ObservationModel.hpp"
using namespace std;

const double LOGTINY=-100.0;
const double EPS=1e-10;
// simple HMM inference algorithm for observation of single read given from a true underlying haplotype
class ObservationModelFB
{
public:
	ObservationModelFB() {};
	virtual ~ObservationModelFB();
	MLAlignment calcLikelihood();
	double getLogLikelihood() { calcLikelihood(); return ml.ll; };
	double* getMarginal(int readBase) { assert(readBase<readSize); return mar[readBase]; };
	void getObsVector(int b, double *vec) const;
	/*!
	 * @abstract Change haplotype used for likelihood computation, which
	 * can be useful in the EM algorithm
	 * @params newHap new haplotype
	 * @discussion newHap must have same length as previous haplotype
	 */
	void changeHaplotype(const Haplotype & newHap);
	void printMarginals();
	void computeMarginals();
	void computeXMarginals();
	void printAlignment(size_t hapScrPos);
	void printStatistics();
	vector<double> getOffHapLik() const { return likOffHap; };
	const MLAlignment &  getMLAlignment() const { return ml; }

protected:
	virtual void Init(const Haplotype & _hap, uint32_t hapStart);
	void forceOnHap();
	void setupReadObservationPotentials();
	virtual void setupTransitionProbs();
	virtual void initHMM();
	//void runHMM();
	virtual void allocateMemory();
	virtual void deleteMemory();
	virtual void computeBMidPrior(vector<double> & _prior, double mapQual);

	void passMessageOneInc(double *alpha_l, const double *alpha_l_1,  const double *obs_l_1);
	void passMessageOneDec(double *alpha_l, const double *alpha_l_1, const double *obs_l_1);
	virtual void passMessageTwoInc(double *beta_l, const double *beta_l_1,const double *obs_l_1);
	virtual void passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1);

	virtual void computeForwardMessages();
	void computeBackwardMessages();
	virtual void calcLikelihoodFromLastSlice();
	void printMarginalsInt(const vector<double*> & pot);
	bool hasErrors();
	bool _badValue(double v);

	// maximum number of bases that may be skipped by the sequencer
	// exponential decay rate of skip probability


	double logLikelihood, logLikelihoodNoPrior;
	MLAlignment ml;

	Haplotype hap;
	Read read;
	vector<double> prior;
	vector<double> logPTrans;
	vector<double> likOffHap;
	vector<double> priorOffHap;



	bool obsInitialized, memAllocated, HMMInitialized, HMMConsistent, likelihoodComputed, forwardDone, backwardDone, marginalsComputed;

	bool makeObsVector;

	// potentials are stored as log-values
	// observation potentials
	vector<double*> obs;

	// posterior marginals given _hap
	vector<double*> mar, xmar;
	// forward and backward messages
	vector<double*> alpha, beta;
	double *obsVector;

	vector<double> logProbError, logProbNoError;

	// structure of read-base variable

	// {LeftOfHaplotype, Hap_1, Hap_2, Hap_3, ..., Hap_length, RightOfHaplotype}

	// HMM internal variables

	int hapSize, readSize, ROState, bMid, hapStart, numT, numS;

	double logpLOgLO, logpFirstgLO;
	double logpInsgIns, logpInsgNoIns, logpNoInsgNoIns, logpNoInsgIns;

public:

	ObservationModelParameters params;
	ObservationModelFB(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & params);


};

class ObservationModelFBMax : public ObservationModelFB
{
public:
	void printMarginals() { throw string("Not possible for this model"); }
	void computeMarginals() { throw string("Not possible for this model"); }
	void computeXMarginals() { throw string("Not possible for this model"); }
	void reportVariants();
	void printAlignment(size_t hapScrPos);
	MLAlignment calcLikelihood();
	ObservationModelFBMax(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & params);
	ObservationModelFBMax() {};

	~ObservationModelFBMax() { deleteMemory(); };
	vector<int> getMapState() { runHMM(); return mapState; };
protected:
	void passMessageOneInc(double *alpha_l, const double *alpha_l_1,  const double *obs_l_1);
	void passMessageOneDec(double *alpha_l, const double *alpha_l_1, const double *obs_l_1);
	virtual void passMessageTwoInc(double *beta_l, const double *beta_l_1,const double *obs_l_1, int * bt_l);
	virtual void passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1, int * bt_l);
	void allocateMemory();
	void deleteMemory();
	void runHMM();
	void computeForwardMessages();
	void computeBackwardMessages();
	void calcLikelihoodFromLastSlice();
	void computeMAPState();
	inline void updateMax(double & destValue, int & destIdx, const double newValue, const int newIdx);
	vector<int *> btf, btb;
	vector<int> mapState; // MAP state for HMM
	//static double EPS=1e-10;
};

class ObservationModelFBMaxErr : public ObservationModelFBMax
{
public:
	ObservationModelFBMaxErr(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & params);

protected:
	void setupTransitionProbs();
	void passMessageTwoInc(double *beta_l, const double *beta_l_1,const double *obs_l_1, int * bt_l);
	void passMessageTwoDec(double *beta_l, const double *beta_l_1,const double *obs_l_1, int * bt_l);
};


#endif /*ObservationModelFB_HPP_*/
