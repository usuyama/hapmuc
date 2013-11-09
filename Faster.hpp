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
 * Faster.hpp
 *
 *  Created on: Feb 24, 2009
 *      Author: caa
 */

#ifndef FASTER_HPP_
#define FASTER_HPP_
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
#include "Read.hpp"
#include "MLAlignment.hpp"
#include "ObservationModel.hpp"
using namespace std;

class ObservationModelS
{


protected:
	class UngappedAlignment
	{
	public:
		UngappedAlignment()
		{
			ll=-HUGE_VAL;
			relPos=-10000;
			numMismatch=10000;
		}
		UngappedAlignment(double _ll, int _relPos, int _numMismatch)
		{
			ll=_ll;
			relPos=_relPos;
			numMismatch=_numMismatch;
		}
		double ll;
		int relPos;
		int numMismatch;
	};
public:
	ObservationModelS() {};
	ObservationModelS(const Haplotype & _hap, const Read & r, uint32_t hapStart, const ObservationModelParameters & _params);
	virtual ~ObservationModelS();
	MLAlignment align(const HapHash & hash);

	//MLAlignment calcLikelihood();
	//double getLogLikelihood() { calcLikelihood(); return ml.ll; };
	// void changeHaplotype(const Haplotype & newHap);
	void printAlignment(size_t hapScrPos);
	void printStatistics();
	ObservationModelParameters params;

protected:
	void computeBMid();
	void setupReadLikelihoods();
	void Align();
	void reportVariants();
	inline void doTransition(int cs, int nr, const vector<int> & state, vector<double> & alpha, vector<double> & bt, const vector<double> & tr, const int & S);
	inline void doTransitionNF(int cs, int nr, const vector<int> & state, vector<double> & alpha, vector<double> & bt, const vector<double> & tr, const int & S);
	void SStateHMM(vector<int> & relPos);
	void AlignHash(const HapHash & hash);
	MLAlignment ml;
	vector<double> logMatch, cumLogMatch, logMismatch;
	vector<int> mapState;
	double llMatch; //log likelihood when all bases in the read match
	int bMid, hlen, rlen;
	double llOff, llOffHMQ, pOffFirst, pOffFirstHMQ;

	const Haplotype *hap_ptr;
	const Read *read_ptr;
	size_t hapStart;
	bool  likelihoodComputed, bMidError;
};


#endif /* FASTER_HPP_ */
