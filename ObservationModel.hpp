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
 * ObservationModel.hpp
 *
 *  Created on: Apr 2, 2009
 *      Author: caa
 */

#ifndef OBSERVATIONMODEL_HPP_
#define OBSERVATIONMODEL_HPP_
#include <iostream>
#include <string>

using namespace std;
class ObservationModelParameters
{
public:
		ObservationModelParameters() {
			modelType="probabilistic";
			setDefaultValues();
		}
		ObservationModelParameters(const string & modelType) {
			if (modelType=="threshold" || modelType=="probabilistic") this->modelType=modelType; else throw string("Model not supported.");
			setDefaultValues();
		}
		void setDefaultValues()
		{
			pError=1e-4;
			baseQualThreshold=0.995;
			fixedBaseQual=0.99;
			maxLengthIndel=10;
			mapQualThreshold=100.0;
			capMapQualFast=40.0;
			scaleErr=0.95;
			numE=3;
			pMut=1e-4;
			minOverlap=0;
			numIndels=1;
			indelDist="exponential";
			maxLengthDel=maxLengthIndel;
			pFirstgLO=0.01;
			checkBaseQualThreshold = 0.95;

			bMid=-1;
			forceReadOnHaplotype=false;
			mapUnmappedReads = false;

			padCover=5;
			maxMismatch=1;
			maxTryHash=5;
		}

		void print()
		{
			cout << "\tmodelType: " << modelType << endl;
			cout << "\tmaxLengthIndel: " << maxLengthIndel << " pError: " << pError << endl;
			cout << "\tbaseQualThreshold: " << baseQualThreshold << " fixedBaseQual: " << fixedBaseQual << endl;
			cout << "\tmapQualThreshold: " << mapQualThreshold << endl;
			cout << "\tcapMapQualFast: " << capMapQualFast << endl;
			cout << "\tminOverlap: " << minOverlap << endl;
			cout << "\tscaleError: " << scaleErr << endl;
			cout << "\tnumE: " << numE << endl;
			cout << "\tpMut: " << pMut << endl;
			cout << "\tnumIndels: " << numIndels << endl;
			cout << "\tindelDistribution: " << indelDist << endl;
			cout << "\tmaxLengthDel: " << maxLengthDel << " pError: " << pError << endl;
			cout << "\tpFirstgLO: " << pFirstgLO << endl;
			cout << "\tpadCover: " << padCover << endl;
			cout << "\tmaxMismatch: " << maxMismatch << endl;
			cout << "\tmaxTryHash: " << maxTryHash << endl;
			cout << "\tcheckBaseQualThreshold: " << checkBaseQualThreshold << endl;
			cout << "\tmapUnmappedReads: " << mapUnmappedReads << endl;
			//cout << "\tlogLikThreshod: " << logLikThreshold << endl;
		}
		double pError, baseQualThreshold, fixedBaseQual, mapQualThreshold, capMapQualFast, scaleErr, pMut;
		int maxLengthIndel, numE, minOverlap, numIndels, bMid;
		double checkBaseQualThreshold;

		string modelType, indelDist;
		int maxLengthDel,maxTryHash;
		bool forceReadOnHaplotype, mapUnmappedReads;
		double pFirstgLO;
		int padCover, maxMismatch;


};



#endif /* OBSERVATIONMODEL_HPP_ */
