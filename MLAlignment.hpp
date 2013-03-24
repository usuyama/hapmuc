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
 * MLAlignment.hpp
 *
 *  Created on: Apr 2, 2009
 *      Author: caa
 */

#ifndef MLALIGNMENT_HPP_
#define MLALIGNMENT_HPP_
#include <vector>
#include <string>
#include <map>
#include "Variant.hpp"
class MLAlignment
{
public:
	static const int INS=-1;
	static const int DEL=-2;
	static const int LO=-3;
	static const int RO=-4;
	MLAlignment()
	{
		relPos=-1;
		ll=0.0;
		llOn=0.0;
		llOff=0.0;
		offHap=false;
		offHapHMQ=false;
		numIndels=0;
		numMismatch=0;
		hr=-1;
		hl=-1;

	}

	int relPos; // relative position of read wrt haplotype



	int firstBase, lastBase; //first and last base of haplotype covered by the read
	map<int, AlignedVariant> indels, snps;

	map<int, bool> hapIndelCovered, hapSNPCovered; // indels and snps in the _haplotype_ covered by the read

	double ll; // loglikelihood
	double llOn, llOff; // without priors/mapping qualities taken into account
	bool offHap, offHapHMQ; // haplotype is mapped outside haplotype window with read-mapping quality and an artificial high-mapping-quality respectively
	int hl, hr; // left and rightmost base on haplotype covered by the read
	int numIndels, numMismatch;
	int nBQT, nmmBQT; // number of aligned bases and number mismatching above threshold
	double mLogBQ; // mean log BaseQuality

	int nMMLeft, nMMRight;

	string align;
	vector<int> hpos;
	operator double() const { return ll; };
	void print()
	{
		cout << "relPos: " << relPos << " offHap: " << (int) offHap << " ll: " << ll << endl;
	}
};

#endif /* MLALIGNMENT_HPP_ */
