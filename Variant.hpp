/*
 * Dindel
 * http://www.sanger.ac.uk/resources/software/dindel/
 *
 * Copyright 2010, Kees Albers
 * Released under the GPL Version 3.
 */
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
 * Variant.hpp
 *
 *  Created on: Aug 27, 2009
 *      Author: caa
 */

#ifndef VARIANT_HPP_
#define VARIANT_HPP_
#include <string>
#include <vector>

using namespace std;


class Variant
{
public:
	Variant() {};
	Variant(const string & _str) { initFromString(_str); }

	int size() const { return length; };
	typedef enum { INS, DEL, SNP, REF } Type;
	const string & getString() const { return str; };
	const string & getSeq() const { return seq; }
	Type getType() const { return type; }
	bool isIndel() const { if (type==INS || type==DEL) return true; else return false; };
	bool isSNP() const { if (type==SNP) return true; else return false; };
	bool isRef() const { if (type==REF) return true; else return false; };
protected:
	void initFromString(const string & str)
	{
		int ok=1;
		if (str.size()>1) {
			if (str[0]=='-') {
				// deletion
				length=int(str.size())-1;
				seq=str.substr(1, length);
				type=Variant::DEL;
			} else if (str[0]=='+') {
				length=(int(str.size())-1);
				seq=str.substr(1, length);
				type=Variant::INS;
			} else if (str.size()==4 && str[1]=='=' && str[2]=='>') {
					type=Variant::SNP;
					seq=str;
					length=1;
			} else if (str=="*REF") {
					type=Variant::REF;
					seq=string("*REF");
					length=1;
			} else ok=0;
		} else ok=0;
		if (!ok) { cout << "input string: " << str << endl; throw string("Unrecognized variant"); }
		this->str=str;
	}
	Type type;
	string seq;
	string str;
	int length;
};


class VariantInfo
{
public:
    string chr;
    int start, end;
    string ref, obs;
    string ref_count_tumor, obs_count_tumor;
    string missrate_tumor, strandrate_tumor;
    string ref_count_normal, obs_count_normal;
    string missrate_normal, strandrate_normal;
    string fisher_score;
    VariantInfo() {};
};


class AlignedVariant : public Variant
{
public:
	AlignedVariant() {};
    VariantInfo info;
	AlignedVariant(const string & _str, int _startHap, int _endHap, int _startRead, int _endRead)
	{
		initFromString(_str);
		startHap=_startHap;
		endHap=_endHap;
		startRead=_startRead;
		endRead=_endRead;

		leftFlankHap = startHap;
		leftFlankRead = startRead;

		rightFlankHap = endHap;
		rightFlankRead = endRead;


		freq = -1.0;
		addComb = false;
	}
	AlignedVariant(const string & _str, int canonicalPos, double _freq=-1.0, bool _addComb = false)
	{
		initFromString(_str);
		startHap = canonicalPos;
		if (type==DEL) {
			endHap = startHap+length-1;
		} else {
			endHap = startHap;
		}
		startRead=-1;
		endRead = -1;

		leftFlankHap = startHap;
		leftFlankRead = startRead;

		rightFlankHap = endHap;
		rightFlankRead = endRead;


		freq=_freq;
		addComb = _addComb;
	}
    
    string print() const {
        cout << getStartHap()<< getType() <<  getString() << endl;
    }

	bool isCovered(int pad, int firstBase, int lastBase) const
	{
		if (firstBase+pad<=startRead && lastBase-pad>=endRead) return true; else return false;
	}

	bool operator<(const AlignedVariant & v) const
	{
		if (startHap!=v.startHap) return startHap<v.startHap; else return this->getString()<v.getString();
	}
	bool isEqual(int pos, int type, const string & str) const {
		if (this->type == type && this->startHap == pos) {
			if (type == AlignedVariant::SNP) {
				if (str.substr(1,3)==this->str.substr(1,3)) return true; else return false;
			} else {
				if (type == AlignedVariant::INS) {
					if (this->getString() == str) return true; else return false;
				} else if (type == AlignedVariant::DEL) {
					if (this->getString().size()==str.size()) return true; else return false;
				}
			}
		} else return false;
		return false;
	}
    bool isEqual(const AlignedVariant & v) const {
        cout << v.getStartHap()<< v.getType() <<  v.getString() << endl;
        cout << getStartHap()<< getType() <<  getString() << endl;
        
		return isEqual(v.getStartHap(), v.getType(), v.getString());
	}
	int getStartRead() const { return startRead; };
	int getStartHap() const { return startHap; };
	int getEndHap() const { return endHap; };
	double getFreq() const  { return freq; };
	bool getAddComb() const { return addComb; };

	int getLeftFlankHap() const { return leftFlankHap; }
	int getRightFlankHap() const { return rightFlankHap; }
	int getLeftFlankRead() const { return leftFlankRead; }
	int getRightFlankRead() const { return rightFlankRead; }

	void setFlanking(int _leftFlankHap, int _rightFlankHap, int _leftFlankRead, int _rightFlankRead)
	{
		leftFlankRead  = _leftFlankRead;
		rightFlankRead = _rightFlankRead;

		leftFlankHap  = _leftFlankHap;
		rightFlankHap = _rightFlankHap;
	}

protected:
	int startHap, endHap;     // position of variant in the haplotype the read is aligned to.
	int startRead, endRead;   // position of variant in the read aligned to the haplotype
	int leftFlankHap, rightFlankHap; // position of left and right base flanking the indel in the _haplotype_ (ie the target sequence)
	int leftFlankRead, rightFlankRead; // position of left and right base flanking the indel in the _read_ (ie the sequence aligned to the target sequence)
	double freq;
	bool addComb;              // add combinatorially in generation of candidate haplotypes?
};



#endif /* VARIANT_HPP_ */
