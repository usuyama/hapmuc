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
 * VariantFile.hpp
 *
 *  Created on: Sep 9, 2009
 *      Author: caa
 */

#ifndef VARIANTFILE_HPP_
#define VARIANTFILE_HPP_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Utils.hpp"
#include "Variant.hpp"
using namespace std;

/*
 * class used as input for realignment
 */

class AlignedCandidates
{
public:
	AlignedCandidates()
	{
		tid="";
	}
	AlignedCandidates(const string & _tid, const vector<AlignedVariant> & _variants, int _leftPos, int _rightPos)
	{
		tid=_tid;
		variants=_variants;
		leftPos = _leftPos;
		rightPos = _rightPos;
		computePositions();
	}
	vector<AlignedVariant> variants;
	string tid;
	int centerPos, leftPos, rightPos;
    const void printAll() const {
        cout << " ******** " << endl;
		for (size_t x=0;x<variants.size();x++) {
			cout << "variants[x].pos: " << variants[x].getStartHap() << " " << variants[x].getString() << endl;
		}
        cout << " ******** " << endl;

    }
	const AlignedVariant * findVariant(int pos, int type, const string & str) const
	{
		//cout << " ******** " << endl;
		for (size_t x=0;x<variants.size();x++) {
			//cout << "pos: " << pos << " variants[x].pos: " << variants[x].getStartHap() << " " << variants[x].getString() << endl;
			if (variants[x].isEqual(pos, type, str)) return (const AlignedVariant *) &(variants[x]);
		}
	//	cout << " ******** " << endl;
		return NULL;
	}
private:
	void computePositions()
	{
		centerPos = leftPos+(rightPos-leftPos)/2;
	}
};


class VariantFile
{
public:
	class Candidates
	{
	public:
		Candidates()
		{
			tid="";
			pos=0;
			prior=-1.0;
		}
		Candidates(const string & _tid, uint32_t _pos, double _prior, const vector<Variant> & _variants, const vector<double> & _freqs)
		{
			tid=_tid;
			pos=_pos;
			prior=_prior;
			variants=_variants;
			freqs=_freqs;
		}
		vector<Variant> variants;
		vector<double> freqs;
		double prior;
		string tid;
		uint32_t pos;
	};


public:
	VariantFile(const string & fileName)
	{
		index=0;
		isOpen=false;
		open(fileName);
	}

	int open(const string & fileName)
	{
		fin.open(fileName.c_str());
		if (!fin.is_open()) throw string("Cannot open variant file ").append(fileName);
		isOpen=true;
		return 0;
	}

	bool eof() { if (isOpen) return fin.eof(); else return true; };

	Candidates getLine(bool isOneBased=false)
	{
		if (!isOpen) return empty;

		uint32_t pos;
		string tid;
		double prior=-1.0;

		string line;
		getline(fin, line);
		if (line.empty()) return empty;

		istringstream is(line);

		index++;

		if (!is.eof()) is >> tid; else return empty;
		if (!is.eof()) is >> pos; else return empty;

		// convert to zero-based coordinates
		if (isOneBased) pos--;

		// get variants from line
		vector<Variant> variants;
		vector<double> freqs;

		string col;

		try {
			while (!is.eof()) {
				is >> col;
				if ( col.size() && ( (col[0]!='-' && col[0] != '+' && col[0] != 'A' && col[0] != 'C' && col[0]!='G' && col[0]!='T' && col[0]!='R') ) ) break;
				Variant variant(col);
				if (variant.getSeq().size()!=0) variants.push_back(variant);
			}
		} catch (string err) {
			cerr <<  "Could not parse variants in line " << index << " in variants file." << endl;
			return empty;
		}

		if (col.find('#') != string::npos) return Candidates(tid, pos, prior, variants, freqs);

		prior=from_string<double>(prior, col, std::dec);

		bool error=false;
		while (!is.eof()) {
			string in;
			is >> in;

			if (in.find('#')!=string::npos) break;
			double freq;
			if (!from_string<double>(freq, in, std::dec)) { error=true; break; };
			freqs.push_back(freq);
		}

		if (error || (!error && freqs.size()>0 && freqs.size()!=variants.size())) {
			freqs.clear();
			cerr << "Could not parse all frequencies in line " << index << " in variants file." << endl;
		}

		if (variants.size()==0) {
			cerr << "Could not parse any variants in line: " << index << " SKIPPING." << endl;
			return empty;
		}

		return Candidates(tid, pos, prior, variants, freqs);
	}

	AlignedCandidates getLineVector(bool isOneBased=false)
	{
		if (!isOpen) return aligned_empty;

		uint32_t pos;
		int leftPos, rightPos;
		string tid;

		string line;
		getline(fin, line);
		if (line.empty()) return aligned_empty;

		istringstream is(line);

		index++;

		if (!is.eof()) is >> tid; else return aligned_empty;
		if (!is.eof()) {
			string str;
			is >> str;
			if (!from_string<int>(leftPos, str, std::dec)) throw string("Cannot read left boundary of region.");

		} else return aligned_empty;
		if (!is.eof()) {
			string str;
			is >> str;
			if (!from_string<int>(rightPos, str, std::dec)) throw string("Cannot read left boundary of region.");

		} else return aligned_empty;

		//cout << "leftPos: " << leftPos << " rightPos: " << rightPos << endl;

		// get variants from line
		vector<AlignedVariant> variants;
		vector<double> freqs;

		string col;

		try {
			while (!is.eof()) {
				string pvf_str;
				if (!is.eof()) is >> pvf_str;
				//cout << "pvf_str: " << pvf_str << endl;

				if (pvf_str.empty()) break;
				if (pvf_str[0]=='#' || pvf_str[0]=='%') break;

				vector<string> els;
				int lastpos=0;
				for (int x=0;x<int(pvf_str.size());x++) {
					if ((pvf_str[x]==';' || pvf_str[x]==',') && x-lastpos>0) {
						els.push_back(pvf_str.substr(lastpos,x-lastpos));
					//	cout << "els " << x << " : " << els[els.size()-1] << endl;
						lastpos=x+1;
					}
				}
				els.push_back(pvf_str.substr(lastpos, pvf_str.size()-lastpos));

				if (els.size()<2) {
					cerr << "Error reading line in variantfile!\n";
				} else {
					double freq=-1.0;
					bool addComb=false;
					if (!from_string<uint32_t>(pos, els[0], std::dec)) throw string("Cannot read position");
				// convert to zero-based coordinates
					if (isOneBased) pos--;

					string & col = els[1];
					if ( col.size()==0 || ( (col[0]!='-' && col[0] != '+' && col[0] != 'A' && col[0] != 'C' && col[0]!='G' && col[0]!='T' && col[0]!='R') ) ) throw string("Unrecognized variant");

					if (els.size()>2) {
						if (!from_string<double>(freq, els[2],std::dec)) throw string("Cannot read prior/frequency");
					}
					if (els.size()>3) {
						int addc;
						if (!from_string<int>(addc, els[3],std::dec)) throw string("Cannot add_combinatorial");
						if (addc) addComb = true;
					}
                    cout << col << pos << freq << addComb << endl;
					AlignedVariant variant(col,pos, freq, addComb);
					if (variant.getSeq().size()!=0) {
						variants.push_back(variant);
					}
				}
				// split into pos, var, col


			}
		} catch (string err) {
			cerr <<  "Could not parse variants in line " << index << " in variants file." << endl;
			cerr <<  "Error: " << err << endl;
			return aligned_empty;
		}

		if (variants.size()==0) {
			cerr << "Could not parse any variants in line: " << index << " SKIPPING." << endl;
			return aligned_empty;
		}

		return AlignedCandidates(tid, variants, leftPos, rightPos);
	}

	~VariantFile()
	{
		fin.close();
	}

protected:
	ifstream fin;
	bool isOpen;
	Candidates empty;
	AlignedCandidates aligned_empty;
	int index;
};


#endif /* VARIANTFILE_HPP_ */
