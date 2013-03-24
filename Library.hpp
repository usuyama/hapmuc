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
 * MyBam.hpp
 *
 *  Created on: Aug 27, 2009
 *      Author: caa
 */
#ifndef LIBRARY_HPP
#define LIBRARY_HPP
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <string>
#include "bam.h"
#include "Utils.hpp"
#include "StringHash.hpp"

using namespace std;

class Library
{
public:
	Library()
	{
	}
	Library(int type)
	{
		if (type == 0) {
			maxins = 2000;
			vector<double> counts = vector<double>(maxins, 1.0);
			calcProb(counts);
		} else throw string("Library type not recognized.");
	}
	Library(const vector<double> & counts)
	{
		// read library from histogram file
		if (counts.size() == 0) {
			cout << "HUH LIB" << endl;
		}
		calcProb(counts);
	}

	int getMaxInsertSize() const { return maxins; };
	int getModus() const { return modeInsertSize; };
	double getProb(int x) const  {
		if (x<0)  x = -x;
		if (x>=maxins) x = maxins-1;
		return probs[x];
	}
	double getNinetyFifthPctProb() const { return ninetyfifth_pct_prob; };
	void print() const 
	{
		for (size_t x=0;x<probs.size();x++) {
			cout << " " << probs[x];
			if (log(probs[x])>0) {
				cout << "LOG>0: " << x << endl;
			}
		}
		cout << endl;
	}

protected:
	int modeInsertSize;
	void calcProb(const vector<double> & counts)
	{
		int max_isize = 2000;
		int max_count = -1;
		// get insertsize at maximum

		//cout << "CALCPROB: " << endl;
		for (int s=0;s<int(counts.size());s++) {
			if (counts[s]>=max_count) {
				max_count = counts[s];
				max_isize = s;
			}
		}

		
		maxins = 25*max_isize;
		if (maxins>int(counts.size())) { maxins=int(counts.size()); }
		probs=vector<double>(maxins, 0.0);
		double z=0.0, max=-1.0;
		modeInsertSize=0;
		for (int d=0;d<maxins;d++) {
			probs[d] = counts[d];
			if (probs[d]>max) {
				max=probs[d];
				modeInsertSize=d;
			}
			z+=probs[d];
		}
		for (int d=0;d<maxins;d++) {
			probs[d] /=z;
			if (probs[d]<1e-10) probs[d]=1e-10;
		//	cout << "probs[" << d<<"]: " << probs[d] << endl;
		}

		sortProbs = probs;
		std::sort(sortProbs.begin(), sortProbs.end());
		double sum = 0.0;
		for (int x=sortProbs.size()-1;x>0;x--) {
			sum+=sortProbs[x];
			if (sum>0.95) {
				ninetyfifth_pct_prob = sortProbs[x];
				break;
			}
		}
		cout << "max: " << max/z << " ninetyfifth_pct_prob: " << ninetyfifth_pct_prob << endl;

	}
	vector<double> probs, sortProbs;
	int maxins;
	double ninetyfifth_pct_prob;
};

class LibraryCollection : public string_hash<Library>
{
public:
	LibraryCollection()
	{
		(*this)["single_end"]=Library(0);
	}

	void addFromFile(const string & fileName)
	{
		ifstream fin;
		fin.open(fileName.c_str());
		if (!fin.is_open()) throw string("Cannot open variant file ").append(fileName);

		int numLibs = 0;
		int numLines = 0;

		vector<double> counts;

		string libName;

		int prev=-1;
		int max_count = -1; // highest count
		int max_isize = -1; // insert size corresponding to highest count
		while (!fin.eof()) {
			string line;
			getline(fin, line);
			numLines++;
			if (line.empty()) {
				break;
			}
			istringstream is(line);
			string isize_str, count_str;
			int isize=-1;
			double count=-1;
			is >> isize_str;
			if (isize_str == "#LIB") {
	
				if (counts.size()>0 && !libName.empty()) {
					// next library!!!!!
					cout << "Storing library: " << libName << endl;


					LibraryCollection::const_iterator it = this->begin();
					it = this->find(libName);
					if (it != this->end()) {
						cerr << "Error: libName: " << libName << endl;
						cerr << "Duplicate library IDs." << endl;
						throw string("Library error");
					} else {
						(*this)[libName]=Library(counts);
						cout << "Number of inserts: " << counts.size() << endl;
						numLibs++;
					}
					counts.clear();
					prev=-1;
				}

				string label;
				is >> label;
				libName = label;
				if (label.empty()) {
					cerr << "Cannot read library name in line " << numLines << " of " << fileName << endl;
					throw string("Cannot read library name ");
				}


				goto nextline;

			} else {
				if (!from_string<int>(isize, isize_str, std::dec)) { cerr << "Error reading from library file" << endl; }
				is >> count_str; if (!from_string<double>(count, count_str, std::dec)) { cerr << "Error reading from library file" << endl; }

				if (isize!=prev+1) {
					cout << "isize: " << isize << " prev: " << prev << endl;
					cerr << "Library insert sizes must be consecutive" << endl;
					throw string("Library error.");
				}
				if (count<0) {
					cerr << "Library insert size count is negative.." << endl;
					throw string("Library error.");
				}

				counts.push_back(count);
				prev=isize;
			}
			nextline:
			line.clear();
		}

		// store last library

		LibraryCollection::const_iterator it = this->begin();
		it = this->find(libName);
		if (it != this->end()) {
			cerr << "Duplicate library IDs." << endl;
			throw string("Library error");
		} else {
			(*this)[libName]=Library(counts);
			numLibs++;
			cout << "Library " << libName << " loaded with " << counts.size() << " insert sizes." << endl;
		}

		if (numLibs==0) {
			cerr << "Could not find any libraries. Are the headers specified correctly?" << endl;
		}

		fin.close();


	}

	double getMaxInsertSize() const
	{
		double max = -HUGE_VAL;
		for (LibraryCollection::const_iterator it = this->begin(); it!= this->end(); it++) {
			if (it->second.getMaxInsertSize()>max) max = it->second.getMaxInsertSize();
		}
		return max;
	}
	//const bam_header_t *getBamHeader() const { return bamHeader; }
protected:
	//const bam_header_t *bamHeader;

};

#endif
