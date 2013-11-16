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
 * OutputData.hpp
 *
 *  Created on: Sep 7, 2009
 *      Author: caa
 */

#ifndef OUTPUTDATA_HPP_
#define OUTPUTDATA_HPP_
#include <vector>
#include "StringHash.hpp"
#include <string>
#include <iostream>
#include "foreach.hpp"
using namespace std;
using namespace boost;

class OutputData
{
	friend class Line;
public:
	OutputData(ostream &_out) { out = &_out; numLines=0;};

	OutputData(ostream & _out, const vector<string> & _columnLabels)
	{
		out=&_out;
		BOOST_FOREACH(string label, _columnLabels) (*this)(label);
		numLines=0;
	}
	OutputData & operator() (const string & label)
	{
		HashIt it=labelToColumn.find(label);
		if (it!=labelToColumn.end()) throw string("Duplicate label ").append(label);
		labelToColumn[label]=int(labels.size());
		labels.push_back(label);
		return *this;
	}
	string headerString() const
	{
		stringstream out;
		if (labels.size()>0) {
			out << labels[0];
			for (size_t x=1;x<labels.size();x++) out << " " << labels[x];
		}
		return out.str();
	}
	template<class T> void outputLine(T x)
	{
		*out << x << endl;
	}
	class Line
	{
	public:
		Line(const OutputData  & od) : lineData(od.labelToColumn.size(),"NA")
		{
			labelToColumnPtr = & od.labelToColumn;

		}
		string get(const string & columnLabel) const
		{
			string_hash<int>::const_iterator it=labelToColumnPtr->find(columnLabel);
			if (it != labelToColumnPtr->end()) {
				return lineData[it->second];
			} else throw string("Column label ").append(columnLabel).append(" not found!");
		}
		template<class T> Line & set(const string & columnLabel, T x)
		{
			string_hash<int>::const_iterator it=labelToColumnPtr->find(columnLabel);
			if (it != labelToColumnPtr->end()) {
				stringstream os;
				os << x;
				lineData[it->second]=os.str();
			} else throw string("Column label ").append(columnLabel).append(" not found!");
			return *this;
		}
		string toString() const
		{
			stringstream out;
			if (lineData.size()>0) {
				out << lineData[0];
				for (size_t x=1;x<lineData.size();x++) out << "\t" << lineData[x];
			}
			return out.str();
		}
		vector<string> lineData;
	protected:
		const string_hash<int> *labelToColumnPtr;
	};
	void output(const OutputData::Line & line)
	{
		numLines++;
		*out << line.toString() << endl;
	}

	ostream *out;
protected:
	typedef string_hash<int>::iterator HashIt;
	string_hash<int> labelToColumn;

	vector<string> labels;
	int numLines;
};




#endif /* OUTPUTDATA_HPP_ */
