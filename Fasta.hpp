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
 * Fasta.hpp
 *
 *  Created on: May 27, 2009
 *      Author: caa
 */

#ifndef FASTA_HPP_
#define FASTA_HPP_

#include <string>
#include "bam.h"
#include "faidx.h"


namespace fasta { 
using namespace std;
class Fasta {
public:
	Fasta()
	{
		fai=NULL;
	}
	Fasta(const string & fileName)
	{
		fai=NULL;
		fai = fai_load(fileName.c_str());
		if (!fai) {
			throw string("Fasta: cannot open reference file.");
		}
	}

	string getSequence(const string & tid, int start, int end)
	{
		char *str;
		char *ref;
		str = (char*)calloc(strlen(tid.c_str()) + 30, 1);
		sprintf(str, "%s:%d-%d", tid.c_str(), start, end);
		int len;
		ref = fai_fetch(fai, str, &len);
		if (len==0) throw string("faidx error, len==0");

		string result(ref);
		transform(result.begin(), result.end(), result.begin(), ::toupper);
		free(str);
		free(ref);
		return result;
	}
	~Fasta()
	{
		if (fai) fai_destroy(fai);
	}

protected:
	faidx_t *fai;
};
}
#endif /* FASTA_HPP_ */
