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
 * GetCandidates.hpp
 *
 *  Created on: Aug 27, 2009
 *      Author: caa
 */

#ifndef GETCANDIDATES_HPP_
#define GETCANDIDATES_HPP_
#include <map>
#include <ext/hash_map>
#include <vector>
#include "MyBam.hpp"
#include "Fasta.hpp"
#include "Variant.hpp"
using __gnu_cxx::hash;
namespace std { using namespace __gnu_cxx; }
// generic class for generating candidates
class GetCandidates
{
public:
	GetCandidates() {};
	GetCandidates(const string & bamFile);

	vector<AlignedVariant> get(const string tid, uint32_t start, uint32_t end);
	void get(const string & outputFileName); // outputs directly to filename of the whole BAMfile
	void outputToFile(const string & fileName);
protected:
	map<string, vector<AlignedVariant> > candidates; // candidates for every chromosome/reference sequence
	MyBam bam;
	virtual ~GetCandidates() {};
};

class GetCandidatesFromCIGAR : public GetCandidates
{
public:
	class Params
	{
	public:
		Params() {
			alignWindow=100;
			refPad=10;
		}
		int alignWindow, refPad;
	} params;
protected:
	class CIGARindel
	{
	public:
		CIGARindel(const uint32_t _refpos, int _len, const string _seq)
		{
			refpos=_refpos;
			len=_len;
			seq=_seq;
		}
		bool operator<(const CIGARindel & c) const
		{
			if (refpos==c.refpos) { if (seq!=c.seq) return seq<c.seq; else return len<c.len; } else return refpos<c.refpos;
		}
		uint32_t refpos;
		int len;
		string seq;
	};
	typedef hash_map<int, map<CIGARindel, int> > HMap;

	class CFFData
	{
	public:
		 HMap hmap;
	};

public:
	GetCandidatesFromCIGAR();
	static int getIndelFromCIGARFetchFunc(const bam1_t *b, void *data);
	void getIndelFromCIGARRegion(const vector<MyBam *> & myBams, const string & tid, int start, int end, const string & outputFileName, fasta::Fasta & fa);
	void realignCandidateFile(const string & _varFile, bool isOneBased, const string & outputFileName, const string & fastaName);
	void get(const string & bamFile, const string & outputFileName);
	void get(const string & bamFile, const string & outputFileName, const string & fastaName);
	~GetCandidatesFromCIGAR();
protected:
	vector<AlignedVariant> alignCIGAR(const string & tid, const CIGARindel & id, fasta::Fasta & fa);
	static void getIndelFromCIGAR(const bam1_t *b, vector<CIGARindel> & indels);
	void outputIndels(const string & tid, const HMap & hmap, ofstream & ofile, fasta::Fasta & fa, int outputType);
	typedef hash_map<int, long int> InsertSizes;
	typedef string_hash<InsertSizes> LibInsertSize;
	typedef LibInsertSize::iterator LibIterator;
	typedef InsertSizes::iterator InsIterator;

	void outputLibraries(LibInsertSize & libInsertSize, const string & outputFile);

};

#endif /* GETCANDIDATES_HPP_ */
