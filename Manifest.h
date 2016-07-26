//
// Manifest.h
//
// Author: Jennifer Liddle (js10)
//
// $Id: Manifest.h 1354 2010-11-11 16:20:09Z js10 $
//
// Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, 
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 3. Neither the name of Genome Research Ltd nor the names of the 
// contributors may be used to endorse or promote products derived from 
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR 
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
// IN NO EVENT SHALL GENOME RESEARCH LTD. BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _MANIFEST_H
#define _MANIFEST_H

#include <map>
#include <string>
#include <vector>

#ifndef SWIG
#include <unordered_map>  // SGI extension to C++ STL standard
#endif

using namespace std;

class snpClass {

 public:

  snpClass(void) {

	index = -1;
	name = "?";
	chromosome = "?";
	position = -1;
	score = -1;
	snp[0] = '?';
	snp[1] = '?';
	iStrand = '?';
	cStrand = '?';
	normId = -1;
	BeadSetID = -1;
	converted = false;

  };

  string toString();
  string strandToString(char strand, bool converted);

  int index;
  string name;
  string chromosome;
  long position;
  float score;
  char snp [2];		// 'GT' or 'AG' for example. Defined to be A/B for TOP strand. 
  char iStrand;		// Original Illumina Strand: B (Bottom) or T (Top)
  char cStrand;		// Original Customer Strand: B (Bottom) or T (Top)
  int normId;		    // index into normalisation table
  int BeadSetID;      // Only available from "wide format" manifest files (.csv not .bpm.csv)

  bool converted; // has SNP been converted from original to ILMN top strand?
  // converted=true changes output of toString(), but not the iStrand and cStrand variables (which represent the original values)

  bool operator<(snpClass other) const {
    if (chromosome.compare(other.chromosome)) {
      return (chromosome.compare(other.chromosome) < 0);
    }
	if (position != other.position) {
      return (position < other.position);
    }
	return (name.compare(other.name) < 0);
  }
};


struct eqstr
{
  bool operator()(string s1, string s2) const
    {
        return s1.compare(s2) == 0;
    }
};
	   

class Manifest {
public:
	Manifest();
	void open(char *filename, bool wide = false);
	void open(string filename, bool wide = false);
	void open (string filename, string chromosome, bool wide = false);

	void order_by_position();
    void order_by_locus();

	string filename;

	vector<snpClass> snps; // If your code removes elements from this vector after
		 		// it's been populated, for now you must not use code
				// which uses the hash_map snpNames for lookup.
	map<int,int> normIdMap;
	void dump(void);

	snpClass* lookup_SNP_by_name (string snpname);  // Caller must delete object returned.
	// snpClass* lookup_SNP_by_position (long pos);    // Ditto.

	snpClass* findSNP(string snpname) { 
		int n = snp2idx((char*)snpname.c_str()); 
		return ((n == -1) ? NULL : &(snps[n])); 
	}

	string get_chromosome_for_SNP (string snpname);

	int snp2idx(char *snp) { return snpNames.find(snp) == snpNames.end() ? -1 : snpNames[snp]; }

	void exclude_cnvs();

	void write(string outpath); // write normalized .bpm.csv to file

 protected:
   void populate_hashmap();
	void convert (snpClass* snip, std::string input_snp); // Convert BOT SNPs to TOP format. 
	void test_convert();
    
	bool find_wide_columns(string& linestr, map<string,int>& dict); 

	int get_map_value (map<string, int>& mymap, const char* const treasure);

	// For fast lookup by SNP name. Key = name of SNP. Value = index in "snps" vector.
#ifdef SWIG
	map<string, int> snpNames;
#else
	unordered_map<string, int, hash<string>, eqstr> snpNames;
#endif

	string selectedChromosome; // Used if we are only concerned with a specific chromosome.

	bool EXCLUDE_CNVS;

	bool hasBeadSetID; // is extra BeadSetID column present?

};


#endif // _MANIFEST_H
