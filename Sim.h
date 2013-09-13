//
// Sim.h
//
// Header file for Sim.cpp
//
// Author: Jennifer Liddle (js10)
//
// $Id: Sim.h 1354 2010-11-11 16:20:09Z js10 $
//

// Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
//
// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
// 3. Neither the name of the Genome Research Ltd nor the names of its contributors 
// may be used to endorse or promote products derived from software without specific
// prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR WARRANTIES, 
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. EVENT SHALL GENOME RESEARCH LTD. BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _SIM_H
#define _SIM_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <stdint.h>

using namespace std;

class Sim {
public:
	// Number formats
	static const int FLOAT = 0;
	static const int INTEGER = 1;
	static const int SCALED_INTEGER = 2;

	static const int VERSION = 1;
	static const int SAMPLE_NAME_SIZE = 255;
	static const int HEADER_LENGTH = 16;
	
public:
	Sim();
	void openInput(string filename);
	void openLowLevel(char *f);
	void close(void);
	void reset(void);
	string dump(void);
	const char *dumpc(void) { return dump().c_str(); }
	void createFile(string filename);
	void writeHeader(uint32_t _numSamples, uint32_t _numProbes, uint8_t _numChannels=2, uint8_t _numberFormat=INTEGER);
	void write(void *buffer, int length);
	string errorMsg;
	string filename;
	string magic;			// expected to be "sim"
	uint8_t version;			// file version (expected to be 1)
	uint16_t sampleNameSize;		// sample name
	uint32_t numSamples;
	uint32_t numProbes;
	uint8_t numChannels;
	uint8_t numberFormat;
	int recordLength;		// calculated when file opened and header read
	int numericBytes; // record size of each number in file
	int sampleIntensityTotal; // number of intensities for each sample

	// These inline functions are for the use of SWIG and Perl
	const char *getFilename(void) { return filename.c_str(); }
	const char *getErrorMsg(void) { return errorMsg.c_str(); }
	int getVersion(void) { return version; }
	const char *getMagic(void) { return magic.c_str(); }

	void getNextRecord(char *sampleName, uint16_t *intensity);
	void getNextRecord(char *sampleName, float *intensity);


private:
	ostream *outfile;
	ofstream fout;
	string inPath;
	FILE *inFile; // low-level file access for greater speed
	map<string,long> sampleIndex;
	void __openout(ostream &f);
	void _openOut(string fname);
};
#endif	// _SIM_H

