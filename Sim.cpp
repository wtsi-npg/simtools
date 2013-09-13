// Sim.cpp
//
// Author: Jennifer Liddle (js10)
//
// $Id: Sim.cpp 1354 2010-11-11 16:20:09Z js10 $
//
// Copyright (c) 2009 - 2010 Genome Research Ltd. 
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
//
#include "Sim.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <errno.h>
#include <stdio.h> // added for low-level file reading by ib5

using namespace std;

Sim::Sim(void) 
{
	version=0;
	inPath = "";
	filename = "";
	errorMsg="";
}

void Sim::openInput(string fname) 
{
        inPath = fname;
        char *f = new char[fname.length()+1];
	strcpy(f, fname.c_str());
        openLowLevel(f);
}

void Sim::openLowLevel(char *fname) {
  // open C-style file descriptor and read header data
  if (strcmp(fname, "-")==0) inFile = stdin;
  else inFile = fopen(fname, "r");
  char *magicChars = (char *) calloc(4, sizeof(char));
  int size_t;
  size_t = fread(magicChars, 1, 3, inFile);
  magic = string(magicChars);
  if (strcmp(magicChars, "sim") || size_t != 3) {
    throw("File " + string(fname) + " has invalid header [" 
	  + string(magic) + "]");
  }
  size_t = fread(&version, 1, 1, inFile);
  if (size_t != 1) throw("Error reading .sim file header version");
  size_t = fread(&sampleNameSize, 2, 1, inFile);
  if (size_t != 1) throw("Error reading .sim file header name size");
  size_t = fread(&numSamples, 4, 1, inFile);
  if (size_t != 1) throw("Error reading .sim file header samples");
  size_t = fread(&numProbes, 4, 1, inFile);
  if (size_t != 1) throw("Error reading .sim file header probes");
  size_t = fread(&numChannels, 1, 1, inFile);
  if (size_t != 1) throw("Error reading .sim file header channels");
  size_t = fread(&numberFormat, 1, 1, inFile);
  if (size_t != 1) throw("Error reading .sim file header numeric format");

  if (ferror(inFile)!=0) {
    throw("Error reading header from .sim file: [" + string(fname) + "]");
  }

  // calculate and store record length
  switch (numberFormat) {
  case FLOAT:	      numericBytes= 4;	break;
  case INTEGER:	      numericBytes = 2; break;
  case SCALED_INTEGER:  numericBytes = 2; break;
  default:	cerr << "Invalid number format " << numberFormat;
    exit(1);
  }
  sampleIntensityTotal = numProbes * numChannels;
  recordLength = numProbes * numChannels;
  recordLength *= numericBytes;
  recordLength += sampleNameSize;
}

void Sim::close(void) {
  // close input and output files (if open, and not equal to stdin or stdout)
  if (filename != "" && filename !="-") {
    fout.close();
  }
  if (inPath !="" && inPath!="-") { 
    if (ferror(inFile)) {
      cerr << "Input file is in error state!" << endl;
      exit(1);
    }
    int status;
    status = fclose(inFile);
    if (status!=0) { 
      cerr << "Failed to close input file" << endl;
      exit(1);
    }
  }
}

void Sim::openOutput(string fname) {
  filename = fname;
  _openOut(filename);

}

void Sim::reset(void)
{
  // return to starting point of .sim file (immediately after header)
  // use for multiple passes through .sim file in QC metrics
  if (filename == "-") {
    throw "Cannot reset file position in standard input!";
  }
  fseek(inFile, HEADER_LENGTH, 0);

}

void Sim::writeHeader(uint32_t _numSamples, uint32_t _numProbes, 
		      uint8_t _numChannels, uint8_t _numberFormat)
{
        if (filename=="") {
	  cerr << "Output not defined; need to call Sim::openOutput()" << endl;
	  exit(1);
	}
	version = VERSION;
	sampleNameSize = SAMPLE_NAME_SIZE;
	numSamples = _numSamples;
	numProbes = _numProbes;
	numChannels = _numChannels;
	numberFormat = _numberFormat;

	outfile->write("sim", 3);
	outfile->write((char*)&version, sizeof(version));
	outfile->write((char*)&sampleNameSize, sizeof(sampleNameSize));
	outfile->write((char*)&numSamples, sizeof(numSamples));
	outfile->write((char*)&numProbes, sizeof(numProbes));
	outfile->write((char*)&numChannels, sizeof(numChannels));
	outfile->write((char*)&numberFormat, sizeof(numberFormat));

	outfile->flush();
}

string Sim::dump(void)
{
	ostringstream s;
	s << "Number of SNPs:   " << numSamples << endl
	  << "Magic:            " << magic << endl
	  << "version:          " << version << endl
	  << "Sample Name Size: " << sampleNameSize << endl
	  << "Num Probes:       " << numProbes << endl
	  << "Num Channels:     " << numChannels << endl
	  << "Number Format:    " << numberFormat << endl
	  << endl;
	return s.str();
}


void Sim::__openout(ostream &f) 
{
	outfile = &f;
}

void Sim::_openOut(string filename)
{
  // open filestream for output
  if (filename == "-") {
    __openout(cout); 
  } else {
    fout.open(filename.c_str(),
	      ios::binary | ios::trunc | ios::in | ios::out); 
    __openout(fout); 
  }
  if (!outfile || !fout) {
    cerr << "Can't open " << filename << " for writing : " 
	 << strerror(errno) << endl;
  }
}

void Sim::getNextRecord(char *sampleName, uint16_t *intensity) {
  char *status = fgets(sampleName, sampleNameSize+1, inFile);
  if (status==NULL) throw("Error reading sample name from .sim file!");
  int items = fread(intensity, numericBytes, sampleIntensityTotal, inFile);
  if (items != sampleIntensityTotal || ferror(inFile)) {
    throw("Error reading intensities from .sim file!");
  }
}

void Sim::getNextRecord(char *sampleName, float *intensity) {
  char *status = fgets(sampleName, sampleNameSize+1, inFile);
  if (status==NULL) throw("Error reading sample name from .sim file!");
  int items = fread(intensity, numericBytes, sampleIntensityTotal, inFile);
  if (items != sampleIntensityTotal || ferror(inFile)) {
    throw("Error reading intensities from .sim file!");
  }
}

void Sim::write(void *buffer, int length)
{
        if (filename=="") {
	  cerr << "Output not defined; need to call Sim::openOutput()" << endl;
	  exit(1);
	}
	outfile->write((const char*)buffer,length);
}

