//
// Egt.h
//
// Author: Iain Bancarz <ib5@sanger.ac.uk>
//
// Copyright (c) 2014 Genome Research Ltd.
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

#ifndef _EGT_H
#define _EGT_H

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

union numericConverter {
  float ncFloat;
  int ncInt;
  char ncChar[sizeof(int)];
};

class Egt {

 public:
  Egt(bool verbose=false);
  //~Egt();
  void getClusters(long index, float params[]);
  void getMeanR(long index, float means[]);
  void getMeanTheta(long index, float means[]);
  void open(char *filename);
  void open(string filename);
  void printHeader();
  void printPreface();
  string filename;
  bool verbose;
  // constants
  int NUMERIC_BYTES;
  int ENTRIES_IN_RECORD;
  int BYTES_IN_RECORD;
  int GENOTYPES_PER_SNP;
  int ENTRIES_TO_USE;
  int PARAMS_PER_SNP;
  // EGT header fields
  long fileVersion;
  string gcVersion;
  string clusterVersion;
  string callVersion;
  string normalizationVersion;
  string dateCreated;
  char mode;
  string manifest;
  // EGT 'file preface' fields
  long dataVersion;
  string opa;
  long snpTotal;
  // arrays for numerical data
  int *counts;
  float *params;
  // array for SNP names
  string *snpNames;

private:
  int* bytesToInts(char block[], int start, int end);
  float* bytesToFloats(char block[], int start, int end);
  numericConverter getNextConverter(ifstream &file);
  void readHeader(ifstream &file);
  int readInteger(ifstream &file);
  float readFloat(ifstream &file);
  void readPreface(ifstream &file);
  void readSNPNames(ifstream &file, string names[]);
  string readString(ifstream &file);

};

#endif	// _EGT_H
