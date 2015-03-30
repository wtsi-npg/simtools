//
// Fcr.h
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

#ifndef _FCR_H
#define _FCR_H

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "Egt.h"
#include "Gtc.h"
#include "Manifest.h"

using namespace std;

class FcrWriter {

 public:
  FcrWriter();
  double BAF(double theta, Egt egt, long snpIndex);
  void compareNumberOfSNPs(Manifest *manifest, Gtc *gtc);
  void illuminaCoordinates(double x, double y, double &theta, double &r);
  string createHeader(string content, int samples, int snps);
  double logR(double theta, double r, Egt egt, long snpIndex);
  void write(Egt *egt, Manifest *manifest, ostream *outStream, vector<string> infiles, vector<string> sampleNames);

};

class FcrReader {
 // Container for data in an FCR file
 // Intended only for running tests
 // "Real" FCR files are typically too big to slurp into memory
 public:
  int totalPairs; // number of SNP/sample pairs
  string timeStampKey;
  string fileKey;
  vector<string> headerKeys;
  map<string, string> header;
  FcrReader(string infile);
  vector<string> snps;
  vector<string> samples;
  vector<string> alleles_a;
  vector<string> alleles_b;
  vector<double> gcScore;
  vector<double> theta;
  vector<double> radius;
  vector<double> x;
  vector<double> y;
  vector<int> x_raw;
  vector<int> y_raw;
  vector<double> logR;
  vector<double> baf;
  bool equivalent(FcrReader other, bool verbose=true);

 private:
  bool equivalentHeaders(FcrReader other, bool verbose=true);
  map<string, string> parseHeader(vector<string> header);
  vector<string> splitByWhiteSpace(string str);


};

#endif	// _FCR_H
