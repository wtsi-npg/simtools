
// Egt.cpp
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

/*
 * Fcr is a class to generate Final Call Report (FCR) files.
 * Output is similar to the FCR files produced by GenomeStudio.
 *
 */

#include <ctime>
#include <cstdio>
#include <string> 
#include "Fcr.h"

using namespace std;

Fcr::Fcr() {
  // placeholder for constructor

}

void Fcr::cartesianToPolar(double x, double y, double &theta, double &r) {
  // convert (x,y) cartesian coordinates to (r, theta) polar
  theta = atan2(y, x);
  r = sqrt(pow(y,2) + pow(x,2));

}

string Fcr::createHeader(string content, int samples, int snps) {
  // generate standard FCR header
  // content argument is typically the manifest name
  // includes data set summary, and column heads for main body
  string header = "[Header]\n";
  header += "GSGT Version\tsimtools\n";
  // create a time structure and add to output
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  char *buffer = new char[20];
  // use date format from Illumina FCR files
  strftime(buffer, 20, "%m/%d/%Y %I:%M %p", timeinfo);
  header += "Processing Date\t"+string(buffer)+"\n";
  delete buffer;
  header += "Content\t"+content+"\n";
  // to_string is not working because of compiler issues
  buffer = new char[50];  
  sprintf(buffer, "%d", snps);
  header += "Num SNPs\t"+string(buffer)+"\n";
  header += "Total SNPs\t"+string(buffer)+"\n";
  buffer = new char[50];  
  sprintf(buffer, "%d", samples);
  header += "Num Samples\t"+string(buffer)+"\n";
  header += "Total Samples\t"+string(buffer)+"\n";
  header += "File\t1 of 1\n"; // no split across files
  header += "[Data]\n";
  // now add column headers for man body
  header += "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tGC Score\tChr\tPosition\tTheta\tR\tX\tY\tX Raw\tY Raw\tB Allele Freq\tLog R Ratio\n";
  //cerr << header;
  return header;

}
