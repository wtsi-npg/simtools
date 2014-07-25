// Author: Iain Bancarz <ib5@sanger.ac.uk>
//
// Copyright (c) 2013 Genome Research Ltd. All rights reserved.
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

// Class to provide a high-level interface to simtools functionality
// Parse command-line options in simtools.cpp, input to methods in this class

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>

#include "Sim.h"
#include "Gtc.h"
#include "Egt.h"
#include "QC.h"
#include "Manifest.h"
#include "json/json.h"

using namespace std;

class SNPSorter {
  // compare SNP classes by position
 public:
  bool operator() (const snpClass &snp1, const snpClass &snp2);
};

class Commander {

 public:

  Commander();

  void compareNumberOfSNPs(Manifest *manifest, Gtc *gtc);
  void loadManifest(Manifest *manifest, string manfile);
  void parseInfile(string infile, vector<string> &sampleNames, vector<string> &infiles);
  //void readGTC(string gtcfile, Manifest manifest, vector<double> &intensities, vector<float> &scores, string &name);
  void normalizeIntensity(double x_raw, double y_raw, double &x_norm, double &y_norm, unsigned int norm_id, Gtc *gtc);
  void commandView(string infile, bool verbose);
  void commandCreate(string infile, string outfile, bool normalize, string manfile, bool verbose);
  void commandFCR(string infile, string outfile, string manfile, string egtfile, int start_pos, int end_pos, bool verbose);
  void commandIlluminus(string infile, string outfile, string manfile, int start_pos, int end_pos, bool verbose);
  void commandGenoSNP(string infile, string outfile, string manfile, int start_pos, int end_pos, bool verbose);
  void commandQC(string infile, string magnitude, string xydiff, bool verbose);




};
