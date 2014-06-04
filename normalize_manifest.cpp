//
// normalize_manifest
//
// Copyright (c) 2013 Genome Research Ltd.
//

// Author: Iain Bancarz <ib5@sanger.ac.uk>
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


// Executable to normalize a .bpm.csv manifest to the Illumina TOP strand
// Write normalized manifest to given output

// Example: /nfs/new_illumina_geno04/call/HumanOmniExpress-12v1_A.bpm.csv

#include <getopt.h>
#include <iostream>
#include "Manifest.h"

using namespace std;

static struct option long_options[] = {
                   {"infile", 1, 0, 0},
                   {"outfile", 1, 0, 0},
                   {"verbose", 0, 0, 0},
                   {"help", 0, 0, 0},
                   {0, 0, 0, 0}
               };


void showUsage(char *argv[]) {

  cout << "Usage:   " << argv[0] << " [options]" << endl << endl;

  cout << "Options: --infile <filename>    Input path to raw (unnormalized) .bpm.csv manifest" << endl;
  cout << "         --outfile <filename>   Output path for normalized .bpm.csv manifest" << endl;
  cout << "         --verbose              Show progress messages to STDERR" << endl;
  cout << "         --help                 Display this help text and exit" << endl;

  exit(0);

}


int main(int argc, char *argv[])
{
  string infile = ""; 
  string outfile = "";
  bool verbose = false;
  bool help = false;
  int option_index = -1;
  int c;

  // check command array
  if (argc < 2) showUsage(argv);

  // Get options
  while ( (c=getopt_long(argc,argv,"h?",long_options,&option_index)) != -1) {
    if (c == '?') showUsage(argv);	// unknown option
    if (option_index > -1) {
      string option = long_options[option_index].name;
      if (option == "infile") infile = optarg;
      if (option == "outfile") outfile = optarg;
      if (option == "verbose") verbose = true;
      if (option == "help") help = true;
    }
  }

  // TODO add an option to warn of any non-canonical SNPs (eg. T/C on TOP)

  if (help) {
    showUsage(argv);
  } else if (infile == "") {
    cerr << "Must specify an input file!" << endl;
    exit(1);
  } else if (outfile == "") {
    cerr << "Must specify an output file!" << endl;
    exit(1);
  }

  Manifest *manifest = new Manifest();
  manifest->open(infile);
  if (verbose) cerr << "Finished reading manifest: " << infile << endl;
  manifest->write(outfile);
  if (verbose) cerr << "Finished writing normalized manifest: " << outfile << endl;

}

