
// Fcr.cpp
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
 * Includes the logR and BAF statistics defined in:
 * Peiffer, Daniel A., et al. "High-resolution genomic profiling of chromosomal aberrations using Infinium whole-genome genotyping." Genome research 16.9 (2006): 1136-1148.
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <string>
#include "Egt.h"
#include "Fcr.h"
#include "Gtc.h"

using namespace std;

const static double pi = 3.141593;

FcrWriter::FcrWriter() {
  // empty constructor
}

double FcrWriter::BAF(double theta, Egt egt, long snpIndex) {
  // estimate the B allele frequency by interpolating between known clusters
  float *meanTheta = new float[egt.GENOTYPES_PER_SNP];
  egt.getMeanTheta(snpIndex, meanTheta);
  double baf;
  if (theta < meanTheta[0]) {
    baf = 0.0;
  } else if (theta > meanTheta[2]) {
    baf = 1.0;
  } else if (theta < meanTheta[1]) {
    baf = ((theta - meanTheta[0])/(meanTheta[1] - meanTheta[0]))*0.5;
  } else {
    baf = 0.5 + ((theta - meanTheta[1])/(meanTheta[2] - meanTheta[1]))*0.5;
  }
  delete [] meanTheta;
  return baf;
}

// sanity check on manifest and gtc file
void FcrWriter::compareNumberOfSNPs(Manifest *manifest, Gtc *gtc) {
  if (manifest->snps.size() != gtc->xRawIntensity.size()) {
    ostringstream msg;
    msg << "Size mismatch: Manifest contains " << manifest->snps.size() 
        << " probes, but GTC " << gtc->filename << " contains " 
        << gtc->xRawIntensity.size() << " probes.";
    cerr << msg.str() << endl;
    throw msg.str();
  }
}

void FcrWriter::illuminaCoordinates(double x, double y, double &theta, double &r) {
  // convert (x,y) cartesian coordinates to Illumina coordinates (theta, r)
  // these are ***NOT*** standard polar coordinates!
  // the angle theta is rescaled s.t. pi/2 radians = 1 "Illumina angular unit"
  // r = x+y instead of r = sqrt(x**2 + y**2)
  theta = atan2(y, x)/(pi/2);
  r = x + y;

}

string FcrWriter::createHeader(string content, int samples, int snps) {
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

  // use date format from Illumina FCR files
  char *buffer = new char[20];
  strftime(buffer, 20, "%m/%d/%Y %I:%M %p", timeinfo);
  header += "Processing Date\t"+string(buffer)+"\n";
  delete [] buffer;

  header += "Content\t"+content+"\n";
  // to_string is not working because of compiler issues

  buffer = new char[50];
  sprintf(buffer, "%d", snps);
  header += "Num SNPs\t"+string(buffer)+"\n";
  header += "Total SNPs\t"+string(buffer)+"\n";
  delete [] buffer;

  buffer = new char[50];
  sprintf(buffer, "%d", samples);
  header += "Num Samples\t"+string(buffer)+"\n";
  header += "Total Samples\t"+string(buffer)+"\n";
  delete [] buffer;

  header += "File\t1 of 1\n"; // no split across files
  header += "[Data]\n";

  // now add column headers for man body
  header += "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tGC Score\tTheta\tR\tX\tY\tX Raw\tY Raw\tB Allele Freq\tLog R Ratio\n";
  return header;
}

double FcrWriter::logR(double theta, double r, Egt egt, long snpIndex) {
  // calculate the LogR metric for given (theta, r) of sample
  // snpIndex is position in the manifest (starting from 0)
  // get (theta, R) for AA, AB, BB from EGT and interpolate
  // find intersection with (theta, R) of sample to get R_expected
  float *meanR = new float[egt.GENOTYPES_PER_SNP];
  float *meanTheta = new float[egt.GENOTYPES_PER_SNP];
  egt.getMeanR(snpIndex, meanR);
  egt.getMeanTheta(snpIndex, meanTheta);
  double rExpected = 0.0;
  for (int i=1; i<egt.GENOTYPES_PER_SNP; i++) {
    if (theta < meanTheta[i] or i+1 == egt.GENOTYPES_PER_SNP) {
      // m = gradient of interpolated line
      double m = (meanR[i] - meanR[i-1])/(meanTheta[i] - meanTheta[i-1]);
      rExpected = m * (theta - meanTheta[i-1]) + meanR[i-1];
      break;
    }
  }
  delete [] meanR;
  delete [] meanTheta;
  return log2(r/rExpected);
}

void FcrWriter::write(Egt *egt, Manifest *manifest, ostream *outStream,
              vector<string> infiles, vector<string> sampleNames) {
  // 'main' method to generate FCR and write to given output stream
 Gtc *gtc = new Gtc();
 string header = createHeader(manifest->filename, infiles.size(),
                              manifest->snps.size());
 *outStream  << header;
 double epsilon = 1e-6;
 // control determines which fields are read from GTC binary
 int control =  Gtc::XFORM | Gtc::INTENSITY | Gtc::SCORES | Gtc::BASECALLS;
 for (unsigned int i = 0; i < infiles.size(); i++) {
    gtc->open(infiles[i], control);
    compareNumberOfSNPs(manifest, gtc);
    string sampleName;
    if (i < sampleNames.size()) sampleName = sampleNames[i];
    else sampleName = gtc->sampleName;
    for (unsigned int j = 0; j < manifest->snps.size(); j++) {
      string snpName = manifest->snps[j].name;
      unsigned short x_raw = gtc->xRawIntensity[j];
      unsigned short y_raw = gtc->yRawIntensity[j];
      float score = gtc->scores[j];
      double x_norm;
      double y_norm;
      unsigned int norm_id = manifest->normIdMap[manifest->snps[j].normId];
      XFormClass xf = gtc->XForm[norm_id];
      xf.normalize(x_raw, y_raw, x_norm, y_norm);
      // correction of negative intensities, for consistency with GenomeStudio
      if (x_norm < epsilon) { x_norm = 0.0; }
      if (y_norm < epsilon) { y_norm = 0.0; }
      char buffer[500] = { }; // initialize to null values
      if (abs(x_raw) < epsilon || abs(y_raw) < epsilon ){
        // (effectively) zero intensity; set other fields to NaN
        string format = string("%s\t%s\t-\t-\tNaN\tNaN\tNaN\tNaN\tNaN")+
          string("\t%d\t%d\tNaN\tNaN\n");
        sprintf(buffer, format.c_str(), snpName.c_str(), sampleName.c_str(),
                int(x_raw), int(y_raw));
      } else {
        // output metrics to correct precision
        double theta;
        double r;
        this->illuminaCoordinates(x_norm, y_norm, theta, r);
        double logR = this->logR(theta, r, *egt, j);
        double baf = this->BAF(theta, *egt, j);
        string format = string("%s\t%s\t%c\t%c\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f")+
          string("\t%d\t%d\t%.4f\t%.4f\n");
        sprintf(buffer, format.c_str(), snpName.c_str(), 
                sampleName.c_str(), gtc->baseCalls[j].a, gtc->baseCalls[j].b,
                score, theta, r, x_norm, y_norm,
                int(x_raw), int(y_raw), baf, logR);
      }
      *outStream << string(buffer);
    }
  }
 // delete gtc;
}

FcrReader::FcrReader(string infile) {
  ifstream inStream;
  string line;
  bool body = false;
  inStream.open(infile.c_str());
  totalPairs = 0;
  timeStampKey = "Processing Date";
  fileKey = "File";
  unsigned int fieldsExpected = 13;
  vector<string> header_lines;
  // populate the list of header prefixes, used as hash keys
  headerKeys.push_back("[Header]");
  headerKeys.push_back("GSGT Version");
  headerKeys.push_back(timeStampKey);
  headerKeys.push_back("Content");
  headerKeys.push_back("Num SNPs");
  headerKeys.push_back("Total SNPs");
  headerKeys.push_back("Num Samples");
  headerKeys.push_back("Total Samples");
  headerKeys.push_back(fileKey);
  headerKeys.push_back("[Data]");
  while (getline(inStream, line)) {
    if (line.compare(0, 8, "SNP Name")==0) {
      // column titles are first line of body
      body = true;
      continue;
    }
    if (body) {
      totalPairs += 1;
      // if length of tokens != 13, raise error
      vector<string> tokens = splitByWhiteSpace(line); 
      if (tokens.size() != fieldsExpected) {
        cerr << "Wrong number of fields in FCR line: Expected " <<
          fieldsExpected << ", found " << tokens.size() << endl;
        throw 1;
      }
      snps.push_back(tokens[0]);
      samples.push_back(tokens[1]);
      alleles_a.push_back(tokens[2]);
      alleles_b.push_back(tokens[3]);      
      gcScore.push_back(atof(tokens[4].c_str()));
      theta.push_back(atof(tokens[5].c_str()));
      radius.push_back(atof(tokens[6].c_str()));
      x.push_back(atof(tokens[7].c_str()));
      y.push_back(atof(tokens[8].c_str()));
      x_raw.push_back(atoi(tokens[9].c_str()));
      y_raw.push_back(atoi(tokens[10].c_str()));
      logR.push_back(atof(tokens[11].c_str()));
      baf.push_back(atof(tokens[12].c_str()));
    } else {
      header_lines.push_back(line);
    }
  }
  if (body==false) {
    cerr << "Body of FCR file not found!" << endl;
    throw 1;
  }
  inStream.close();
  this->header = parseHeader(header_lines);
}

bool FcrReader::equivalent(FcrReader other, bool verbose) {
  // check for equality on data fields with another FcrReader object
  bool equal = true;
  double epsilon = 1e-5;
  if (totalPairs != other.totalPairs) {
    equal = false;
    if (verbose) {
      cerr << "Number of (snp, sample) pairs is not equal" << endl;
    }
  } else if (!equivalentHeaders(other)) {
    equal = false;
    if (verbose) {
      cerr << "FCR headers are not equivalent" << endl;
    }
  } else {
    for (int i=0; i<totalPairs; i++) {
      if (snps[i] != other.snps[i]) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal SNPs at position " << i << endl; 
        }
        break;
      } else if (samples[i] != other.samples[i]) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal sample names at position " << i << endl; 
        }
        break;
      } else if (alleles_a[i] != other.alleles_a[i]) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal A alleles at position " << i << endl; 
        }
        break;
      } else if (alleles_b[i] != other.alleles_b[i]) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal B alleles at position " << i << endl; 
        }
        break;
      } else if (abs(gcScore[i] - other.gcScore[i]) > epsilon) {
        cerr << "ABS: " << abs(gcScore[i] - other.gcScore[i]) << endl;
        equal = false;
        if (verbose) { 
          cerr << "Unequal GC scores at position " << i << endl; 
        }
        break;
      } else if (abs(theta[i] - other.theta[i]) > epsilon) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal theta at position " << i << endl; 
        }
        break;
      } else if (abs(radius[i] - other.radius[i]) > epsilon) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal radius at position " << i << endl; 
        }
        break;
      } else if (abs(x[i] - other.x[i]) > epsilon) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal normalised x intensity at position " << i << endl; 
        }
        break;
      } else if (abs(y[i] - other.y[i]) > epsilon) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal normalised y intensity at position " << i << endl; 
        }
        break;
      } else if (x_raw[i] != other.x_raw[i]) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal raw x intensity at position " << i << endl; 
        }
        break;
      } else if (y_raw[i] != other.y_raw[i]) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal raw y intensity at position " << i << endl; 
        }
        break;
      } else if (abs(logR[i] - other.logR[i]) > epsilon) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal logR value at position " << i << endl; 
        }
        break;
      } else if (abs(baf[i] - other.baf[i]) > epsilon) {
        equal = false;
        if (verbose) { 
          cerr << "Unequal B allele frequency at position " << i << endl; 
        }
        break;
      }
    }
  }
  return equal;
}

bool FcrReader::equivalentHeaders(FcrReader other, bool verbose) {
  map<string, string> myHead = this->header;
  map<string, string> otherHead = other.header;
  bool equivalent = true;
  for (unsigned int i=0; i<headerKeys.size(); i++) {
    string myVal = myHead[headerKeys[i]];
    string otherVal = otherHead[headerKeys[i]];
    if (headerKeys[i].compare(fileKey) == 0 || 
	headerKeys[i].compare(timeStampKey) == 0) {
      continue; // ignore the timestamp, and "File K of N" lines
    } else if (myVal.compare(otherVal)!=0) {
      equivalent = false;
      if (verbose) {
	cerr << "Differing values in FCR headers: " << myVal << ", " 
	     << otherVal << endl;
      }
      break;
    }
  }
  return equivalent;
}

map<string, string> FcrReader::parseHeader(vector<string> input) {
  // parse header fields
  vector<unsigned int> keyLengths(headerKeys.size(), 0);
  for (unsigned int i=0; i<headerKeys.size(); i++) {
    keyLengths[i] = headerKeys[i].size();
  }
  map<string, string> header;
  for (unsigned int i=0; i<input.size(); i++) {
    for (unsigned int j=0; j<headerKeys.size(); j++) {
      if (input[i].compare(0, keyLengths[j], headerKeys[j])==0) {
        // remove prefix string from the map value
	// also remove the following tab character (if any)
	// [Header], [Data] and Content have empty strings as values
	int start;
	if (input[i].size() == keyLengths[j]) {  start = keyLengths[j]; }
	else { start = keyLengths[j] + 1; } // remove the tab
        header[headerKeys[j]] = input[i].substr(start);
        break;
      }
    }
  }
  // the File line is optional; all others should have values
  if (header.size() < headerKeys.size() -1) {
    cerr << "Insufficient lines parsed in header: Expected minimum " <<
      headerKeys.size() -1 << ", found " << header.size() << endl;
    throw 1;
  }
  return header;
}

vector<string> FcrReader::splitByWhiteSpace(string str) {
  // split line into tokens by iterating over a stringstream
  string buffer;
  stringstream ss(str); // Insert the string into a stream
  vector<string> tokens; // Create vector to hold our words
  while (ss >> buffer) {
    tokens.push_back(buffer);
  }
  return tokens;
}
