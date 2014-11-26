
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

Fcr::Fcr() {
  // empty constructor
}

double Fcr::BAF(double theta, Egt egt, long snpIndex) {
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
  delete meanTheta;
  return baf;
}

// sanity check on manifest and gtc file
void Fcr::compareNumberOfSNPs(Manifest *manifest, Gtc *gtc) {
  if (manifest->snps.size() != gtc->xRawIntensity.size()) {
    ostringstream msg;
    msg << "Size mismatch: Manifest contains " << manifest->snps.size() 
        << " probes, but GTC " << gtc->filename << " contains " 
        << gtc->xRawIntensity.size() << " probes.";
    cerr << msg.str() << endl;
    throw msg.str();
  }
}

void Fcr::illuminaCoordinates(double x, double y, double &theta, double &r) {
  // convert (x,y) cartesian coordinates to Illumina coordinates (theta, r)
  // these are ***NOT*** standard polar coordinates!
  // the angle theta is rescaled s.t. pi/2 radians = 1 "Illumina angular unit"
  // r = x+y instead of r = sqrt(x**2 + y**2)
  theta = atan2(y, x)/(pi/2);
  r = x + y;

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
  header += "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tGC Score\tTheta\tR\tX\tY\tX Raw\tY Raw\tB Allele Freq\tLog R Ratio\n";
  return header;
}

double Fcr::logR(double theta, double r, Egt egt, long snpIndex) {
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
  delete meanR;
  delete meanTheta;
  return log2(r/rExpected);
}


void Fcr::write(Egt *egt, Manifest *manifest, ostream *outStream, 
              vector<string> infiles, vector<string> sampleNames) {
  // 'main' method to generate FCR and write to given output stream
 Gtc *gtc = new Gtc();
 string header = createHeader(manifest->filename, infiles.size(), 
                              manifest->snps.size());
 *outStream  << header;
 double epsilon = 1e-6;
 for (unsigned int i = 0; i < infiles.size(); i++) {
    // TODO is SCORES flag necessary?
    gtc->open(infiles[i], Gtc::XFORM | Gtc::INTENSITY | Gtc::SCORES);
    compareNumberOfSNPs(manifest, gtc);
    string sampleName;
    if (i < sampleNames.size()) sampleName = sampleNames[i];
    else sampleName = gtc->sampleName;
    for (unsigned int j = 0; j < manifest->snps.size(); j++) {
      string snpName = manifest->snps[j].name;
      double x_raw = gtc->xRawIntensity[j];
      double y_raw = gtc->yRawIntensity[j];
      float score = gtc->scores[j];
      double x_norm;
      double y_norm;
      unsigned int norm_id = manifest->normIdMap[manifest->snps[j].normId];
      char *alleles = manifest->snps[j].snp;
      gtc->normalizeIntensity(x_raw, y_raw, x_norm, y_norm, norm_id);
      // correction of negative intensities, for consistency with GenomeStudio
      if (x_norm < epsilon) { x_norm = 0.0; }
      if (y_norm < epsilon) { y_norm = 0.0; }
      char buffer[500] = { }; // initialize to null values
      if (x_raw == 0 || y_raw == 0){
        // zero intensity; set other fields to NaN
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
                sampleName.c_str(), alleles[0], alleles[1],
                score, theta, r, x_norm, y_norm,
                int(x_raw), int(y_raw), baf, logR);
      }
      *outStream << string(buffer);
    }
  }
 delete gtc;
} 
