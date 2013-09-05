//
// QC.cpp
//
// Class to compute QC metrics on .sim files
// Metrics include: Sample magnitude (normalized by probe), sample xydiff
//
//
// Author: Iain Bancarz <ib5@sanger.ac.uk>
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


#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdlib.h>  
#include <string.h>

#include "QC.h"
#include "Sim.h"


using namespace std;

QC::QC(string simPath) {

  qcsim = new Sim();
  if (!qcsim->errorMsg.empty()) {
    cout << qcsim->errorMsg << endl;
    exit(1);
  }
  qcsim->open(simPath);
  cout << "Opened .sim file " << simPath << endl;

}


void QC::writeMagnitude(string outPath) {
  // compute normalized magnitudes by sample, write to given file
  FILE *outFile = fopen(outPath.c_str(), "w");
  qcsim->reset(); // return read position to first sample
  float *probeMagArray;
  probeMagArray = (float *) calloc(qcsim->numProbes, sizeof(float));
  magnitudeByProbe(probeMagArray);
  float *sampleMagArray;
  sampleMagArray = (float *) calloc(qcsim->numSamples, sizeof(float));
  char sampleNames[qcsim->numSamples][Sim::SAMPLE_NAME_SIZE+1];
  qcsim->reset(); 
  magnitudeBySample(sampleMagArray, probeMagArray, sampleNames);
  for (unsigned int i=0; i<qcsim->numSamples; i++) {
    // use fprintf to control number of decimal places
    fprintf(outFile, "%s\t%.6f\n", sampleNames[i], sampleMagArray[i]);
  }
  fclose(outFile);
}

void QC::writeXydiff(string outPath) {
  // compute XY intensity difference by sample, write to given file
  if (qcsim->numChannels!=2) {
    cerr << "Error: XY intensity difference is only defined for exactly "
      "two intensity channels." << endl;
    exit(1);
  }
  FILE *outFile = fopen(outPath.c_str(), "w");
  qcsim->reset(); // return read position to first sample
  float *xydArray;
  xydArray = (float *) calloc(qcsim->numSamples, sizeof(float));
  char sampleNames[qcsim->numSamples][Sim::SAMPLE_NAME_SIZE+1];
  xydiffBySample(xydArray, sampleNames);
  for (unsigned int i=0; i<qcsim->numSamples; i++) {
    // use fprintf to control number of decimal places
    fprintf(outFile, "%s\t%.6f\n", sampleNames[i], xydArray[i]);
  }
  fclose(outFile);
}

void QC::getNextMagnitudes(float magnitudes[], char *sampleName, Sim *sim) {
  // compute magnitudes for each probe from next sample in .sim input
  // can handle arbitrarily many intensity channels; also reads sample name
  vector<uint16_t> *intensity_int = new vector<uint16_t>;
  vector<float> *intensity_float = new vector<float>;
  if (sim->numberFormat == 0) {
    sim->getNextRecord(sampleName, intensity_float);
  } else {
    sim->getNextRecord(sampleName, intensity_int);
  }
  for (unsigned int i=0; i < sim->numProbes; i++) {
    float total = 0.0; // running total of squared intensities
    for (int j=0; j<sim->numChannels; j++) {
      int index = i*sim->numChannels + j;
      float signal;
      if (sim->numberFormat == 0) signal = intensity_float->at(index);
      else signal = intensity_int->at(index);
      total += signal * signal;
    }
    magnitudes[i] = sqrt(total);
  }
}

void QC::magnitudeByProbe(float magByProbe[]) {
  // iterate over samples; update running totals of magnitude by probe
  // then divide to find mean for each probe
  float *magnitudes;
  magnitudes = (float *) calloc(qcsim->numProbes, sizeof(float));
  char *sampleName; // placeholder; name used in magnitudeBySample
  sampleName = new char[qcsim->sampleNameSize];
  for(unsigned int i=0; i < qcsim->numSamples; i++) {
    getNextMagnitudes(magnitudes, sampleName, qcsim);
    for (unsigned int j=0; j < qcsim->numProbes; j++) {
      magByProbe[j] += magnitudes[j];
    }
  }
  for (unsigned int i=0; i < qcsim->numProbes; i++) {
    magByProbe[i] = magByProbe[i] / qcsim->numSamples;
  }
}

void QC::magnitudeBySample(float magBySample[], float magByProbe[], 
			   char sampleNames[][Sim::SAMPLE_NAME_SIZE+1]) {
  // find mean sample magnitude, normalized for each probe
  // also read sample names
  float *magnitudes;
  magnitudes = (float *) calloc(qcsim->numProbes, sizeof(float));
  for(unsigned int i=0; i < qcsim->numSamples; i++) {
    char *sampleName;
    sampleName = new char[qcsim->sampleNameSize];
    getNextMagnitudes(magnitudes, sampleName, qcsim);
    strcpy(sampleNames[i], sampleName);
    float mag = 0;
    for (unsigned int j=0; j < qcsim->numProbes; j++) {
      mag += magnitudes[j]/magByProbe[j];
    }
    magBySample[i] = mag / qcsim -> numProbes;
  }
}

void QC::xydiffBySample(float xydBySample[], 
			char sampleNames[][Sim::SAMPLE_NAME_SIZE+1]) {
  // find xydiff and read sample names
  for(unsigned int i=0; i < qcsim->numSamples; i++) {
    char *sampleName;
    sampleName = new char[qcsim->sampleNameSize];
    float xydTotal = 0.0; // running total of xy difference
    vector<uint16_t> *intensity_int = new vector<uint16_t>;
    vector<float> *intensity_float = new vector<float>;
    if (qcsim->numberFormat == 0) {
      qcsim->getNextRecord(sampleName, intensity_float);
    } else {
      qcsim->getNextRecord(sampleName, intensity_int);
    }
    strcpy(sampleNames[i], sampleName);
    for (unsigned int j=0; j<qcsim->numProbes; j++) {
      int index = j*qcsim->numChannels;
      float xyd;
      if (qcsim->numberFormat == 0) {
	xyd = intensity_float->at(index+1) - intensity_float->at(index);
      }
      else {
	xyd = intensity_int->at(index+1) - intensity_int->at(index);
      }
      xydTotal += xyd;
    }
    xydBySample[i] = xydTotal / qcsim->numProbes;
  }
}



