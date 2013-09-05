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
  ofstream outFStream;
  ostream *outStream;
  qcsim->reset(); // return read position to first sample
  float *magByProbe;
  magByProbe = (float *) calloc(qcsim->numProbes, sizeof(float));
  magnitudeByProbe(magByProbe, qcsim);
  float *magBySample;
  magBySample = (float *) calloc(qcsim->numSamples, sizeof(float));
  char sampleNames[qcsim->numSamples][Sim::SAMPLE_NAME_SIZE+1];
  qcsim->reset(); 
  magnitudeBySample(magBySample, magByProbe, sampleNames, qcsim);
  if (outPath == "-") {
    outStream = &cout;
  } else {
    outFStream.open(outPath.c_str(), ios::binary | ios::trunc | ios::out);
    outStream = &outFStream;
  }
  for (unsigned int i=0; i<qcsim->numSamples; i++) {
    *outStream << sampleNames[i] << "\t" << magBySample[i] << endl;
  }
}

void QC::writeXydiff(string simPath, string outPath) {
  // compute XY intensity difference by sample, write to given file
  cout << "XYdiff not yet implemented!" << endl;
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

void QC::magnitudeByProbe(float magByProbe[], Sim *sim) {
  // iterate over samples; update running totals of magnitude by probe
  // then divide to find mean for each probe
  float *magnitudes;
  magnitudes = (float *) calloc(sim->numProbes, sizeof(float));
  char *sampleName; // placeholder; name used in magnitudeBySample
  sampleName = new char[sim->sampleNameSize];
  for(unsigned int i=0; i < sim->numSamples; i++) {
    getNextMagnitudes(magnitudes, sampleName, sim);
    for (unsigned int j=0; j < sim->numProbes; j++) {
      magByProbe[j] += magnitudes[j];
    }
  }
  for (unsigned int i=0; i < sim->numProbes; i++) {
    magByProbe[i] = magByProbe[i] / sim->numSamples;
  }
}

void QC::magnitudeBySample(float magBySample[], float magByProbe[], 
			   char sampleNames[][Sim::SAMPLE_NAME_SIZE+1], 
			   Sim *sim) {
  // find mean sample magnitude, normalized for each probe
  // also read sample names
  float *magnitudes;
  magnitudes = (float *) calloc(sim->numProbes, sizeof(float));
  for(unsigned int i=0; i < sim->numSamples; i++) {
    char *sampleName;
    sampleName = new char[sim->sampleNameSize];
    getNextMagnitudes(magnitudes, sampleName, sim);
    strcpy(sampleNames[i], sampleName);
    float mag = 0;
    for (unsigned int j=0; j < sim->numProbes; j++) {
      mag += magnitudes[j]/magByProbe[j];
    }
    magBySample[i] = mag / sim -> numProbes;
  }
}

