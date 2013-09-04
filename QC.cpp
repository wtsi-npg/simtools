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

#include "QC.h"
#include "Sim.h"


using namespace std;

void QC::hello(string simPath)
{
  cout << "Hello, and welcome to genotyping QC!" << endl;
  cout << ".sim path is: " << simPath << endl;
  Sim *sim = new Sim();
  sim->open(simPath);
  if (!sim->errorMsg.empty()) {
    cout << sim->errorMsg << endl;
    exit(1);
  }
  cout << "Magic:     " << sim->magic << endl;
  cout << "Version:   " << (int)sim->version << endl;
  cout << "Name Size: " << sim->sampleNameSize << endl;
  cout << "Samples:   " << sim->numSamples << endl;
  cout << "Probes:    " << sim->numProbes << endl;
  cout << "Channels:  " << (int)sim->numChannels << endl;
  cout << "Format:    " << (int)sim->numberFormat << endl;
  cout << "RecLength: " << (int)sim->recordLength << endl;
  cout << endl;
  
}

void QC::writeMagnitude(string simPath, string outPath) {

  cout << ".sim path is: " << simPath << endl;
  cout << "Magnitude path is: " << outPath << endl;

  Sim *sim = new Sim();
  sim->open(simPath);
  if (sim->numChannels != 2) throw("simtools can only handle SIM files with exactly 2 channels at present");

  float *magByProbe;
  magByProbe = (float *) calloc(sim->numProbes, sizeof(float));
  magnitudeByProbe(magByProbe, sim);
  cout << "Magnitude by probe found." << endl;
  for (int i=0; i<sim->numProbes; i++) cout << magByProbe[i] << " ";
  cout << endl;
  float *magBySample;
  magBySample = (float *) calloc(sim->numSamples, sizeof(float));
  sim->reset();
  magnitudeBySample(magBySample, magByProbe, sim);
  cout << "Magnitude by sample found." << endl;
  for (int i=0; i<sim->numSamples; i++) cout << magBySample[i] << " ";
  cout << endl;

}

void QC::getNextMagnitudes(float magnitudes[], char *sampleName, Sim *sim) {
  // compute magnitudes for each probe from next sample in .sim input
  // also reads sample name
  vector<uint16_t> *intensity_int = new vector<uint16_t>;
  vector<float> *intensity_float = new vector<float>;
  if (sim->numberFormat == 0) {
    sim->getNextRecord(sampleName, intensity_float);
  } else {
    sim->getNextRecord(sampleName, intensity_int);
  }
  float a, b;
  for (int i=0; i < sim->numProbes; i++) {
    int index = i*sim->numChannels;
    // hard-coded to find magnitude in first two channels only
    if (sim->numberFormat == 0) {
      a = intensity_float->at(index);
      b = intensity_float->at(index+1);
    } else {
      a = intensity_int->at(index);
      b = intensity_int->at(index+1);
    }
    magnitudes[i] = sqrt(a*a + b*b);
  }
}

void QC::magnitudeByProbe(float magByProbe[], Sim *sim) {
  // iterate over samples; update running totals of magnitude by probe
  // then divide to find mean for each probe

  float *magnitudes;
  magnitudes = (float *) calloc(sim->numProbes, sizeof(float));
  char *sampleName; // placeholder only; sample name needed elsewhere
  sampleName = new char[sim->sampleNameSize];
  for(unsigned int i=0; i < sim->numSamples; i++) {
    getNextMagnitudes(magnitudes, sampleName, sim);
    for (int j=0; j < sim->numProbes; j++) {
      magByProbe[j] += magnitudes[j];
    }
  }
  for (int i=0; i < sim->numProbes; i++) {
    magByProbe[i] = magByProbe[i] / sim->numSamples;
  }
}

void QC::magnitudeBySample(float magBySample[], float magByProbe[], 
			   Sim *sim) {
  // find mean sample magnitude, normalized for each probe
  // also read sample names
  float *magnitudes;
  magnitudes = (float *) calloc(sim->numProbes, sizeof(float));
  for(unsigned int i=0; i < sim->numSamples; i++) {
    char *sampleName;
    sampleName = new char[sim->sampleNameSize];
    getNextMagnitudes(magnitudes, sampleName, sim);
    float mag = 0;
    for (int j=0; j < sim->numProbes; j++) mag += magnitudes[j]/magByProbe[j];
    magBySample[i] = mag / sim -> numProbes;
  }
}






