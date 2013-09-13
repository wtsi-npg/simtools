//
// QC.h
//
// Header file for QC.cpp
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

#ifndef _QC_H
#define _QC_H

#include <cmath>
#include <ctime>
#include <iostream>
#include <stdlib.h>  

#include "Sim.h"

using namespace std;

class QC {
 public:  
  static const int VERBOSE_FREQ = 1000; // frequency of verbose output
  static const int TIME_BUFFER = 100; // max size in bytes of timestamp string

  QC(string simPath, bool verbose);
  void close(void);
  void writeMagnitude(string outPath, bool verbose);
  void writeXydiff(string outPath, bool verbose);

 private:
  Sim *qcsim;
  uint16_t *intensity_int_array;
  float *intensity_float_array;

  void getNextMagnitudes(float magnitudes[], char* sampleName, Sim *sim);
  void magnitudeByProbe(float magByProbe[], bool verbose);
  void magnitudeBySample(float magBySample[], float magByProbe[], 
			 char sampleNames[][Sim::SAMPLE_NAME_SIZE+1],
			 bool verbose);
  void xydiffBySample(float xydBySample[], 
		      char sampleNames[][Sim::SAMPLE_NAME_SIZE+1]);
  void timeText(char *buffer);

};

#endif	// _QC_H
