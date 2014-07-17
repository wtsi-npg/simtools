
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


#include <iostream>
#include <vector>
#include <map>
//#include <cstdio>
#include <iostream>
#include <fstream>
#include <string> 
#include "Egt.h"

using namespace std;

Egt::Egt(void)
{
  cout << endl << "Class to represent Illumina EGT cluster files" << endl;

  NUMERIC_BYTES = 4;

}

/*
 * Want to port the EGT.py class in zCall to C++
 * 
 * EGT format: Binary file with header, body
 * Data types in binary file:
 * - Bute
 * - 4-byte integers (32-bit signed long int)
 * - 4-byte floats (long double?)
 * - String, encoded with the first byte denoting length of the subsequent string. (Max 256 characters.)
 */




void Egt::open(string filename)
{

  ifstream file;
 
  cout << "Opening file: " << filename << endl;

  file.open(filename.c_str());
  if (!file) {
    cout << "Can't open file: " << filename << endl << flush;
    exit(1);
  }
  // read header data
  // sanity check on EGT version, use to verify little-endianness
  long fileVersionLE = readInteger(file, true);
  if (fileVersionLE >=0 and fileVersionLE <= 1000) {
    little_endian = true;
    fileVersion = fileVersionLE;
  }  else {
    little_endian = false;
    file.seekg(0);
    fileVersion = readInteger(file, LITTLE_ENDIAN);
  }
  gcVersion = readString(file);
  file.close();
}

void Egt::open(char *filename)
{
  string f = filename;
  open(f);
}

int Egt::readInteger(ifstream &file, bool littleEndian)
{
  long result = 0;
  char * buffer;
  buffer = new char[NUMERIC_BYTES];
  file.read(buffer, NUMERIC_BYTES);
  if (littleEndian)
    for (int i = NUMERIC_BYTES-1; i>=0; i--) {
      result = (result << 8) + buffer[i];
    }
  else
     for (int i = 0; i<NUMERIC_BYTES; i++) {
       result = (result << 8) + buffer[i];
     }
  delete buffer;
  return result;
}

string Egt::readString(ifstream &file) {
  // EGT string format is as follows: 
  // - First byte is a *signed* char encoding the string length
  // - Subsequent bytes contain the string
  // Total bytes read is (length encoded in first byte)+1 -- at most 128
  char length = file.get(); // get a single byte
  char * buffer;
  buffer = new char[length+1];
  file.read(buffer, length);
  string result = string(buffer);
  delete buffer;
  return result;
}
