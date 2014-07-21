
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
 * Egt is a class to parse Illumina EGT cluster files.
 *
 * This is a port of the EGT Python class in zCall to C++.
 * There is no public documentation for the EGT binary file format.
 * (However, the zCall EGT class is claimed to originate from Illumina.)
 *
 * Repository for zCall:  https://github.com/wtsi-npg/zCall
 * 
 * Native EGT format stores coordinates as (R, Theta) polar coordinates
 * Python EGT class converts to Cartesian for storage in memory
 * This class stores as polar, but has method for Cartesian conversion
 * (More efficient than converting and storing as Cartesian, since both are 
needed for FCR output)
 */

#include <iostream>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string> 
#include "Egt.h"

using namespace std;



Egt::Egt(void)
{
  verbose = true; // TODO read verbose as an input parameter
  GENOTYPES_PER_SNP = 3; 
  PARAMS_PER_SNP = 12;
  NUMERIC_BYTES = 4;
  ENTRIES_IN_RECORD = 30;
  BYTES_IN_RECORD = NUMERIC_BYTES * ENTRIES_IN_RECORD;
  ENTRIES_TO_USE = 15;
}

void Egt::open(string filename)
{
  this->filename = filename;
  ifstream file;
 
  file.open(filename.c_str());
  if (!file) {
    cout << "Can't open file: " << filename << endl << flush;
    exit(1);
  }
  // read header data
  readHeader(file);
  readPreface(file);
  if (verbose) {
    printHeader();
    printPreface();
  }
  // read cluster data
  counts = new int[GENOTYPES_PER_SNP*snpTotal]; // nAA, nAB, nBB
  // expected cluster positions, in polar coordinates (R, Theta) for (AA,AB,BB)
  // Order of params is: devRAA, devRAB, devRBB, meanRAA, meanRAB, meanRBB, devThetaAA, devThetaAB, devThetaBB, meanThetaAA, meanThetaAB, meanThetaBB
  // Can't use a multidimensional array because we don't know snpTotal at compile time; instead use a pseudo-multidimensional array such that the (i,j)th value is (i*WIDTH + j)
  params = new float[PARAMS_PER_SNP*snpTotal];
  char *block = new char[BYTES_IN_RECORD];
  for (int i=0; i<snpTotal; i++) { 
    file.read(block, BYTES_IN_RECORD);
    int *ints = bytesToInts(block, 0, GENOTYPES_PER_SNP);
    float *floats = bytesToFloats(block, GENOTYPES_PER_SNP, ENTRIES_IN_RECORD);
    for (int j=0;j<GENOTYPES_PER_SNP;j++)
      counts[i*GENOTYPES_PER_SNP + j] = ints[j];
    for (int j=0;j<PARAMS_PER_SNP;j++)
      params[i*PARAMS_PER_SNP + j] = floats[j];
    delete floats;
    delete ints;
  }
  delete block;
  snpNames = readSNPNames(file);
  file.close();
}

void Egt::open(char *filename)
{
  string f = filename;
  open(f);
}


int* Egt::bytesToInts(char block[], int start, int end) {
  // convert a section of a byte array into ints
  // start, end indices refer to positions in the array of ints (not bytes)
  int *results = new int[end - start];
  numericConverter converter;
  for (int i=start;i<end;i++) {
    for (int j=0;j<NUMERIC_BYTES;j++) {
      converter.ncChar[j] = block[i*NUMERIC_BYTES + j];
    }
    results[i-start] = converter.ncInt;
  }
  return results;
}

float* Egt::bytesToFloats(char block[], int start, int end) {
  // convert a section of a byte array into floats
  // start, end indices refer to positions in the array of floats (not bytes)
  float *results = new float[end - start];
  numericConverter converter;
  for (int i=start;i<end;i++) {
    for (int j=0;j<NUMERIC_BYTES;j++) {
      converter.ncChar[j] = block[i*NUMERIC_BYTES + j];
    }
    results[i-start] = converter.ncFloat;
  }
  return results;
}


numericConverter Egt::getNextConverter(ifstream &file) {
  // convenience method to read the next few bytes into a union
  // can then use the union for numeric conversion
  char * buffer;
  buffer = new char[NUMERIC_BYTES];
  file.read(buffer, NUMERIC_BYTES);
  numericConverter converter;
  for (int i=0; i<NUMERIC_BYTES; i++) {
    converter.ncChar[i] = buffer[i];
  }
  delete buffer;
  return converter;
}

float Egt::readFloat(ifstream &file) {
  float result;
  numericConverter converter = getNextConverter(file);
  result = converter.ncFloat;
  return result;

}

void Egt::readHeader(ifstream &file) 
{
  // populate instance variables with header values
  file.seekg(0); // set read position to zero, if not already there
  fileVersion = readInteger(file);
  gcVersion = readString(file);
  clusterVersion = readString(file);
  callVersion = readString(file);
  normalizationVersion = readString(file);
  dateCreated = readString(file);
  mode = file.get();
  manifest = readString(file);
}

int Egt::readInteger(ifstream &file)
{
  int result;
  numericConverter converter = getNextConverter(file);
  result = converter.ncInt;
  return result; 
}

void Egt::readPreface(ifstream &file) {
  // read the 'preface' from the body of an EGT file
  // assumes file is positioned at start of the body
  dataVersion = readInteger(file);
  opa = readString(file);
  snpTotal = readInteger(file);
}

string* Egt::readSNPNames(ifstream &file) {
  // read SNP names from an EGT file
  // assumes file is positioned at end of cluster (mean, sd) data
  string *names = new string[snpTotal];
  int pos = file.tellg();
  file.seekg(pos + 13 * snpTotal); // skip SNP quality scores
  for (int i=0;i<snpTotal;i++) {
    // skip genotype scores
    // length of strings is unknown, so must read each one and discard it
    readString(file);
  }
  for (int i=0;i<snpTotal;i++) {
    names[i] = readString(file);
  }
  return names;
}

string Egt::readString(ifstream &file) {
  // EGT string format is as follows: 
  // - First byte is a *signed* char encoding the string length
  // - Subsequent bytes contain the string
  // Total bytes read is (length encoded in first byte)+1 -- at most 128
  char length = file.get(); // get a single byte
  if (length <= 0)
    throw("Illegal string length in EGT file");
  char * buffer;
  buffer = new char[length+1];
  file.read(buffer, length);
  string result = string(buffer);
  delete buffer;
  return result;
}

void Egt::printHeader() {
  // convenience method to print file header
  cout << "FILE_VERSION " << fileVersion << endl;
  cout << "GC_VERSION " << gcVersion << endl;
  cout << "CLUSTER_VERSION " << clusterVersion << endl;
  cout << "CALL_VERSION " << callVersion << endl;
  cout << "NORMALIZATION_VERSION " << normalizationVersion << endl;
  cout << "DATE_CREATED " << dateCreated << endl;
  cout << "MODE " << (int) mode << endl;
  cout << "MANIFEST " << manifest << endl;
}

void Egt::printPreface() {
  // convenience method to print preface of "main" data section
  cout << "DATA_VERSION " << dataVersion << endl;
  cout << "OPA " << opa << endl;
  cout << "TOTAL_SNPS " << snpTotal << endl;
}
