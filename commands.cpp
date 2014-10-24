// Author: Iain Bancarz <ib5@sanger.ac.uk>, Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
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

#include "commands.h"
#include "Sim.h"
#include "Gtc.h"
#include "Egt.h"
#include "Fcr.h"
#include "QC.h"
#include "Manifest.h"
#include "json/json.h"

using namespace std;

// class defining a SNP ordering operator
// use class instance as argument to sort()
bool SNPSorter::operator() (const snpClass &snp1, const snpClass &snp2)
{
  if (snp1.chromosome.compare(snp2.chromosome))
    return (snp1.chromosome.compare(snp2.chromosome) < 0);
  if (snp1.position != snp2.position)
    return (snp1.position < snp2.position);
  return (snp1.name.compare(snp2.name) < 0);
}

// constructor
Commander::Commander() {

}

// sanity check on manifest and gtc file
void Commander::compareNumberOfSNPs(Manifest *manifest, Gtc *gtc) {
  if (manifest->snps.size() != gtc->xRawIntensity.size()) {
    ostringstream msg;
    msg << "Size mismatch: Manifest contains " << manifest->snps.size() 
        << " probes, but GTC " << gtc->filename << " contains " 
        << gtc->xRawIntensity.size() << " probes.";
    cerr << msg.str() << endl;
    throw msg.str();
  }
}


// convenience function to load manifest
void Commander::loadManifest(Manifest *manifest, string manfile)
{
  if (manfile == "") throw("No manifest file specified");
  manifest->open(manfile);
}

void Commander::normalizeIntensity(double x_raw, double y_raw,
                                   double &x_norm, double &y_norm,
                                   unsigned int norm_id, Gtc *gtc) {
  // This is the normalization calculation, according to Illumina
  XFormClass *XF = &(gtc->XForm[norm_id]);
  double tempx = x_raw - XF->xOffset;
  double tempy = y_raw - XF->yOffset;
  double cos_theta = cos(XF->theta);
  double sin_theta = sin(XF->theta);
  double tempx2 = cos_theta * tempx + sin_theta * tempy;
  double tempy2 = -sin_theta * tempx + cos_theta * tempy;
  double tempx3 = tempx2 - XF->shear * tempy2;
  double tempy3 = tempy2;
  x_norm = tempx3 / XF->xScale;
  y_norm = tempy3 / XF->yScale;
}
  
// Parse the infile, which is either an ascii list of GTC files, 
// or a JSON format file
// Return an array of filenames and (for a JSON file) a list of sample names
void Commander::parseInfile(string infile, vector<string> &sampleNames, vector<string> &infiles)
{
  Json::Value root;   // will contains the root value after parsing.
  Json::Reader reader;
  ifstream f;
  
  f.open(infile.c_str());
  if (infile.find(".json") != string::npos) {
    // Parse JSON file
    bool parsingSuccessful = reader.parse( f, root );
    if ( !parsingSuccessful ) throw("Could not parse json file "+infile);
    for ( unsigned int index = 0; index < root.size(); ++index ) {
      sampleNames.push_back(root[index]["uri"].asString());
      infiles.push_back(root[index]["result"].asString());
    }
  } else {
    // simple ascii text file
    string filename;
    while (f >> filename) infiles.push_back(filename);
  }
  f.close();
  // do simple validation
  if (infiles.size() == 0) throw("No GTC files are specified in the infile");
  Gtc *gtc = new Gtc();
  for (unsigned int i = 0; i < infiles.size(); i++) {
    gtc->open(infiles[i],0);
    if (gtc->errorMsg.length()) throw gtc->errorMsg;
  }
  delete gtc;
}

// View the header of a .sim file and (optionally) its contents
void Commander::commandView(string infile, bool verbose)
{
  Sim *sim = new Sim();
  
  cout << endl << "Commander::commandView" << endl;
  cout << endl << "Reading SIM file: " << infile << endl;
  sim->openInput(infile);
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
  
  char *sampleName = new char[sim->sampleNameSize+1];
  uint16_t *intensity_int = 
    (uint16_t *) calloc(sim->sampleIntensityTotal, sizeof(uint16_t));
  float *intensity_float = 
    (float *) calloc(sim->sampleIntensityTotal, sizeof(float));
  int i;
  for (unsigned int n = 0; n < sim->numSamples; n++) {
    if (sim->numberFormat == 0) sim->getNextRecord(sampleName,intensity_float);
    else                        sim->getNextRecord(sampleName,intensity_int);
    cout << sampleName << "\t: ";
    if (verbose) {	// dump intensities as well as sample names
			// there *must* be a better way of doing this...
      if (sim->numberFormat == 0) {
	for (i=0; i<sim->sampleIntensityTotal; i++) {
	  cout << intensity_float[i] << " ";
	}
      } else {
	for (i=0; i<sim->sampleIntensityTotal; i++) {
	  cout << intensity_int[i] << " ";
	}
      }
    }
    cout << endl;
  }
  sim->reportNonNumeric();
  free(intensity_int);
  free(intensity_float);
  delete sampleName;
  sim->close();
  delete sim;
}

//
// Create a SIM file from one or more GTC files
//
// infile      a file containing either a simple list of GTC files, or a list in JSON format
// outfile     the name of the SIM file to create, or '-' to write to stdout
// normalize   if true, normalize the intensities, else store the raw values in the SIM file
// manfile     the name of the manifest file
// verbose     boolean (default false)
//
// Note the the SIM file is written with the intensities sorted into position order, as given
// by the manifest file.
//
void Commander::commandCreate(string infile, string outfile, bool normalize, string manfile, bool verbose)
{
  vector<string> sampleNames;	// list of sample names from JSON input file
  vector<string> infiles;	// list of GTC files to process
  Sim *sim = new Sim();
  Gtc *gtc = new Gtc();
  Manifest *manifest = new Manifest();
  int numberFormat = normalize ? 0 : 1;
  
  //
  // First, get a list of GTC files. and possibly sample names
  //
  if (infile == "") throw("commandCreate(): infile not specified");
  parseInfile(infile,sampleNames,infiles);
  
  // We need a manifest file to sort the SNPs and to normalise the intensities (if required)
  loadManifest(manifest, manfile);
  // Sort the SNPs into position order
  sort(manifest->snps.begin(), manifest->snps.end(), SNPSorter());
  
  // Create the SIM file and write the header
  sim->openOutput(outfile);
  sim->writeHeader(infiles.size(), manifest->snps.size(), 2, numberFormat);
  
  // For each GTC file, write the sample name and intensities to the SIM file
  for (unsigned int n = 0; n < infiles.size(); n++) {
    gtc->open(infiles[n], Gtc::XFORM | Gtc::INTENSITY);
    compareNumberOfSNPs(manifest, gtc);
    char *buffer = new char[sim->sampleNameSize+1];
    memset(buffer,0,sim->sampleNameSize);
    // if we have a sample name from the json file, use it
    if (n < sampleNames.size()) { strcpy(buffer, sampleNames[n].c_str()); }
    else                        { strcpy(buffer,gtc->sampleName.c_str()); }
    sim->write(buffer, sim->sampleNameSize);
    if (verbose) {
      cerr << "Gtc file " 
	   << n+1
	   << " of " 
	   << infiles.size()
	   << "  File: "
	   << infiles[n]
	   << "  Sample: "
	   << buffer
	   << endl;
    }
    // Note that we write the intensities in SNP order, sorted by position
    for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
      double xn;
      double yn;
      int idx = snp->index - 1;   // index is zero based in arrays, but starts from 1 in the map file
      double x_raw = gtc->xRawIntensity[idx];
      double y_raw = gtc->yRawIntensity[idx];
      unsigned int norm_id = manifest->normIdMap[snp->normId];
      if (normalize) {
        normalizeIntensity(x_raw, y_raw, xn, yn, norm_id, gtc);
      } else {
	xn = gtc->xRawIntensity[idx];
	yn = gtc->yRawIntensity[idx];
      }
      if (numberFormat == 0) {
	float v;
	v = xn; sim->write(&v,sizeof(v));
	v = yn; sim->write(&v,sizeof(v));
      } else {
	uint16_t v;
	v = xn; sim->write(&v,sizeof(v));
	v = yn; sim->write(&v,sizeof(v));
      }
    }
    
  }
  sim->close();
  delete sim;
}

// write a Final Call Report (FCR) file
//
// FCR consists of header and body
// Fields for each line in body: (snp_name, sample_id, allele_A, allele_B, 
// score, chr, pos, theta, R, X_normalized, Y_normalized, X_raw, Y_raw, 
// BAF, logR)
// 

void Commander::commandFCR(string infile, string outfile, string manfile, string egtfile, int start_pos, int end_pos, bool verbose)
{
  vector<string> sampleNames;	// list of sample names from JSON input file
  vector<string> infiles;	// list of GTC files to process
  Gtc *gtc = new Gtc();
  Manifest *manifest = new Manifest();
  Egt *egt = new Egt();
  Fcr *fcr = new Fcr(); // Final Call Report generator
  ofstream outFStream;
  ostream *outStream;
  if (outfile == "-") {
    outStream = &cout;
  } else {
    outFStream.open(outfile.c_str(),ios::trunc | ios::out);
    outStream = &outFStream;
  }
  if (infile == "") throw("commandCreate(): infile not specified");
  parseInfile(infile, sampleNames, infiles);
  loadManifest(manifest, manfile);
  //sort(manifest->snps.begin(), manifest->snps.end(), SNPSorter());
  egt->open(egtfile);
  // now we have output stream, GTC paths, populated manifest and EGT
  // want to read intensities and scores from each GTC file
  // iterate over all (snp, sample) pairs
  // write output to a (gzipped?) FCR file
  // Fields in each FCR line: snp_name, sample_id, allele_A, allele_B, score, chr, pos, theta, R, X_normalized, Y_normalized, X_raw, Y_raw, BAF, logR
  if (end_pos == -1) end_pos = manifest->snps.size();// - 1;
  string header = fcr->createHeader(manifest->filename, infiles.size(), 
                                    manifest->snps.size());
  *outStream  << header;
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
      normalizeIntensity(x_raw, y_raw, x_norm, y_norm, norm_id, gtc);
      double theta;
      double r;
      fcr->illuminaCoordinates(x_norm, y_norm, theta, r);
      double logR = fcr->logR(theta, r, *egt, j);
      double baf = fcr->BAF(theta, *egt, j);
      // produce tab-delimited output
      char buffer[500] = { }; // initialize to null values
      string format = string("%s\t%s\t%c\t%c\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f")+
        string("\t%d\t%d\t%.4f\t%.4f\n");
      sprintf(buffer, format.c_str(), snpName.c_str(), 
              sampleName.c_str(), alleles[0], alleles[1],
              score, theta, r, x_norm, y_norm,
              int(x_raw), int(y_raw), baf, logR);
      *outStream << string(buffer);
    }
  }
  delete gtc;
  delete manifest;
  delete egt;
}


//
// Generate Illuminus output
//
// infile		is a filename or '-' for stdin
// outfile		is a filename or '-' for stdout
// manfile 		is the full path to the manifest file
// start_pos	is the Probe number (starting from 0) to start from
// end_pos		is the Probe number (from 0 to numProbes-1) to end at, or -1
// verbose		if true will display progress messages to stderr
//
void Commander::commandIlluminus(string infile, string outfile, string manfile, int start_pos, int end_pos, bool verbose)
{
  Sim *sim = new Sim();
  ofstream outFStream;
  ostream *outStream;
  char *sampleName;
  vector<vector<float> > SampleArray;
  Manifest *manifest = new Manifest();
  
  if (outfile == "-") {
    outStream = &cout;
  } else {
    outFStream.open(outfile.c_str(),ios::binary | ios::trunc | ios::out);
    outStream = &outFStream;
  }
  
  sim->openInput(infile);
  
  if (sim->numChannels != 2) throw("simtools can only handle SIM files with exactly 2 channels at present");
  uint16_t *intensity_int = 
    (uint16_t *) calloc(sim->sampleIntensityTotal, sizeof(uint16_t));
  float *intensity_float = 
    (float *) calloc(sim->sampleIntensityTotal, sizeof(float));
  sampleName = new char[sim->sampleNameSize+1];
  
  // We need a manifest file to sort the SNPs
  loadManifest(manifest, manfile);
  // Sort the SNPs into position order
  sort(manifest->snps.begin(), manifest->snps.end(), SNPSorter());
  
  if (end_pos == -1) end_pos = sim->numProbes - 1;
  
  // load the (relevant parts of) the SIM file
  if (verbose) cerr << "Reading SIM file " << infile << endl;
  *outStream << "SNP\tCoor\tAlleles";
  for(unsigned int n=0; n < sim->numSamples; n++) {
    vector<float> *s = new vector<float>; 
    if (!s) { cerr << "new s failed" << endl; exit(1); }
    if (sim->numberFormat == 0) sim->getNextRecord(sampleName, intensity_float);
    else                        sim->getNextRecord(sampleName, intensity_int);
    for (int i=start_pos; i <= end_pos; i++) {
      for (int c=0; c < sim->numChannels; c++) {
	float v;
	int k = i * sim->numChannels + c;
	if (sim->numberFormat==0) v = intensity_float[k];
	else                      v = intensity_int[k];
	s->push_back(v);
      }
    }
    SampleArray.push_back(*s);
    // Ooops! This is hardcoded for two channels. To Be Fixed. FIXME
    *outStream << "\t" << sampleName << "A\t" << sampleName << "B";
  }
  *outStream << endl;
  
    // Now write it out in Illuminus format
  if (verbose) cerr << "Writing Illuminus file " << outfile << endl;
  for (int n = start_pos; n <= end_pos; n++) {
    *outStream << manifest->snps[n].name << "\t" << manifest->snps[n].position << "\t" << manifest->snps[n].snp[0] << manifest->snps[n].snp[1];
    for (unsigned int i = 0; i < sim->numSamples; i++) {
      vector<float> s = SampleArray[i];
      for (unsigned int j=0; j < sim->numChannels; j++) {
	int k = (n - start_pos) * sim->numChannels + j;
	*outStream << '\t' << setw(7) << std::fixed << setprecision(3) << s[k];
      }
    }
    *outStream << endl;
  }
  if (verbose) sim->reportNonNumeric();
  free(intensity_int);
  free(intensity_float);
  sim->close();
  delete sim;
}


void Commander::commandGenoSNP(string infile, string outfile, string manfile, int start_pos, int end_pos, bool verbose)
{
  Sim *sim = new Sim();
  ofstream outFStream;
  ostream *outStream;
  
  outStream = &cout;
  if (outfile == "-") {
  } else {
    outFStream.open(outfile.c_str(),ios::binary | ios::trunc | ios::out);
    outStream = &outFStream;
  }
  
  sim->openInput(infile);
  
  if (end_pos == -1) end_pos = sim->numSamples - 1;
  
  char *sampleName = new char[sim->sampleNameSize+1];
  uint16_t *intensity = (uint16_t *) calloc(sim->sampleIntensityTotal, 
					    sizeof(uint16_t));
  for (int n=0; n <= end_pos ; n++) {
    sim->getNextRecord(sampleName, intensity);
    if (n < start_pos) continue;
    *outStream << sampleName << "\t" << sampleName;
    for (int i=0; i<sim->sampleIntensityTotal; i+=2) {
      *outStream << "\t" << std::fixed << setprecision(3) << intensity[i];
      *outStream << " " << std::fixed << setprecision(3) << intensity[i+1];
    }
    *outStream << endl;
  }
  if (verbose) sim->reportNonNumeric();
  free(sampleName);
  free(intensity);
  sim->close();
  delete sim;
}


void Commander::commandQC(string infile, string magnitude, string xydiff, bool verbose)
{
  if (infile == "-") {
    // QC requires multiple passes through the .sim input
    cerr << "Error: QC metrics require a .sim file, cannot accept "
      "standard input." << endl;
    exit(1);
  } else if (magnitude == "" && xydiff == "") {
    cerr << "Error: Must specify at least one of "
      "--magnitude, --xydiff for QC" << endl;
    exit(1);
  }
  QC *qc = new QC(infile, verbose);
  if (magnitude!="") {
    qc->writeMagnitude(magnitude, verbose);
  }
  if (xydiff!="") {
    qc->writeXydiff(xydiff, verbose);
  }
  qc->close();
  delete qc;
  
}


