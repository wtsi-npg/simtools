//
// Gtc.h
//
// Header file for Gtc.cpp
//
// Author: Jennifer Liddle (js10)
//
// $Id: Gtc.h 1354 2010-11-11 16:20:09Z js10 $
//

// Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
//
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
//

#ifndef _GTC_H
#define _GTC_H

#include <string>
#include <vector>

using namespace std;

class XFormClass {
public:
	XFormClass();
	int version;
	float xOffset;
	float yOffset;
	float xScale;
	float yScale;
	float shear;
	float theta;
};

class BaseCallClass {
public:
	BaseCallClass(){};
	char a;
	char b;
};

class Gtc {
public:
	static const int XFORM = 1;
	static const int CONTROL = 2;
	static const int INTENSITY = 4;
	static const int GENOTYPES = 8;
	static const int BASECALLS = 16;
	static const int SCORES = 32;
	static const int ALL = 255;
public:
	Gtc();
	void open(string filename, int control=ALL);
	void open(char *f, int control=ALL);
	string dump(void);
	const char *dumpc(void) { return dump().c_str(); }
	double passRate(double cutoff);
	double correctedPassRate(double cutoff);

        void normalizeIntensity(double x_raw, double y_raw, double &x_norm, double &y_norm, unsigned int norm_id);

	string errorMsg;
	string filename;
	int version;			// file version (expected to be 3
	int numSnps;			// number of SNPs
	string sampleName;		// sample name
	string samplePlate;		// sample plate
	string sampleWell;		// sample well
	string clusterFile;		// cluster file name
	string manifest;		// name of SNP manifest
	string imagingDate;
	string autocallDate;
	string autocallVersion;
	string scannerName;
	int pmtGreen;
	int pmtRed;
	string scannerVersion;
	string imagingUser;
	vector<XFormClass> XForm;
	vector<unsigned short> xRawControl;
	vector<unsigned short> yRawControl;
	vector<unsigned short> xRawIntensity;
	vector<unsigned short> yRawIntensity;
	vector<char> genotypes;
	vector<BaseCallClass> baseCalls;
	vector<float> scores;

	// These inline functions are for the use of SWIG and Perl
	const char *getFilename(void) { return filename.c_str(); }
	const char *getErrorMsg(void) { return errorMsg.c_str(); }
	int getVersion(void) { return version; }
	int getNumSnps(void) { return numSnps; }
	const char * getSampleName(void) { return sampleName.c_str(); }
	const char * getSamplePlate(void) { return samplePlate.c_str(); }
	const char * getSampleWell(void) { return sampleWell.c_str(); }
	const char * getClusterFile(void) { return clusterFile.c_str(); }
	const char * getManifest(void) { return manifest.c_str(); }
	const char * getImagingDate(void) { return imagingDate.c_str(); }
	const char * getAutocallDate(void) { return autocallDate.c_str(); }
	const char * getAutocallVersion(void) { return autocallVersion.c_str(); }
	const char * getScannerName(void) { return scannerName.c_str(); }
	const char * getScannerVersion(void) { return scannerVersion.c_str(); }
	const char * getImagingUser(void) { return imagingUser.c_str(); }

	const char * getBaseCall(int n) { char *s = new char[3]; s[0]=baseCalls[n].a; s[1]=baseCalls[n].b; s[2]=0; return s; }
	float getScore(int n) { return scores[n]; }
	unsigned int getScoresSize(void) { return scores.size(); }
	unsigned short getXRawIntensity(int n) { return xRawIntensity[n]; }
	unsigned short getYRawIntensity(int n) { return yRawIntensity[n]; }

private:
	string readString(ifstream &f, int offset);
	string _readString(ifstream &f);
	void readXForm(ifstream &f, int offset);
	void readShortArray(ifstream &file, int offset, vector<unsigned short> &array);
	void readByteArray(ifstream &file, int offset, vector<char> &array);
	void readBaseCallArray(ifstream &file, int offset, vector<BaseCallClass> &array);
	void readFloatArray(ifstream &file, int offset, vector<float> &array);
};
#endif	// _GTC_H

