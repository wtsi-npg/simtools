//
// Gtc.cpp
//
// Author: Jennifer Liddle (js10)
//
// $Id: Gtc.cpp 1354 2010-11-11 16:20:09Z js10 $
//
// Copyright (c) 2009 - 2010 Genome Research Ltd. 
//
// Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
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
#include "Gtc.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <iomanip>  
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

XFormClass::XFormClass(int version, float xOffset, float yOffset, 
                       float xScale, float yScale, float shear, float theta)
{
  this->version = version;
  this->xOffset = xOffset;
  this->yOffset = yOffset;
  this->xScale = xScale;
  this->yScale = yScale;
  this->shear = shear;
  this->theta = theta;  
}

void XFormClass::normalize(unsigned short x_raw, unsigned short y_raw, 
                           double &x_norm, double &y_norm)
{
  // This is the intensity normalization calculation, according to Illumina
  double tempx = x_raw - xOffset;
  double tempy = y_raw - yOffset;
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double tempx2 = cos_theta * tempx + sin_theta * tempy;
  double tempy2 = -sin_theta * tempx + cos_theta * tempy;
  double tempx3 = tempx2 - shear * tempy2;
  double tempy3 = tempy2;
  x_norm = tempx3 / xScale;
  y_norm = tempy3 / yScale;
}

string XFormClass::toString(void)
{
  int digits = 6;
  stringstream sstream;
  sstream << "Version: " << version << endl;
  sstream << setprecision(digits);
  sstream << "XOffset: " << xOffset << endl;
  sstream << "YOffset: " << yOffset << endl;
  sstream << "XScale: " << xScale << endl;
  sstream << "YScale: " << yScale << endl;
  sstream << "Shear: " << shear << endl;
  sstream << "Theta: " << theta << endl;  
  return sstream.str();
}

Gtc::Gtc(void) 
{
	version=0;
	filename = "";
}

void Gtc::open(char *f, int control)
{
	string s;
	s = f;
	open(s,control);
}

void Gtc::open(string filename, int control) 
{
	char buff[256];
	char c;
	int nEntries;
	ifstream file;

	errorMsg = "";

	// Clear arrays (to help prevent memory leaks)
	XForm.clear();
	xRawControl.clear();
	yRawControl.clear();
	xRawIntensity.clear();
	yRawIntensity.clear();
	genotypes.clear();
	baseCalls.clear();
	scores.clear();

	this->filename = "";
	file.open(filename.c_str(),ios::binary);
	if (!file) {
		errorMsg = "Can't open file " + filename;
		return;
	}
	file.get(buff,4);
	if (strcmp(buff,"gtc")) {
		errorMsg = "File " + filename + " has invalid header";
		return;
	}

	this->filename = filename;
	file.get(c);
	version = (int)c;

/*
	if (sizeof(int) != 4) {
		cout << "Errr...an integer seems to be " << sizeof(int) << " bytes" << endl;
	}
	if (sizeof(short) != 2) {
		cout << "Errr...a short seems to be " << sizeof(short) << " bytes" << endl;
	}
	if (sizeof(float) != 4) {
		cout << "Errr...a float seems to be " << sizeof(float) << " bytes" << endl;
	}
*/
	file.read((char*)&nEntries,4);

	for (int n=0; n<nEntries; n++) {
		short id;
		unsigned int offset;
		file.read((char*)&id,2);
		file.read((char*)&offset,4);
		switch (id) {
			case 1:		// number of SNPs
					numSnps = offset;
					break;

			case 10:	sampleName = readString(file,offset);	break;
			case 11:	samplePlate = readString(file,offset);	break;
			case 12:	sampleWell = readString(file,offset);	break;
			case 100:	clusterFile = readString(file,offset);	break;
			case 101:	manifest = readString(file,offset);		break;
			case 200:	imagingDate = readString(file,offset);	break;
			case 201:	autocallDate = readString(file,offset);	break;
			case 300:	autocallVersion = readString(file,offset);	break;
			case 400:	if (control & XFORM) { readXForm(file,offset); } break;
			case 500:	if (control & CONTROL) { readShortArray(file,offset,xRawControl); } break;
			case 501:	if (control & CONTROL) { readShortArray(file,offset,yRawControl); } break;
			case 1000:	if (control & INTENSITY) { readShortArray(file,offset,xRawIntensity); } break;
			case 1001:	if (control & INTENSITY) { readShortArray(file,offset,yRawIntensity); } break;
			case 1002:	if (control & GENOTYPES) { readByteArray(file,offset,genotypes); } break;
			case 1003:	if (control & BASECALLS) { readBaseCallArray(file,offset,baseCalls); } break;
			case 1004:	if (control & SCORES) { readFloatArray(file,offset,scores); } break;
			case 1005:	// scanner data
			            ios::pos_type pos = file.tellg();
			            file.seekg(offset);
						scannerName = _readString(file);
						file.read((char*)&pmtGreen,4);
						file.read((char*)&pmtRed,4);
						scannerVersion = _readString(file);
						imagingUser = _readString(file);
						file.seekg(pos);
						break;
		}
	}
	file.close();
}

string Gtc::_readString(ifstream &file)
{
	char c;
	char buff[256];

	memset(buff,0,256);
	file.get(c);
	int len = (unsigned char)c;
	file.read(buff,len);
	string *s = new string(buff);
	return *s;
}

string Gtc::readString(ifstream &file, int offset)
{
	ios::pos_type pos = file.tellg();
	file.seekg(offset);
	string s = _readString(file);
	file.seekg(pos);
	return s;
}

string Gtc::dump(void)
{
	ostringstream s;
	s << "Number of SNPs:   " << numSnps << endl
			 << "Sample Name:      " << sampleName << endl
	         << "Sample Plate:     " << samplePlate << endl
			 << "Sample Well:      " << sampleWell << endl
			 << "Cluster File:     " << clusterFile << endl
	         << "Manifest:         " << manifest << endl
			 << "Imaging Date:     " << imagingDate << endl
			 << "Autocall Date:    " << autocallDate << endl
			 << "Autocall Version: " << autocallVersion << endl
			 << "Scanner Name:     " << scannerName << endl
			 << "Pmt Green:        " << pmtGreen << endl
			 << "Pmt Red:          " << pmtRed << endl
			 << "Scanner Version:  " << scannerVersion << endl
			 << "Imaging User:     " << imagingUser << endl
	         << "XForm contains:   " << XForm.size() << " entries" << endl
	         << "xRawControl:      " << xRawControl.size() << " entries" << endl
	         << "yRawControl:      " << yRawControl.size() << " entries" << endl
	         << "xRawIntensity:    " << xRawIntensity.size() << " entries" << endl
	         << "yRawIntensity:    " << yRawIntensity.size() << " entries" << endl
	         << "Genotypes:        " << genotypes.size() << " entries" << endl
	         << "Scores:           " << scores.size() << " entries" << endl
	         << "Base Calls:       " << baseCalls.size() << " entries" << endl
			 << endl;
	return s.str();
}

string Gtc::json_dump(void)
{
	ostringstream s;
	s << "{" 
	  << "\"sample_name\":\"" << sampleName << "\","
	  << "\"sample_plate\":\"" << samplePlate << "\","
	  << "\"sample_well\":\"" << sampleWell << "\","
	  << "\"cluster_file\":\"" << clusterFile << "\","
	  << "\"manifest\":\"" << manifest << "\","
	  << "\"imaging_date\":\"" << imagingDate << "\","
	  << "\"autocall_date\":\"" << autocallDate << "\","
	  << "\"autocall_version\":\"" << autocallVersion << "\","
	  << "\"scanner_name\":\"" << scannerName << "\","
	  << "\"scanner_version\":\"" << scannerVersion << "\","
	  << "\"pmt_green\":\"" << pmtGreen << "\","
	  << "\"pmt_red\":\"" << pmtRed << "\","
	  << "\"imaging_user\":\"" << imagingUser << "\""
	  << "}";
	return s.str();
}

void Gtc::readXForm(ifstream &file, int offset)
{
	int arrayLen;
        int total_reserved = 6; // 'reserved' floats after the xform fields
	float reserved;

        int version;
        float xOffset;
        float yOffset;
        float xScale;
        float yScale;
        float shear;
        float theta;

	ios::pos_type pos = file.tellg();
	file.seekg(offset);
	file.read((char*)&arrayLen,4);

	for (int i=0; i<arrayLen; i++) {
            file.read((char*)&version, 4);
            file.read((char*)&xOffset, 4);
            file.read((char*)&yOffset, 4);
            file.read((char*)&xScale, 4);
            file.read((char*)&yScale, 4);
            file.read((char*)&shear, 4);
            file.read((char*)&theta, 4);
            XFormClass X = XFormClass(version, xOffset, yOffset, xScale, 
                                      yScale, shear, theta);
            this->XForm.push_back(X);
            for (int j=0; j<total_reserved; j++) {
              file.read((char*)&reserved, 4);
            }
	}
	file.seekg(pos);
}

void Gtc::readShortArray(ifstream &file, int offset, vector<unsigned short> &a)
{
	int arrayLen;
	unsigned short x;

	ios::pos_type pos = file.tellg();
	file.seekg(offset);
	file.read((char*)&arrayLen, 4);

	for (int n=0; n<arrayLen; n++) {
		file.read((char*)&x,2);
		a.push_back(x);
	}
	file.seekg(pos);
}

void Gtc::readByteArray(ifstream &file, int offset, vector<char> &a)
{
	int arrayLen;
	char c;

	ios::pos_type pos = file.tellg();
	file.seekg(offset);
	file.read((char*)&arrayLen, 4);

	for (int n=0; n<arrayLen; n++) {
		file.get(c);
		a.push_back(c);
	}
	file.seekg(pos);
}


void Gtc::readBaseCallArray(ifstream &file, int offset, vector<BaseCallClass> &a)
{
	int arrayLen;
	BaseCallClass b;
	char c;

	ios::pos_type pos = file.tellg();
	file.seekg(offset);
	file.read((char*)&arrayLen, 4);

	for (int n=0; n<arrayLen; n++) {
		file.get(c); b.a = c;
		file.get(c); b.b = c;
		a.push_back(b);
	}
	file.seekg(pos);
}


void Gtc::readFloatArray(ifstream &file, int offset, vector<float> &a)
{
	int arrayLen;
	float f;

	ios::pos_type pos = file.tellg();
	file.seekg(offset);
	file.read((char*)&arrayLen, 4);

	for (int n=0; n<arrayLen; n++) {
		file.read((char*)&f,4);
		a.push_back(f);
	}
	file.seekg(pos);
}


double Gtc::passRate(double cutOff) 
{
	int pass=0;
	for (vector<float>::iterator s = scores.begin(); s != scores.end(); s++) {
		if (*s >= cutOff) pass++;
	}
	if (scores.size() == 0) return 0;
	return (double)pass / (double)scores.size() * 100.0;
}

double Gtc::correctedPassRate(double cutOff) 
{
	int pass=0;
	int fail=0;
	for (vector<float>::iterator s = scores.begin(); s != scores.end(); s++) {
		if (*s > 0.00000001) {
			if (*s >= cutOff) pass++;
			else              fail++;
		}
	}
	if (pass+fail == 0) return 0.0;
	return (double)pass / (double)(pass+fail) * 100.0;
}


