// 
// gtc.cpp
//
// Program to read the contents of a GTC file
//
// $Id: gtc.cpp 1510 2013-01-21 11:29:16Z js10 $
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

#include <iostream>
#include <map>

#include "Gtc.h"
#include "Manifest.h"

using namespace std;

int main(int argc, char *argv[])
{
	Gtc *gtc = new Gtc();

		cout << endl << "Reading GTC file: " << argv[1] << endl;
		gtc->open(argv[1], Gtc::ALL);
		if (!gtc->errorMsg.empty()) {
			cout << gtc->errorMsg << endl;
			exit(1);
		}
		cout << gtc->dump() << endl;
	cout << "Pass Rate: " << gtc->passRate(.15) << endl;
	cout << "Corrected Pass Rate: " << gtc->correctedPassRate(.15) << endl;

	cout << "Base calls:" << endl;
	for (int i=85435; i<85445; i++) {
		cout << gtc->baseCalls[i].a << gtc->baseCalls[i].b << " ";
	}
	cout << endl << endl;

#if 0
	cout << "Genotypes:" << endl;
	for (int i=0; i<20; i++) {
		cout << gtc->genotypes[i] << " ";
	}
	cout << endl << endl;
#endif

 	cout << "XForm:" << endl;
	cout << "idx\tversion\txOffset\tyOffset\txScale\tyScale\tshear\ttheta" << endl;
	int n=0;
	vector<XFormClass>::iterator pos;
	for (pos = gtc->XForm.begin(); pos < gtc->XForm.end(); pos++) {
		cout << n << "\t" << pos->version << "\t" << pos->xOffset << "\t" << pos->yOffset << "\t"
		     << pos->xScale << "\t" << pos->yScale << "\t"
		     << pos->shear << "\t" << pos->theta
		     << endl;
		n++;
	}
	cout << endl;

	cout << "Raw Intensities" << endl;
	for (int i=0; i<20; i++) {
		cout << "Idx = " << i << "\tX = " << gtc->xRawIntensity[i] << "\tY = " << gtc->yRawIntensity[i] << endl;
	}
	cout << endl;

	// Read the manifest
	string f = argv[1];
	f = f.substr(0,f.find_last_of('/'));
	f = f + "/../" + gtc->manifest + ".csv";
	Manifest *m = new Manifest();
	try { 
		m->open(f); 
	}
	catch (string s) {
		cout << s << endl;
	}
//	m->dump();

//	int n = m->snp2idx(argv[2]);
//	snpClass *s = m->findSNP(argv[2]);
//	cout << "SNP = " << argv[2] << "  n = " << n << "  s = " << s << endl;

	unsigned int max = 0;
	long p = 0;
	for (vector<snpClass>::iterator snp = m->snps.begin(); snp != m->snps.end(); snp++) {
		if (max < snp->name.length()) max = snp->name.length();
		if (p < snp->position) p = snp->position;
	}

	cout << endl << "Maximum length of SNP name is " << max << endl;
	cout << endl << "Maximum position of SNP is " << p << endl;

	cout << endl;
	cout << "    size = " << gtc->xRawIntensity.size() << endl;
	cout << "max_size = " << gtc->xRawIntensity.max_size() << endl;
	cout << "capacity = " << gtc->xRawIntensity.capacity() << endl;

	gtc->xRawIntensity.reserve(0);
	cout << "    size = " << gtc->xRawIntensity.size() << endl;
	cout << "max_size = " << gtc->xRawIntensity.max_size() << endl;
	cout << "capacity = " << gtc->xRawIntensity.capacity() << endl;

	return 0;
}

