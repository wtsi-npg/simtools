// 
// sim.cpp
//
// Program to read the contents of a GTC file
//
// $Id: sim.cpp 1354 2010-11-11 16:20:09Z js10 $
//
// Copyright (c) 2009 - 2010 Genome Research Ltd.
//
// Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
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

#include <iostream>
#include <vector>
#include <map>

#include "Sim.h"
#include "Manifest.h"

using namespace std;

int main(int argc, char *argv[])
{
	Sim *sim = new Sim();
string m;
		cout << endl << "Reading SIM file: " << argv[1] << endl;
try { sim->open(argv[1]); } 
catch (string m) { cout << m << endl; }

	cout << "Magic:     " << sim->magic << endl;
	cout << "Version:   " << (int)sim->version << endl;
	cout << "Name Size: " << sim->sampleNameSize << endl;
	cout << "Samples:   " << sim->numSamples << endl;
	cout << "Probes:    " << sim->numProbes << endl;
	cout << "Channels:  " << (int)sim->numChannels << endl;
	cout << "Format:    " << (int)sim->numberFormat << endl;

	cout << endl;
	char *sampleName = new char[sim->sampleNameSize+1];
	vector<uint16_t> *intensity = new vector<uint16_t>;;
	for (unsigned int n=0; n < sim->numSamples; n++) {
		intensity->clear();
		sim->getNextRecord(sampleName, intensity);
		cout << sampleName << "\t: ";
		for (vector<uint16_t>::iterator i = intensity->begin(); i != intensity->end(); i++) {
			cout << *i << " ";
		}
		cout << endl;
	}

	return 0;
}

