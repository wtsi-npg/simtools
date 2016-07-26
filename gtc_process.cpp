//
// g2i: Process one or more GTC files to produce a file suitable for being read
//      by illuminus.
//
// Author: Jennifer Liddle (js10)
//
// $Id: g2i.cpp 813 2009-11-20 09:25:49Z js10 $
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

#include <unistd.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

#include "Gtc.h"
#include "Manifest.h"
#include "win2unix.h"

using namespace std;

bool verbose = false;

//
//
double goForIt(Gtc *gtc, Manifest *manifest)
{
	double meanTotal = 0;
	int n = 0;
        double epsilon = 1e-6;

	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		int idx = snp->index - 1;	// index is zero based in arrays, but starts from 1 in the map file
		unsigned int norm = manifest->normIdMap[snp->normId];
		XFormClass *XF = &gtc->XForm[norm];

		// first do the normalisation calculation
		double tempx = gtc->xRawIntensity[idx] - XF->xOffset;
		double tempy = gtc->yRawIntensity[idx] - XF->yOffset;

		double cos_theta = cos(XF->theta);
		double sin_theta = sin(XF->theta);
		double tempx2 = cos_theta * tempx + sin_theta * tempy;
		double tempy2 = -sin_theta * tempx + cos_theta * tempy;

		double tempx3 = tempx2 - XF->shear * tempy2;
		double tempy3 = tempy2;

		if (abs(XF->xScale) > epsilon && abs(XF->yScale) > epsilon) {
                        // intensities are non-zero (for tolerance epsilon)
			double xn = tempx3 / XF->xScale;
			double yn = tempy3 / XF->yScale;
			if (!std::isnan(xn) && !std::isnan(yn)) {
				meanTotal += (yn-xn);
//cout << xn << '\t' << yn << '\t' << tempx3 << '\t' << tempy3 << '\t' << XF->xScale << '\t' << XF->yScale << '\t' << meanTotal << endl;
				n++;
			}
		}
	}

	return n ? meanTotal / n : 0;
}

// Read the manifest
//
// gtcName is the full pathname of the GTC file, which we use to find the correct directory
// manifestName is the name of the manifest file as stored in the GTC file
//
Manifest *loadManifest(string gtcName, string manifestName)
{
	Manifest *manifest = new Manifest();
	string f = gtcName;
	f = f.substr(0,f.find_last_of('/'));
	f = f + "/../" + manifestName + ".csv";
	if (verbose) cout << "Reading manifest: " << f << endl;
	try { 
		manifest->open(f); 
	}
	catch (string s) {
		cout << s << endl;
	}
	return manifest;
}

double getMeanIntensity(Gtc *gtc, Manifest *manifest)
{
	double mean = goForIt(gtc,manifest);
	return mean;
}

double getIlluminaPassrate(double cutOff, Gtc *gtc, Manifest *manifest)
{
	int pass=0;
	int fail=0;
	int cnv=0;

	if (gtc->scores.size() != manifest->snps.size()) {
		cerr << "Mismatch in sizes: scores = " << gtc->scores.size() << "  snps = " << manifest->snps.size() << endl;
		exit(1);
	}

	if (gtc->scores.size() == 0) return 0;

	vector<float>::iterator s;
	vector<snpClass>::iterator snp;

	for (s = gtc->scores.begin(), snp = manifest->snps.begin(); s != gtc->scores.end(); s++, snp++) {
		if (snp->name.find("cnv") == string::npos) {
			if (*s >= cutOff) pass++;
			else              fail++;
		} else {
			cnv++;
		}
	}

	return (double)pass / (double)(pass+fail) * 100.0;
}

#ifdef TEST
int main(int argc, char *argv[])
{
	string manifestFile = "";
	Manifest *manifest = NULL;
	Gtc *gtc = new Gtc();
	verbose = true;

	for (int n=1; n<argc; n++) {
		gtc->open(win2unix(argv[n]),Gtc::XFORM | Gtc::INTENSITY | Gtc::SCORES);
		if (manifestFile != gtc->manifest) manifest = loadManifest(win2unix(argv[n]),gtc->manifest);
		manifestFile = gtc->manifest;
//		double mean = getMeanIntensity(gtc, manifest);
//		printf("Mean Intensity Difference for %s \t= %lf\n", argv[n], mean);
		double passrate = getIlluminaPassrate(0.15, gtc, manifest);
		printf("Illumina Passrate = %lf\n", passrate);
//		passrate = gtc->passRate(0.15);
//		printf("Passrate = %lf\n", passrate);
//		passrate = gtc->correctedPassRate(0.15);
//		printf("Corrected Passrate = %lf\n", passrate);
	}
}
#endif


