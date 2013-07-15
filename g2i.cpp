//
// g2i: Process one or more GTC files to produce a file suitable for being read
//      by illuminus.
//
// Author: Jennifer Liddle (js10)
//
// $Id: g2i.cpp 1511 2013-01-21 15:37:42Z js10 $
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

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <hash_map>

#include "Gtc.h"
#include "Manifest.h"
#include "win2unix.h"
#include "Sim.h"
#include "json/json.h"
#include "plink_binary.h"

#define CACHESIZE 256
#define GCCACHESIZE 32

using namespace std;

Json::Value root;   // will contains the root value after parsing.
Json::Reader reader;

hash_map<string,string> gtcHash;	// <sample_name, filename>
Manifest *manifest = new Manifest();
vector<string> sampleArray;
hash_map<string,float*> cache;
hash_map<string,string> gcCache;
hash_map<string,int> exclusionList;	// List of samples to exclude
typedef hash_map<string,fstream*> FileMap;
FileMap outFile;
FileMap::iterator pos;
vector<string> filenameArray;
vector<string> sampleNames;
vector<string> gender_code;

typedef hash_map<string,ios::pos_type> PosMap;
PosMap filePos;
PosMap::iterator fp;

Gtc gtc;

// command line options
bool verbose=false;
bool bed=false;
bool genoSNP=false;
bool genoCalling = false;
bool simOutput = false;
bool excludeCnv=false;
bool includeCnv=false;
bool excludeExo=false;
bool normalise=true;
string outputFile;
string project;
string manifestDir;
string exclusionFile;
string chrSelect;
string tmpFile;

char * timestamp(void)
{
	static char buffer[64];
	time_t t = time(NULL);
	strftime(buffer,64,"%F %T: ", localtime(&t));
	return buffer;
}

void displayUsage(void)
{
	cout << "g2i: convert Illumina GTC files into a format suitable for loading into illuminus or genoSNP, or export genotype calls, or create BED files"
	     << endl << endl
	     << "Usage:" << endl
	     << "g2i -o <filename> [options] <gtcfile> <gtcfile> ..." << endl
	     << "or" << endl
	     << "g2i -o <filename> -i <filename> [options]" << endl << endl
	     << "Where 'options' are one or more of:" << endl
	     << "-o <prefix>      specifies the output filename to create" << endl
	     << "-t <prefix>      Specify a temporary file location to use" << endl
	     << "-v               verbose mode" << endl
	     << "-h               displays this message" << endl
	     << "-g               produce  genoSNP input rather than Illuminus" << endl
	     << "-s               produce genotype calling data rather than Illuminus" << endl
	     << "-m               produce SIM file rather than Illuminus" << endl
	     << "-b               produce BED file rather than Illuminus" << endl
	     << "-r <chromosome>  select samples for this chromosome only" << endl
	     << "-p <project>     extract data for samples in this project ID from the Illumina LIMS" << endl
	     << "-n               do NOT perform normalisation on intensities" << endl
	     << "-c               exclude any SNPs beginning with cnv" << endl
	     << "-k               include any SNPs beginning with cnv" << endl
	     << "-e               exclude any samples beginning with 'Exo' or 'EXO'" << endl
	     << "-x <filename>    exclude any samples contained in <filename>" << endl
	     << "-i <filename>    specifies a file containing a list of GTC files to process" << endl
         << "-d <directory>   specifies an alternative manifest file location" << endl
	     << "<gtcfile>...     is a list of GTC files to process" << endl << endl;
	exit(0);
}


void loadExclusionFile(string fname) 
{
	ifstream f; 
	string s;
	f.open(fname.c_str());
	while (f >> s) {
		exclusionList[s] = 1;
		if (verbose) cout << timestamp() << "Excluding: " << s << endl;
	}
	f.close();
}

// Sort function to sort SNPs by position
bool SortByPosition(const snpClass &snp1, const snpClass &snp2)
{
	if (snp1.chromosome.compare(snp2.chromosome))
		return (snp1.chromosome.compare(snp2.chromosome) < 0);
	if (snp1.position != snp2.position)
		return (snp1.position < snp2.position);
	return (snp1.name.compare(snp2.name) < 0);
}

// Sort function to sort SNPs by position for BED file
bool SortByPositionBED(const snpClass &snp1, const snpClass &snp2)
{

	int c1 = atoi(snp1.chromosome.c_str());
	int c2 = atoi(snp2.chromosome.c_str());
	if (c1 && c2 && (c1!=c2))
		return (c1 < c2);
	if (snp1.chromosome.compare(snp2.chromosome))
		return (snp1.chromosome.compare(snp2.chromosome) < 0);
	if (snp1.position != snp2.position)
		return (snp1.position < snp2.position);
	return (snp1.name.compare(snp2.name) < 0);
}

void flushCache(int cacheIndex)
{
	if (verbose) cout << timestamp() << "Flushing cache..." << endl;
	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
		if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
		// look up the file and position for this SNP
		fstream *f = outFile[snp->chromosome];
		f->seekp(filePos[snp->name]);

		ostringstream os;
		for (int i=0; i<cacheIndex; i++) {
			if (cache[snp->name][i] < 0) cache[snp->name][i] = 0; 
			os << '\t' << setw(7) << std::fixed << setprecision(3) << cache[snp->name][i];
		}
		f->write(os.str().c_str(), os.str().size());
		filePos[snp->name] += os.str().size();
	}
}

void flushgcCache(int cacheIndex)
{
	if (verbose) cout << timestamp() << "Flushing cache..." << endl;
	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
		if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
		// look up the file and position for this SNP
		fstream *f = outFile[snp->chromosome];
		f->seekp(filePos[snp->name]);

		ostringstream os;
		os << gcCache[snp->name];
		f->write(os.str().c_str(), os.str().size());
		filePos[snp->name] += os.str().size();
		gcCache[snp->name] = "";
	}
}

//
// We've read the Manifest and all the GTC files
// Now it's time to create the output files
//
void goForIt(string fname)
{
	int recordLength = gtcHash.size() * 10 * 2;
	char *buffer = new char[recordLength];
	memset(buffer,' ',recordLength);
	buffer[recordLength-1] = '\n';

	// Sort the SNPs into position order
	sort(manifest->snps.begin(), manifest->snps.end(), SortByPosition);

	// Create lockfile
	string lockFileName = fname + ".lock";
	FILE *lockfile = fopen(lockFileName.c_str(), "w");
	if (!lockfile) {
		cerr << "Can't create lock file " << lockFileName << endl;
		cerr << strerror(errno) << endl;
		exit(1);
	}
	fclose(lockfile);

	//
	// Create all of the output files - one for each chromosome
	//
	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
		if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
		fstream *f = outFile[snp->chromosome];
		if (!f) {
			f = new fstream();
			string fullFname = fname + "_intu_" + snp->chromosome + ".txt";
			filenameArray.push_back("_intu_" + snp->chromosome + ".txt");
			if (verbose) cout << timestamp() << "creating file " << fullFname << endl;
			f->open(fullFname.c_str(), ios::in | ios::out | ios::trunc);

			*f << "SNP\tCoor\tAlleles";
			// write sample names from all the gtc files
			for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
				*f << "\t" << i->first << "A\t" << i->first << "B";
			}
			*f << endl;
			// associate the file handle with the chromosome
			outFile[snp->chromosome] = f;
		}
		f = outFile[snp->chromosome];
		*f << snp->name << "\t" << snp->position << "\t" << snp->snp[0] << snp->snp[1];
		filePos[snp->name] = f->tellp();	// store next position to write
		f->write(buffer,recordLength);	// fill with nulls (or spaces)

		cache[snp->name] = new float[CACHESIZE];
	}

	//
	// Process each GTC file in turn
	//
	int n=1;
	int cacheIndex = 0;
	for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
		if (verbose) cout << timestamp() << "Processing GTC file " << n++ << " of " << gtcHash.size() << endl;
		gtc.open(i->second,Gtc::XFORM | Gtc::INTENSITY);	// reload GTC file to read XForm and Intensity arrays

		for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
			if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
			if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
			int idx = snp->index - 1;	// index is zero based in arrays, but starts from 1 in the map file
			unsigned int norm = manifest->normIdMap[snp->normId];
			XFormClass *XF = &gtc.XForm[norm];

			double xn, yn;
			if (normalise) {
				// first do the normalisation calculation
				double tempx = gtc.xRawIntensity[idx] - XF->xOffset;
				double tempy = gtc.yRawIntensity[idx] - XF->yOffset;

				double cos_theta = cos(XF->theta);
				double sin_theta = sin(XF->theta);
				double tempx2 = cos_theta * tempx + sin_theta * tempy;
				double tempy2 = -sin_theta * tempx + cos_theta * tempy;

				double tempx3 = tempx2 - XF->shear * tempy2;
				double tempy3 = tempy2;

				xn = tempx3 / XF->xScale;
				yn = tempy3 / XF->yScale;
			} else {
				xn = gtc.xRawIntensity[idx];
				yn = gtc.yRawIntensity[idx];
			}

			cache[snp->name][cacheIndex] = xn;
			cache[snp->name][cacheIndex+1] = yn;
		}
		cacheIndex += 2;
		if (cacheIndex == CACHESIZE) { flushCache(cacheIndex); cacheIndex=0; }
	}

	flushCache(cacheIndex);

	// close all of the files
	for (pos = outFile.begin(); pos != outFile.end(); pos++) {
		pos->second->close();
	}

	// delete lockfile and create donefile
	string doneFileName = lockFileName;
	string::size_type dot = doneFileName.find(".lock");
	doneFileName.replace(dot, 5, ".g2i");
	rename(lockFileName.c_str(), doneFileName.c_str());
	if (verbose) cout << timestamp() << "Renamed " << lockFileName << " to " << doneFileName << endl;
}

// Read the manifest
//
// manifestPath is the full pathname of the manifest directory, either
// calculated relative to the first GTC file (the default) or suppled
// on the command line. manifestName is the name of the manifest file
// as stored in the GTC file
//
void loadManifest(string manifestPath, string manifestName)
{
    string mfile = manifestPath + '/' + manifestName + ".csv";
	if (verbose) cout << timestamp() << "Reading manifest: " << mfile << endl;
	try { 
        manifest->open(mfile); 
	}
	catch (string s) {
		cout << s << endl;
		cerr << s << endl;
	}
}

//
// Return a manifest file path relative to GTC file name gtcName
//
//
string relativeManifestPath(string gtcName)
{
  string rpath = gtcName;
  size_t i =  rpath.find('/');

  if (i == string::npos) {
    rpath = '.';
  }
  else {
    rpath = rpath.substr(0,rpath.find_last_of('/'));
  }

  rpath = rpath + "/..";

  return rpath;
}


void getAutocallInfo(vector<string> &infiles, string project)
{
	char buffer[1024];

        // project may contain space chars ...
        if (project.find("\\") == string::npos) {
            // ... that haven't been escaped
            size_t pos = project.find(" ");
            while (pos != string::npos) {
                // prepend with "\"
                project.replace(pos,1,"\\ ");
                // next search starts after current position
                // (pos + 1), but we've also replace a 1-char
                // string with a 2-char one, so check pos + 2
                pos = project.find(" ", pos + 2);
            }
        }
	string pname = "autocall_report.pl -t -p "+project;
	if (verbose) cout << timestamp() << pname << endl;
	FILE *f = popen(pname.c_str(), "r");
	if (!f) {
		cout << "Oops! " << strerror(errno) << endl;;
		cerr << "Oops! " << strerror(errno) << endl;;
		exit(1);
	}

	while (fgets(buffer,1024,f)) {
                // Anorexia_670    WG0015336-DNA   A02     144693_A02_CCC3_AN1208139       5065997151      R01C01  AutoCall Completed      Pass    08/05/10        /nfs/illumina_geno04/call/20100508/5065997151_R01C01.gtc
		char *p = strtok(buffer,"\t");
		int n = 1;
		while (p) {
			if (n==10 && strlen(p) > 2) {
				string filepath = p;
				filepath = win2unix(filepath);
				infiles.push_back(filepath);
				if (verbose) cout << timestamp() << "Read " << p << " ==> " << filepath << endl;
			}
			n++;
			p = strtok(NULL,"\t");
		}
	}
	
}

void createBedFile(string fname, vector<string>infiles)
{
	vector<int> genotypes;
	plink_binary *pb = new plink_binary();
	pb->bed_mode = 0;	// specify individual major
	pb->open(fname,1);

	// Sort the SNPs into position order
	sort(manifest->snps.begin(), manifest->snps.end(), SortByPositionBED);

	// Load the SNP names into gftools
	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		gftools::snp gfsnp;
		if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
		if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
		gfsnp.name = snp->name;
		gfsnp.chromosome = snp->chromosome;
		if (snp->chromosome == "X") gfsnp.chromosome = "23";
		if (snp->chromosome == "Y") gfsnp.chromosome = "24";
		if (snp->chromosome == "XY") gfsnp.chromosome = "25";
		if (snp->chromosome == "MT") gfsnp.chromosome = "26";
		gfsnp.physical_position = snp->position;
		if (snp->snp[0] == 'N') {
			gfsnp.allele_a = '?';
			gfsnp.allele_b = '?';
		} else {
			gfsnp.allele_a = snp->snp[0];
			gfsnp.allele_b = snp->snp[1];
		}
		pb->snps.push_back(gfsnp);
	}

	if (verbose) cerr << "Pushed " << pb->snps.size() << " SNPs" << endl;
	//
	// Process each GTC file in turn
	//
	for (unsigned int n = 0; n < infiles.size(); n++) {
//	for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
		if (verbose) cout << timestamp() << "Processing GTC file " << n+1 << " of " << infiles.size() << endl << infiles[n] << endl;
		gtc.open(infiles[n],Gtc::GENOTYPES | Gtc::BASECALLS | Gtc::SCORES);	// reload GTC file to read required arrays

		gftools::individual ind;
		if (!sampleNames.empty()) ind.name = sampleNames[n];
		else                      ind.name = gtc.sampleName;
		if (!gender_code.empty()) ind.sex = gender_code[n]; 
		pb->individuals.push_back(ind);

		for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
			if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
			if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
			int idx = snp->index - 1;	// index is zero based in arrays, but starts from 1 in the map file
			char buff[3];
		    sprintf(buff,"%c%c", gtc.baseCalls[idx].a, gtc.baseCalls[idx].b);
			char allele_a = snp->snp[0];
			char allele_b = snp->snp[1];
			int call = -1;
			if (buff[0] == '-' && buff[1] == '-') call=0;
			if (buff[0] == allele_a && buff[1] == allele_a) call = 1;
			if (buff[0] == allele_a && buff[1] == allele_b) call = 2;
			if (buff[0] == allele_b && buff[1] == allele_a) call = 2;
			if (buff[0] == allele_b && buff[1] == allele_b) call = 3;
			if (call == -1) {
				cerr << "malformed data: " << endl;
				cerr << "snp = '" << allele_a << allele_b << "'" << endl;
				cerr << "buff = '" << buff << "'" << endl;
				cerr << "name = " << snp->name << endl;
				cerr << "idx = " << idx << endl;
				exit(1);
			}
			genotypes.push_back(call);
		}
		if (verbose) cerr << "SNP Count = " << pb->snps.size() << endl;
		pb->write_individual(genotypes);
		genotypes.clear();
	}

	pb->close();
}

void createSimFile(string fname)
{
	Sim *sim = new Sim();

	hash_map<string,string>::iterator i = gtcHash.begin();
	gtc.open(i->second, Gtc::INTENSITY);
	sim->createFile(fname);
	sim->writeHeader(gtcHash.size(), gtc.xRawIntensity.size());

	//
	//
	// Process each GTC file in turn
	//
	unsigned int n=1;
	for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
		char *buffer;
		if (verbose) cout << timestamp() << "Processing GTC file " << n << " of " << gtcHash.size() << endl;
		//
		// add sample name to each output file
		// no family info as yet (todo?) - write sample ID twice
//		fn << i->first << endl;

		gtc.open(i->second,Gtc::XFORM | Gtc::INTENSITY);	// reload GTC file to read XForm and Intensity arrays
		buffer = new char[sim->sampleNameSize];
		memset(buffer,0,sim->sampleNameSize);
		// if we have a sample name from the json file, use it
		if (sampleNames.size() > (n-1)) { strcpy(buffer, sampleNames[n-1].c_str()); }
		else                            { strcpy(buffer,gtc.sampleName.c_str()); }
		sim->write(buffer, sim->sampleNameSize);

		for (unsigned int idx = 0; idx < gtc.xRawIntensity.size(); idx++) {
			uint16_t v;
			v = gtc.xRawIntensity[idx];
			sim->write(&v,sizeof(v));
			v = gtc.yRawIntensity[idx];
			sim->write(&v,sizeof(v));
		}
		n++;

#if 0
		for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
			if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
			if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
			int idx = snp->index - 1;	// index is zero based in arrays, but starts from 1 in the map file
			unsigned int norm = manifest->normIdMap[snp->normId];
			XFormClass *XF = &gtc.XForm[norm];

			// first do the normalisation calculation
			double tempx = gtc.xRawIntensity[idx] - XF->xOffset;
			double tempy = gtc.yRawIntensity[idx] - XF->yOffset;

			double cos_theta = cos(XF->theta);
			double sin_theta = sin(XF->theta);
			double tempx2 = cos_theta * tempx + sin_theta * tempy;
			double tempy2 = -sin_theta * tempx + cos_theta * tempy;

			double tempx3 = tempx2 - XF->shear * tempy2;
			double tempy3 = tempy2;

			double xn = tempx3 / XF->xScale;
			double yn = tempy3 / XF->yScale;

			// add raw/norm x/y to .raw and .nor files
//			fn << "\t" << std::fixed << setprecision(3) << xn << " " << yn;
		}
#endif

	}
	sim->close();
}


void createGenoSNP(string fname)
{
	// Sort the SNPs into position order
	sort(manifest->snps.begin(), manifest->snps.end(), SortByPosition);

	// Create lockfile
	string lockFileName = fname + ".lock";
	FILE *lockfile = fopen(lockFileName.c_str(), "w");
	fclose(lockfile);

	//
	// Create all of the output files
	//
	string rawFname = fname + ".raw.txt";
	string snpFname = fname + ".snp.txt";
	string norFname = fname + ".nor.txt";
	filenameArray.push_back(".raw.txt");
	filenameArray.push_back(".snp.txt");
	filenameArray.push_back(".nor.txt");
	
	fstream fn;
	fstream fs;
	fstream fr;
	fr.open(rawFname.c_str(), ios::in | ios::out | ios::trunc);
	fn.open(norFname.c_str(), ios::in | ios::out | ios::trunc);
	fs.open(snpFname.c_str(), ios::in | ios::out | ios::trunc);
	if (verbose) cout << timestamp() << "creating file " << rawFname << ", " << norFname << " and " << snpFname << endl;

        // Write SNP list to .snp
	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
		if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
		fs << snp->name << '\t'
		   << (snp->normId % 100) + 1 << '\t'
		   << snp->snp[0] << " " << snp->snp[1] << endl;
	}
        fs.close();
        if (fs.bad())
            cout << "Error closing snp file" << endl;

	//
	// Process each GTC file in turn
	//
	int n=1;
	for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
		if (verbose) cout << timestamp() << "Processing GTC file " << n++ << " of " << gtcHash.size() << endl;
		//
		// add sample name to each output file
		// no family info as yet (todo?) - write sample ID twice
		fn << i->first << "\t" << i->first;
		fr << i->first << "\t" << i->first;

		gtc.open(i->second,Gtc::XFORM | Gtc::INTENSITY);	// reload GTC file to read XForm and Intensity arrays

		for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
			if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
			if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
			int idx = snp->index - 1;	// index is zero based in arrays, but starts from 1 in the map file
			unsigned int norm = manifest->normIdMap[snp->normId];
			XFormClass *XF = &gtc.XForm[norm];

			// first do the normalisation calculation
			double tempx = gtc.xRawIntensity[idx] - XF->xOffset;
			double tempy = gtc.yRawIntensity[idx] - XF->yOffset;

			double cos_theta = cos(XF->theta);
			double sin_theta = sin(XF->theta);
			double tempx2 = cos_theta * tempx + sin_theta * tempy;
			double tempy2 = -sin_theta * tempx + cos_theta * tempy;

			double tempx3 = tempx2 - XF->shear * tempy2;
			double tempy3 = tempy2;

			double xn = tempx3 / XF->xScale;
			double yn = tempy3 / XF->yScale;

			// add raw/norm x/y to .raw and .nor files
			fr << "\t" << std::fixed << setprecision(3) << gtc.xRawIntensity[idx] << " " << gtc.yRawIntensity[idx];
			fn << "\t" << std::fixed << setprecision(3) << xn << " " << yn;
		}
		fn << endl;
		fr << endl;
	}

	fn.close();
	if (fn.bad()) cout << "Error closing nor file" << endl;
	fr.close();
	if (fr.bad()) cout << "Error closing raw file" << endl;

	// delete lockfile and create donefile
	string doneFileName = lockFileName;
	string::size_type dot = doneFileName.find(".lock");
	doneFileName.replace(dot, 5, ".g2i");
	rename(lockFileName.c_str(), doneFileName.c_str());
	if (verbose) cout << timestamp() << "Renamed " << lockFileName << " to " << doneFileName << endl;
}

void createGenoCalling(string fname)
{
	int recordLength = gtcHash.size() * 18;
	char *buffer = new char[recordLength];
	memset(buffer,' ',recordLength);
	buffer[recordLength-1] = '\n';

	// Sort the SNPs into position order
	sort(manifest->snps.begin(), manifest->snps.end(), SortByPosition);

	// Create lockfile
	string lockFileName = fname + ".lock";
	FILE *lockfile = fopen(lockFileName.c_str(), "w");
	if (!lockfile) throw (strerror(errno));
	fclose(lockfile);

	//
	// Create all of the output files - one for each chromosome
	//
	for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
		if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
		if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
		fstream *f = outFile[snp->chromosome];
		if (!f) {
			f = new fstream();
			string fullFname = fname + "_gtu_" + snp->chromosome + ".txt";
			filenameArray.push_back("_gtu_" + snp->chromosome + ".txt");
			
			if (verbose) cout << timestamp() << "creating file " << fullFname << endl;
			f->open(fullFname.c_str(), ios::in | ios::out | ios::trunc);

			// write sample names from all the gtc files
			for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
				*f << "\t" << i->first;
			}
			*f << endl;
			// associate the file handle with the chromosome
			outFile[snp->chromosome] = f;
		}
		f = outFile[snp->chromosome];
		*f << snp->name;
		filePos[snp->name] = f->tellp();	// store next position to write
		f->write(buffer,recordLength);	// fill with nulls (or spaces)

		gcCache[snp->name] = "";
	}

	//
	// Process each GTC file in turn
	//
	int n=1;
	int cacheIndex = 0;
	for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
		if (verbose) cout << timestamp() << "Processing GTC file " << n++ << " of " << gtcHash.size() << endl;
		gtc.open(i->second,Gtc::GENOTYPES | Gtc::BASECALLS | Gtc::SCORES);	// reload GTC file to read required arrays

		for (vector<snpClass>::iterator snp = manifest->snps.begin(); snp != manifest->snps.end(); snp++) {
			if (excludeCnv && snp->name.find("cnv") != string::npos) continue;
			if (chrSelect.size() && chrSelect.compare(snp->chromosome)) continue;
			int idx = snp->index - 1;	// index is zero based in arrays, but starts from 1 in the map file
			char buffer[128];
                        if (gtc.genotypes[idx] > 3)
			   cout << "Unknown genotype value: " << gtc.genotypes[idx] << endl;
                        if (gtc.genotypes[idx] == 0)
			    sprintf(buffer,"\tNN;%f", gtc.scores[idx]);
                        else
			    sprintf(buffer,"\t%c%c;%f", gtc.baseCalls[idx].a, gtc.baseCalls[idx].b, gtc.scores[idx]);
			gcCache[snp->name] += buffer;
		}
		cacheIndex++;
		if (cacheIndex == GCCACHESIZE) { flushgcCache(cacheIndex); cacheIndex=0; }
	}

	flushgcCache(cacheIndex);

	// close all of the files
	for (pos = outFile.begin(); pos != outFile.end(); pos++) {
		pos->second->close();
	}

	// delete lockfile and create donefile
	string doneFileName = lockFileName;
	string::size_type dot = doneFileName.find(".lock");
	doneFileName.replace(dot, 5, ".g2i");
	rename(lockFileName.c_str(), doneFileName.c_str());
	if (verbose) cout << timestamp() << "Renamed " << lockFileName << " to " << doneFileName << endl;
}


int main(int argc, char *argv[])
{
	char c;
	int nBadFiles = 0;
	vector<string> infiles;

	while ((c = getopt (argc, argv, "neckmsvbw?hgo:i:x:p:d:r:t:")) != -1) {
		switch (c) {
				case 'v':	verbose=true; break;
				case 'b':	bed=true; break;
				case 'c':	excludeCnv=true; break;
				case 'k':	includeCnv=true; break;
				case 'e':	excludeExo=true; break;
				case 'g':	genoSNP=true; break;
				case 's':	genoCalling=true; break;
				case 'm':	simOutput=true; break;
				case 'n':	normalise=false; break;
				case 'h':
				case '?':	displayUsage(); break;
				case 'o':	outputFile = optarg;	break;
				case 't':	tmpFile = optarg;	break;
				case 'p':	project = optarg; break;
                case 'd':   manifestDir = optarg; break;
                case 'x':   exclusionFile = optarg; break;
				case 'r':	chrSelect = optarg; break;
				case 'i':	ifstream f; string s;
							f.open(optarg);
							if (strstr(optarg,".json")) {
								bool parsingSuccessful = reader.parse( f, root );
								if ( !parsingSuccessful ) {
									std::cout  << "Failed to parse json file\n" << reader.getFormatedErrorMessages() << endl;
									return 0;
								}
								for ( unsigned int index = 0; index < root.size(); ++index ) { // Iterates over the sequence elements.
									sampleNames.push_back(root[index]["uri"].asString());
//									cerr << "Gender Code = [" << root[index]["gender_code"] << "]" << endl;
									char buffer[8];
//									if (root[index]["gender_code"]) sprintf(buffer,"%d",root[index]["gender_code"].asInt());
//									else                            strcpy(buffer,"-9");
									sprintf(buffer,"%d",root[index]["gender_code"].asInt());
									gender_code.push_back(buffer);
									infiles.push_back(root[index]["result"].asString());
cerr << root[index]["uri"].asString() << "\t" << root[index]["result"].asString() << "\t" << root[index]["gender_code"].asInt() << endl;
								}
							} else {
								while (f >> s) infiles.push_back(s);
							}
							f.close();
							break;
		}
	}

	if (outputFile.size() == 0) {
		cout << "No output file specified" << endl;
		displayUsage();
	}

	if (tmpFile.size() == 0) {
		tmpFile = outputFile;
	}

	if (project.size()) {
		getAutocallInfo(infiles,project);
	}

	if (exclusionFile.size() != 0) {
		loadExclusionFile(exclusionFile);
	}

	// read the rest of the command line
	for (int n = optind; n < argc; n++) {
		infiles.push_back(argv[n]);
	}

	if (infiles.size() == 0) {
		cout << "No input files specified" << endl;
		displayUsage();
	}

	if (excludeCnv && includeCnv) {
		cout << "You can't specify -c AND -k options together! That doesn't make sense!" << endl;
		displayUsage();
	}

	if (!excludeCnv && !includeCnv) excludeCnv=true;	// default to 'exclude' if nothing specified

	// sanity check the GTC files
	if (verbose) cout << timestamp() << "Sanity checking";

	for (vector<string>::iterator f = infiles.begin(); f != infiles.end(); f++) {
		bool badFile = false;
		if (verbose) cout << '.';
		try {
			gtc.open(*f,0);
			if (gtc.errorMsg.length()) throw gtc.errorMsg;
		}
		catch (string s) {
			cout << s << endl;
			badFile = true;
			nBadFiles++;
		}
		if (excludeExo && ( (gtc.sampleName.find("Exo") == 0) || 
		                    (gtc.sampleName.find("EXO") == 0))) {
			badFile=true;
			if (verbose) cout << "X";
		}
		if (exclusionList[gtc.sampleName]) {
			badFile = true;
			if (verbose) cout << "X";
		}
		if (!badFile) {
			gtcHash[gtc.sampleName] = *f;
			sampleArray.push_back(gtc.sampleName);
		}
	}
	if (verbose) cout << endl;

	// a little validation...
	if (gtcHash.size() == 0) {
		cout << "No gtc files specified\n";
		exit(1);
	}

	string manifestName = "";
	for (hash_map<string,string>::iterator i = gtcHash.begin(); i != gtcHash.end(); i++) {
		gtc.open(i->second,0);
		if (manifestName != "" && gtc.manifest != manifestName) {
			cout << "GTC files do not all have the same manifest" << endl;
			exit(1);
		}
		manifestName = gtc.manifest;
	}

	if (verbose) {
		cout << timestamp() << "About to process " << gtcHash.size() << " GTC files" << endl;
		if (nBadFiles) cout << "Cannot process " << nBadFiles << " GTC files" << endl;
	}

	if (!simOutput) {
	if (verbose) cout << "Loading manifest " << manifestName << endl;
    
    if (manifestDir.size()) {
      loadManifest(manifestDir, manifestName);
    } else {
      loadManifest(relativeManifestPath(infiles[0]), manifestName);
    }
	}

	try {
	if (genoSNP)          { createGenoSNP(tmpFile); }
	else if (genoCalling) { createGenoCalling(tmpFile);   }
	else if (simOutput)   { createSimFile(tmpFile);   }
	else if (bed)         { createBedFile(tmpFile,infiles);   }
	else                  { goForIt(tmpFile);       }
	}
	catch (char *s) {
		cerr << timestamp() << "Caught fatal error: " << s << endl;
		return 1;
	}

	// Copy temporary files to final output files
	if (tmpFile.compare(outputFile)) {
		for (vector<string>::iterator i = filenameArray.begin(); i != filenameArray.end(); i++) {
			string command = "cp " + tmpFile + *i + " " + outputFile + *i;
			if (verbose) cout << timestamp() << "Sending command: '" << command << "'" << endl;
			int ret = system(command.c_str());
			if (ret) {
				cerr << "system(" << command << ") returned value " << ret;
				exit(1);
			}
                        if (remove((tmpFile + *i).c_str()) != 0)
				cerr << "error cleaning up temporary file" << endl;
		}
	}

	return 0;
}

