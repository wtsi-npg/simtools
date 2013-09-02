// 
// sim_tools.cpp
//
// Program to read and process contents of a SIM file
//
// $Id: sim.cpp 1354 2010-11-11 16:20:09Z js10 $
//
// Copyright (c) 2009 - 2013 Genome Research Ltd.
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

#include <getopt.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>

#include "Sim.h"
#include "Gtc.h"
#include "Manifest.h"
#include "json/json.h"

using namespace std;

static struct option long_options[] = {
                   {"infile", 1, 0, 0},
                   {"outfile", 1, 0, 0},
                   {"man_file", 1, 0, 0},
                   {"man_dir", 1, 0, 0},
                   {"normalize", 0, 0, 0},
                   {"verbose", 0, 0, 0},
                   {"start", 1, 0, 0},
                   {"end", 1, 0, 0},
                   {0, 0, 0, 0}
               };


// Sort function to sort SNPs by position
bool SortByPosition(const snpClass &snp1, const snpClass &snp2)
{
	if (snp1.chromosome.compare(snp2.chromosome))
		return (snp1.chromosome.compare(snp2.chromosome) < 0);
	if (snp1.position != snp2.position)
		return (snp1.position < snp2.position);
	return (snp1.name.compare(snp2.name) < 0);
}


void showUsage(int argc, char *argv[])
{
	string command = "";

	if (argc > 2) command = argv[2];

	cout << endl;
	if (command == "create") {
		cout << "Usage:   " << argv[0] << " create [options]" << endl << endl;
		cout << "Create a SIM file from a list of GTC files" << endl<< endl;
		cout << "Options: --infile <filename>    File containing list of GTC files to process" << endl;
		cout << "         --outfile <filename>   Name of SIM file to create or '-' for STDOUT" << endl;
		cout << "         --man_file <dirname>    Directory to look for Manifest file in" << endl;
		cout << "         --normalize            Normalize the intensities (default is raw values)" << endl;
		cout << "         --verbose              Show progress messages to STDERR" << endl;
		exit(0);
	}

	if (command == "illuminus") {
		cout << "Usage:   " << argv[0] << " illuminus [options]" << endl << endl;
		cout << "Create an Illuminus file from a SIM file" << endl<< endl;
		cout << "Options: --infile <filename>    Name of SIM file to provess or '-' for STDIN" << endl;
		cout << "         --outfile <filename>   Name of Illuminus file to create or '-' for STDOUT" << endl;
		cout << "         --man_file <dirname>    Directory to look for Manifest file in" << endl;
		cout << "         --start <index>        Which SNP to start processing at (default is to start at the beginning)" << endl;
		cout << "         --end <index>          Which SNP to end processing at (default is to continue until the end)" << endl;
		cout << "         --verbose              Show progress messages to STDERR" << endl;
		exit(0);
	}

	if (command == "genosnp") {
		cout << "Usage:   " << argv[0] << " genosnp [options]" << endl << endl;
		cout << "Create a GenoSNP file from a SIM file" << endl<< endl;
		cout << "Options: --infile   The name of the SIM file or '-' for STDIN" << endl;
		cout << "         --outfile  Name of GenoSNP file to create or '-' for STDOUT" << endl;
		cout << "         --start <index>        Which sample to start processing at (default is to start at the beginning)" << endl;
		cout << "         --end <index>          Which sample to end processing at (default is to continue until the end)" << endl;
		cout << "         --verbose              Show progress messages to STDERR" << endl;
		exit(0);
	}

	if (command == "view") {
		cout << "Usage:   " << argv[0] << " view [options]" << endl << endl;
		cout << "Display some information from a SIM file" << endl<< endl;
		cout << "Options: --infile   The name of the SIM file or '-' for STDIN" << endl;
		cout << "         --verbose  Display intensities as well as header information" << endl;
		exit(0);
	}

	if (command == "qc") {
		cout << "Usage:   " << argv[0] << " qc [options]" << endl << endl;
		cout << "Compute genotyping QC metrics and write to text files" << endl<< endl;
		cout << "Options: --infile      The name of the SIM file or '-' for STDIN" << endl;
		cout << "         --magnitude   Output file for sample magnitude (normalised by SNP)" << endl;
		cout << "         --xydiff      Output file for XY intensity difference" << endl;
		cout << "         --verbose     Show progress messages to STDERR" << endl;
		exit(0);
	}

	cout << "Usage:   " << argv[0] << " <command> [options]" << endl;
	cout << "Command: view        Dump SIM file to screen" << endl;
	cout << "         create      Create a SIM file from GTC files" << endl;
	cout << "         illuminus   Produce Illuminus output" << endl;
	cout << "         genosnp     Produce GenoSNP output" << endl;
	cout << "         qc          Produce QC metrics" << endl;
	cout << "         help        Display this help. Use 'help <command>' for more help" << endl;
	exit(0);
}

void loadManifest(Manifest *manifest, string manfile)
{
	if (manfile == "") throw("No manifest file specified");
	manifest->open(manfile);
}

void commandView(string infile, bool verbose)
{
	Sim *sim = new Sim();

	cout << endl << "Reading SIM file: " << infile << endl;
	sim->open(infile);
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

	char *sampleName = new char[sim->sampleNameSize];
	vector<uint16_t> *intensity_int = new vector<uint16_t>;
	vector<float> *intensity_float = new vector<float>;
	for (unsigned int n = 0; n < sim->numSamples; n++) {
		intensity_int->clear();
		intensity_float->clear();
		if (sim->numberFormat == 0) sim->getNextRecord(sampleName,intensity_float);
		else                        sim->getNextRecord(sampleName,intensity_int);
		cout << sampleName << "\t: ";
		if (verbose) {	// dump intensities as well as sample names
			// there *must* be a better way of doing this...
			if (sim->numberFormat == 0) {
				for (vector<float>::iterator i = intensity_float->begin(); i != intensity_float->end(); i++) {
					cout << *i << " ";
				}
			} else {
				for (vector<uint16_t>::iterator i = intensity_int->begin(); i != intensity_int->end(); i++) {
					cout << *i << " ";
				}
			}
		}
		cout << endl;
	}
	delete intensity_int;
	delete intensity_float;
	delete sampleName;
}

//
// Parse the infile, which is either an ascii list of GTC files, or a JSON format file
// Return an array of filenames and (for a JSON file) a list of sample names
//
void parseInfile(string infile, vector<string> &sampleNames, vector<string> &infiles)
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
void commandCreate(string infile, string outfile, bool normalize, string manfile, bool verbose)
{
	vector<string> sampleNames;		// list of sample names from JSON input file
	vector<string> infiles;			// list of GTC files to process
	Sim *sim = new Sim();
	Gtc *gtc = new Gtc();
	Manifest *manifest = new Manifest();
	int numberFormat = normalize ? 0 : 1;

	//
	// First, get a list of GTC files. and possibly sample names
	//
	if (infile == "") throw("commandCreate(): infile not specified");

	parseInfile(infile,sampleNames,infiles);
	if (infiles.size() == 0) throw("No GTC files are specified in the infile");

	// Let's check the GTC files, shall we?
	for (unsigned int n = 0; n < infiles.size(); n++) {
		gtc->open(infiles[n],0);
		if (gtc->errorMsg.length()) throw gtc->errorMsg;
	}

	// We need a manifest file to sort the SNPs and to normalise the intensities (if required)
	loadManifest(manifest, manfile);
	// Sort the SNPs into position order
	sort(manifest->snps.begin(), manifest->snps.end(), SortByPosition);

	// Create the SIM file and write the header
	sim->createFile(outfile);
	sim->writeHeader(infiles.size(),gtc->numSnps, 2, numberFormat);

	// For each GTC file, write the sample name and intensities to the SIM file
	for (unsigned int n = 0; n < infiles.size(); n++) {
		gtc->open(infiles[n], Gtc::XFORM | Gtc::INTENSITY);
		if (manifest->snps.size() != gtc->xRawIntensity.size()) {
			ostringstream msg;
			msg << "Size mismatch: Manifest contains " << manifest->snps.size() << " probes, but " 
			    << infiles[0] << " contains " << gtc->xRawIntensity.size() << " probes.";
			throw msg.str();
		}
		char *buffer = new char[sim->sampleNameSize];
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
			if (normalize) {
				// This is the normalization calculation, according to Illumina
				unsigned int norm = manifest->normIdMap[snp->normId];
				XFormClass *XF = &(gtc->XForm[norm]);
				double tempx = gtc->xRawIntensity[idx] - XF->xOffset;
				double tempy = gtc->yRawIntensity[idx] - XF->yOffset;
				double cos_theta = cos(XF->theta);
				double sin_theta = sin(XF->theta);
				double tempx2 = cos_theta * tempx + sin_theta * tempy;
				double tempy2 = -sin_theta * tempx + cos_theta * tempy;
				double tempx3 = tempx2 - XF->shear * tempy2;
				double tempy3 = tempy2;
				xn = tempx3 / XF->xScale;
				yn = tempy3 / XF->yScale;
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
void commandIlluminus(string infile, string outfile, string manfile, int start_pos, int end_pos, bool verbose)
{
	Sim *sim = new Sim();
	ofstream outFStream;
	ostream *outStream;
	char *sampleName;
    vector<uint16_t> *intensity_int = new vector<uint16_t>;
    vector<float> *intensity_float = new vector<float>;
	vector<vector<float> > SampleArray;
	Manifest *manifest = new Manifest();

	if (outfile == "-") {
		outStream = &cout;
	} else {
		outFStream.open(outfile.c_str(),ios::binary | ios::trunc | ios::out);
		outStream = &outFStream;
	}

	sim->open(infile);

	if (sim->numChannels != 2) throw("simtools can only handle SIM files with exactly 2 channels at present");

	sampleName = new char[sim->sampleNameSize];

	// We need a manifest file to sort the SNPs
	loadManifest(manifest, manfile);
	// Sort the SNPs into position order
	sort(manifest->snps.begin(), manifest->snps.end(), SortByPosition);

	if (end_pos == -1) end_pos = sim->numProbes - 1;

	// load the (relevant parts of) the SIM file
	if (verbose) cerr << "Reading SIM file " << infile << endl;
	*outStream << "SNP\tCoor\tAlleles";
	for(unsigned int n=0; n < sim->numSamples; n++) {
		vector<float> *s = new vector<float>; 
		if (!s) { cerr << "new s failed" << endl; exit(1); }
		intensity_float->clear();
		intensity_int->clear();
		if (sim->numberFormat == 0) sim->getNextRecord(sampleName, intensity_float);
		else                        sim->getNextRecord(sampleName, intensity_int);
		for (int i=start_pos; i <= end_pos; i++) {
			for (int c=0; c < sim->numChannels; c++) {
				float v;
				int k = i * sim->numChannels + c;
				if (sim->numberFormat==0) v = intensity_float->at(k);
				else                      v = intensity_int->at(k);
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

}

void commandGenoSNP(string infile, string outfile, string manfile, int start_pos, int end_pos, bool verbose)
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

	sim->open(infile);

	if (end_pos == -1) end_pos = sim->numSamples - 1;

	char *sampleName = new char[sim->sampleNameSize];
    vector<uint16_t> *intensity = new vector<uint16_t>;;
    for (int n=0; n <= end_pos ; n++) {
        intensity->clear();
        sim->getNextRecord(sampleName, intensity);
	if (n < start_pos) continue;
		*outStream << sampleName << "\t" << sampleName;
		for (vector<uint16_t>::iterator i = intensity->begin(); i != intensity->end(); i+=2) {
			*outStream << "\t" << std::fixed << setprecision(3) << *i;
			*outStream << " " << std::fixed << setprecision(3) << *(i+1);
		}
		*outStream << endl;
	}

	delete sim;
}

void commandQC(string infile, string magnitude, string xydiff, bool verbose)
{
  //Sim *sim = new Sim();
  cout << "QC functionality not yet operational!" << endl;

}

int main(int argc, char *argv[])
{
	string infile = "-";
	string outfile = "-";
	string manfile = "";
	string magnitude = "-";
	string xydiff = "-";
	bool verbose = false;
	bool normalize = false;
	int start_pos = 0;
	int end_pos = -1;
	int option_index = -1;
	int c;

	// Get command
	if (argc < 2) showUsage(argc,argv);
	string command = argv[1];
	if (command == "help") showUsage(argc, argv);

	// get options
	while ( (c=getopt_long(argc,argv,"h?",long_options,&option_index)) != -1) {
		if (c == '?') showUsage(argc,argv);	// unknown option
		if (option_index > -1) {
			string option = long_options[option_index].name;
			if (option == "infile") infile = optarg;
			if (option == "outfile") outfile = optarg;
			if (option == "verbose") verbose = true;
			if (option == "man_file") manfile = optarg;
			if (option == "man_dir") manfile = optarg;
			if (option == "normalize") normalize = true;
			if (option == "start") start_pos = atoi(optarg);
			if (option == "end") end_pos = atoi(optarg);
			if (option == "magnitude") magnitude = optarg;
			if (option == "xydiff") xydiff = optarg;
		}
	}

	// Allow stdin and stdout as synonyms for "-"
	if (infile == "stdin") infile = "-";
	if (outfile == "stdout") outfile = "-";

	// Process the command
	try {
	     if (command == "view")      commandView(infile, verbose);
	else if (command == "create")    commandCreate(infile, outfile, normalize, manfile, verbose);
	else if (command == "illuminus") commandIlluminus(infile, outfile, manfile, start_pos, end_pos, verbose);
	else if (command == "genosnp")   commandGenoSNP(infile, outfile, manfile, start_pos, end_pos, verbose);
	else if (command == "qc")        commandQC(infile, magnitude, xydiff, verbose);
	else {
		cerr << "Unknown command '" << command << "'" << endl;
		showUsage(argc,argv);
	}
	} catch (const char *error_msg) {
		cerr << error_msg << endl << endl;
		exit(1);
	}
	catch (string error_msg) {
		cerr << error_msg << endl << endl;
		exit(1);
	}
	return 0;
}


