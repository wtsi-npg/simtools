// 
// sim_tools.cpp
//
// Program to read and process contents of a SIM file
//
// $Id: sim.cpp 1354 2010-11-11 16:20:09Z js10 $
//
//
// Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>, Iain Bancarz <ib5@sanger.ac.uk>
//
// Copyright (c) 2009 - 2013 Genome Research Ltd. All rights reserved.
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


#include <getopt.h>
#include <iostream>
#include "commands.h"

using namespace std;

static struct option long_options[] = {
                   {"infile", 1, 0, 0},
                   {"outfile", 1, 0, 0},
                   {"man_file", 1, 0, 0},
                   {"man_dir", 1, 0, 0},
                   {"egt_file", 1, 0, 0},
                   {"normalize", 0, 0, 0},
                   {"verbose", 0, 0, 0},
                   {"start", 1, 0, 0},
                   {"end", 1, 0, 0},
                   {"magnitude", 1, 0, 0},
                   {"xydiff", 1, 0, 0},
                   {0, 0, 0, 0}
               };

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
          cout << "         --man_file <dirname>   Directory to look for Manifest file in" << endl;
          cout << "         --normalize            Normalize the intensities (default is raw values)" << endl;
          cout << "         --verbose              Show progress messages to STDERR" << endl;
          exit(0);
	}

        if (command == "fcr") {
          cout << "Usage:   " << argv[0] << " fcr [options]" << endl << endl;
          cout << "Create a FCR file from a list of GTC files" << endl<< endl;
          cout << "Options: --infile <filename>    File containing list of GTC files to process" << endl;
          cout << "         --outfile <filename>   Name of FCR file to create or '-' for STDOUT" << endl;
          cout << "         --man_file <dirname>   Path to bpm.csv manifest file" << endl;
          cout << "         --egt_file <dirname>   Path to EGT binary cluster file" << endl;
          cout << "         --verbose              Show progress messages to STDERR" << endl;
          exit(0);
        }

	if (command == "illuminus") {
          cout << "Usage:   " << argv[0] << " illuminus [options]" << endl << endl;
          cout << "Create an Illuminus file from a SIM file" << endl<< endl;
          cout << "Options: --infile <filename>    Name of SIM file to process or '-' for STDIN" << endl;
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
          cout << "Options: --infile      The name of the SIM file (cannot accept STDIN)" << endl;
          cout << "         --magnitude   Output file for sample magnitude (normalised by SNP); cannot use STDOUT" << endl;
          cout << "         --xydiff      Output file for XY intensity difference; cannot use STDOUT" << endl;
          cout << "         --verbose     Show progress messages to STDERR" << endl;
          exit(0);
	}

	cout << "Usage:   " << argv[0] << " <command> [options]" << endl;
	cout << "Command: view        Dump SIM file to screen" << endl;
	cout << "         create      Create a SIM file from GTC files" << endl;
	cout << "         fcr         Create a FCR file from GTC files" << endl;
	cout << "         illuminus   Produce Illuminus output" << endl;
	cout << "         genosnp     Produce GenoSNP output" << endl;
	cout << "         qc          Produce QC metrics" << endl;
	cout << "         help        Display this help. Use 'help <command>' for more help" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{
	string infile = "-";
	string outfile = "-";
	string manfile = "";
        string egtfile = "";
	string magnitude = "";
	string xydiff = "";
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
			if (option == "egt_file") egtfile = optarg;
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

	Commander *commander = new Commander();
	// Process the command
	try {
	  if (command == "view") {
	    commander->commandView(infile, verbose);
	  } else if (command == "create") {
	    commander->commandCreate(infile, outfile, normalize, 
				     manfile, verbose);
	  } else if (command == "fcr") {
            commander->commandFCR(infile, outfile, manfile, egtfile, verbose);
          } else if (command == "illuminus") {
	    commander->commandIlluminus(infile, outfile, manfile, 
					start_pos, end_pos, verbose);
	  } else if (command == "genosnp") {
	    commander->commandGenoSNP(infile, outfile, manfile, 
				      start_pos, end_pos, verbose);
	  } else if (command == "qc") {
	    commander->commandQC(infile, magnitude, xydiff, verbose);
	  } else {
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
	delete commander;
	return 0;
}


