// b2g.cpp
//
// Matthew Gillman, WTSI, 13th October 2010.
//
// $Id: b2g.cpp 1331 2010-10-25 16:07:24Z mg10 $
//
// b2g class implementation.
//


#include "b2g.h"
#include <fstream>

using namespace std;



/////////////////////////////////////////////////////////////////
//
// void b2g::initialise_colpos ()
//
// Sets all _required_ values in colpos to -1.
//
// NB it is important that other, non-required values are NOT initialised.
//
/////////////////////////////////////////////////////////////////

void b2g::initialise_colpos () {

  colpos.clear();

  colpos[SNP_name] = -1;
  colpos[SampleID] = -1;
  colpos[XRaw] = -1;
  colpos[YRaw] = -1;

}



/////////////////////////////////////////////////////////////////
//
// void b2g::write_SNPs_file (...)
//
// final: the name of the SNPs file to produce
//
// Write out SNPs file for input to GenoSNP.
// We just write them out in the same order in which they are stored
// in Manifest class's vector. 
//
/////////////////////////////////////////////////////////////////

void b2g::write_SNPs_file (std::string& final) {

  // Create temp file, write temp file, then copy file across to final resting place and
  // check its hashvalue to make sure it was not corrupted duing the copy! Then delete tempfile.
  
  string tempname = "/tmp/b2i_snpsfile_temp_XXXXXX";
 
  // This gets name of temp file and also opens it:
  char* tempfile = create_tempfile (string("write_GenoSNP_SNPs_file"), tempname);

  // Associate file stream with the file.
  FILE* pFile = fopen (tempfile, "w");
  
  if ( NULL == pFile ) {
    cout << "\n\nError in opening output file. \n\n" << flush;
    exit(1);
  }

  int numprobes = mf.snps.size();

  for ( int i = 0; i < numprobes; i++ ) {
	
	snpClass s = mf.snps[i];
	std::string name = s.name;
	
	// Now get the properties for this probe.

	if ( ! s.position || ! s.snp ) {
	  cout << "\nWARNING: insufficient data found in list for SNP " 
		   << name << "." << flush;
	  continue;
	}

	// OK, we have data; let's output it. Probe details first:
	int beadsetID = s.BeadSetID;
	if ( -1 == beadsetID ) {
	  beadsetID = (s.normId % 100) + 1;
	}

	if ( 0 > fprintf (pFile, "%s\t%i\t%c %c\n", name.c_str(), beadsetID,
					  s.snp[0], s.snp[1]) ) {
	  cout << "\nOutput error!\n" << flush;
	  exit(1);
	}

  }

  old_copy_file((const char* const) tempfile, (const char* const) final.c_str());

  delete[] tempfile;

}



/////////////////////////////////////////////////////////////////
//
// FILE* b2g::create_samples_file (string& tempname)
//
// Creates temp file for GenoSNP sample results. Stores the name
// of the temp file created in tempname. returns the pointer to
// the file stream.
//
/////////////////////////////////////////////////////////////////

FILE* b2g::create_samples_file ( string& tempname ) {

  string mytempname = "/tmp/b2i_samplesfile_temp_XXXXXX";
		   
  char* tempfile = create_tempfile (string("write_GenoSNP_Samples_file"), mytempname);
		   
  tempname = tempfile;

  // Associate file stream with the file.
  FILE* pFile = fopen (tempfile, "w");
  
  if ( NULL == pFile) {
	cout << "\n\nError in opening output file.\n\n" << flush;
	exit(1);
  }


  return pFile;

}



/////////////////////////////////////////////////////////////////
//
// void b2g::process_input_file (...)
//
// Extracts data from a FCR file to be processed.
//
// input: name of the file to process.
// outFile: pointer to filestream to the output file.
//
/////////////////////////////////////////////////////////////////

void b2g::process_input_file ( std::string input, 
							   FILE* outFile ) {


  // Reset stored column numbers in case we are processing multiple input files
  // (they could vary between files).
  initialise_colpos();

  ifstream filestr;
  char* in_string;
  try {
	in_string = new char[input.size() + 1];
  }
  catch (bad_alloc& ba) {
    cerr << "\nERROR: b2g::process_input_file(), in_string, caught bad_alloc: "
		 << ba.what() << "\n" << flush;
    exit(1);
  }
  strcpy(in_string, input.c_str());
  filestr.open ( in_string, fstream::in );
  delete[] in_string;
  
  if ( filestr.fail() ) {
    cout << "\n\nError in opening input file " << input << ".\n\n" << flush;
    exit(1);
  }
  
  
  // Skip over header block (keep looping until "SNP Name" row found).
  // good() checks that none of the bits eofbit, failbit or badbit are set.
  // Note that this means it checks no EOF encountered yet.
  string str;
  getline (filestr, str);
  
  if ( ! filestr.good() ) {
	cout << "\nI/O file error in b2g::process_input_file() (loop top).\n" 
		 << flush;
	exit(1);
  }

  while ( string::npos == (size_t) str.find ("SNP Name") ) { // Equal until found.
    //cout << "\nRead:" << str << "\n";
    getline (filestr, str);
	if ( ! filestr.good() ) {
	  cout << "\nI/O file error in b2g::process_input_file() (loop body).\n" 
		   << flush;
	  exit(1);
	}
  }
  
  // Remove any nasty Windows carriage return characters at end of string:
  size_t where = str.find_last_of ('\r');
  if ( where != str.npos ) {
    str.erase (where, 1);
  }
  
  // We now have the line containing the column titles for the data. We need to find
  // which columns have the required data in. In an attempt to be robust to any
  // future changes in input file format, e.g. addition or deletion of columns, 
  // rather than hard-coding in the required column numbers we locate them:
  
  std::string delimiter = "\t";
  
  if ( ! find_column_numbers(str, delimiter) ) {
    
    delimiter = ",";
    
    if ( ! find_column_numbers(str, delimiter) ) {
      cout << "\n\nERROR I haven't found all required columns. Either they or not all present or the input file has a delimiter which I don't know about. Goodbye...\n\n" << flush;
      exit(1);
    }

  }
  
  char delim = delimiter.at(0); // Note that this assumes the delimiter is a single char.
  cout << "\nOK, the input file is delimited by " << delimiter << "." << flush;



  string pile_of_junk = "this is a pile of junk!";  
  string current_sample = pile_of_junk;


  //
  // OK, we can now cycle through the input file and read in the data values...
  //

  // Map to store current sample's values:
  map<string, string> mymap; // Key = probe name; value = readings.


  while ( getline(filestr, str) ) {

    // cout << "\nDEBUG: Read:" << str << "\n";



    if ( ! filestr.good() ) {

	  cout << "\nLast line read in input file was:\n" << str << "\n" << flush;


      // Check. If the eof bit has been set this is fine,
      // as long as it's the last line in the file!
	  // If one or both of the other bits have been set, that's bad.

	  if ( filestr.fail() ) {
		cout << "\nERROR: badbit and/or failbit were set\n" << flush;
		if ( filestr.bad() ) {
		  cout << "\nThe bad bit was definitely set. I don't know about the failbit. \n" 
               << flush;
		}
		exit(1);		
	  }

	  else {
		if ( filestr.eof() ) {
		  cout << "\nReached end of file." << flush;

		  // Clear error state
		  filestr.clear();
		}
		else {
		  cout << "\nI'm in b2g::process_input_file(), but in an impossible place!\n" 
			   << flush;
		  exit(1);
		}

	  }

    } // End if ( ! filestr.good() )


    // Remove any nasty Windows carriage return characters at end of string:
    size_t where = str.find_last_of ('\r');
    if ( where != str.npos ) {
      str.erase (where, 1);
    }

    vector<string> line_things;
    char *line = strdup(str.c_str());
	if ( NULL == line ) {
	  cout << "\nERROR: unable to strdup() memory in read_input_file()\n" 
		   << flush;
	  exit(1);
	}
     
    int length = strlen(line);
    //cout <<"\n\nlength = " << length << "\n";
    int delimnum = 0;

    // For this line, locate and store the positions of each delimiter:
    unsigned int* positions;
	try {
	  positions = new unsigned int [LAST_COL_NUM + 1];
	}
   	catch (bad_alloc& ba) {
      cerr << "\nERROR: in b2g::process_input_file(), positions, caught bad_alloc: " 
		   << ba.what() << "\n" << flush;
	  exit(1);
	}
    for (int pos = 0; pos < length; pos++ ) {      
      if ( line[pos] == delim ) {
		//cout << "\n at " << pos;
        *(positions + delimnum) = pos;
        delimnum++;
      }      
    }
	// Now delimnum stores the number of entries in the positions array.

	//cout << "\nDEBUG: delimnum = " << delimnum << endl << flush;

    // Now pull the info out of "line":
    
    std::string a_snpname, a_sample, a_rawX, a_rawY;
    
    extract_datum_from_string (positions, delimnum, line, 
							   length, a_snpname, colpos[SNP_name]);
    extract_datum_from_string (positions, delimnum, line, 
							   length, a_sample, colpos[SampleID]);
	extract_datum_from_string (positions, delimnum, line, 
							   length, a_rawX, colpos[XRaw]);
	extract_datum_from_string (positions, delimnum, line, 
							   length, a_rawY, colpos[YRaw]);
    	    
    delete[] positions;
 
    free(line);


	// If the sample ID has changed, we need to write out
    // the data of (what was) the current sample:

	if ( a_sample != current_sample ) {

	  if ( current_sample != pile_of_junk ) {

		write_line_to_samples_file (/*mf,*/ outFile, mymap, current_sample);
		mymap.clear();

	  }

      // As we have changed samples, update things:
	  current_sample = a_sample;
	 
	  // cout << "\nCurrent sample is " << current_sample << "\n" << flush;

    } // End if ( a_sample != current_sample )


    // Store the data for this probe:
    mymap[a_snpname] = a_rawX + " " + a_rawY;


  } // End while ( getline(filestr, str) )

  // Now write out the data for the final sample in the file:
  write_line_to_samples_file (/*mf,*/ outFile, mymap, current_sample);
  mymap.clear();

  filestr.close();


} // End of read_input_file()



/////////////////////////////////////////////////////////////////
//
// void b2g::write_line_to_samples_file(...)
//
// Writes one sample's one-and-only data line to the Samples file.
//
// outFile: pointer to output filestream
// mymap: holds data for this sample (key = SNP name, value = data)
// sampleID: the ID of this sample.
//
/////////////////////////////////////////////////////////////////

void b2g::write_line_to_samples_file ( FILE* outFile,
									   map<string, string>& mymap,
									   string& sampleID ) {

  // Dump out the data of the old current_sample.
  // Manifest's vector of snps is now in ascending position order.
  
  if ( 0 > fprintf (outFile, "%s\t%s", sampleID.c_str(), sampleID.c_str()) ) {
	cout << "\nERROR writing sampleID to file for sample "
		 << sampleID <<"\n" << flush;
	exit(1);
  }
  
  
  int siz = mf.snps.size();
  
  //cout << "\nb2g DBG: manifest size = " << siz << endl << flush;
  
  for (int i = 0; i < siz; i++) {
	snpClass s = mf.snps[i];
	string snpName = s.name;
	string data;
	
	// Do not use map's [] operator for access; it creates a blank
	// entry if the data are not present in the map.
	map<string, string>::iterator it;
	it = mymap.find(snpName);
	if ( it == mymap.end() ) {
	  data = "NaN NaN";
	}
	else {
	  data = it->second;
	}
	
	if ( 0 > fprintf (outFile, "\t%s", data.c_str()) ) {
	  cout << "\nERROR writing data to file for sample "
		   << sampleID <<"\n" << flush;
	  exit(1);
	}
	
  } // End (for)
  
  if ( 0 > fprintf (outFile, "\n") ) {
	cout << "\nERROR writing newline to file for sample "
		 << sampleID <<"\n" << flush;
	exit(1);
  }
  

}



/////////////////////////////////////////////////////////////////
//
// void b2g::usage ()
//
// Prints help message.
//
/////////////////////////////////////////////////////////////////

void b2g::usage () {

  cout << "\nb2g: produce GenoSNP input files from BeadStudio FCR file.";

  cout << "\n\nUsage:"
	   << "\n b2g -m manifest [-i input_file | -l list_of_files] -s snps -d data [-w] [-k] "
	   << "\nwhere:"
	   << "\nmanifest (aka snps file): list of probes in FCR file(s)."
	   << "\ninput_file: a single BeadStudio Final Call Report (FCR) file."
	   << "\nlist_of_files: a file listing a number of such files (1 per line)"
	   << "\nsnps: name of GenoSNP SNPs file produced by b2g"
	   << "\ndata: name of GenoSNP Samples file produced by b2g"
	   << "\n\nUsually, manifest files will have the extension .bpm.csv and"
	   << "\nhave 9 columns. The -w option allows the \"wide format\" manifest"
	   << "\nfile to be used instead. This has many more fields; it also"
	   << "\nstores BeadSetID values, so should be more accurate. Recommended."
	   << "\n\n(Not recommended) The -k option keeps CNV probes (default: remove).";

  cout << "\n\nDo not use spaces in filenames.";

  cout << "\n\n";

}



/////////////////////////////////////////////////////////////////
//
// void b2g::process_command_line(...)
//
// Deals with the optins specified at program start.
//
// argc, argv: passed from main().
// strings: needed for initialisation for main()
// wf: this function sets this to true if using wide-format snp file
//
/////////////////////////////////////////////////////////////////

void b2g::process_command_line (int argc, char** argv,
								string& manifest_file,
								string& input,
								string& listing,
								string& snps,
								string& data,
								bool& wf) {

  char c;

  bool keep_cnvs = false;

  while ( (c = getopt (argc, argv, "m:i:l:s:d:wkh")) != -1 ) {

    switch (c) {
	case 'm':	
	  manifest_file = optarg; 
	  break;
	case 'i':	// Single input file.
	  input = optarg;
	  break;
	case 'l':   // List of input file(s).
	  listing = optarg;
	  break;
	case 's':
	  snps = optarg;
	  break;
	case 'd':
	  data = optarg;
	  break;
	case 'w':
	  wf = true;
	  break;
	case 'k':
	  keep_cnvs = true;
	  break;
	case 'h':
	  usage();
	  exit(1);
	  break;
	}

  }

  if ( manifest_file.empty() ) {
	cout << "\nPlease specify a manifest file.\n";
	exit(1);
  }

  if ( input.empty() && listing.empty() ) {
	cout << "\nPlease specify either a single input file or a file listing one or more files.\n";
	exit(1);
  }

  if ( ( ! input.empty() ) && ( ! listing.empty() ) ) {
	cout << "\nPlease use EITHER -i OR -l option; not both.\n";
	exit(1);
  }

  if ( snps.empty() ) {
	cout << "\nPlease specify name of GenoSNP SNPs file to produce.\n";
	exit(1);
  }

  if ( data.empty() ) {
	cout << "\nPlease specify name of GenoSNP Samples file to produce.\n";
	exit(1);
  }

  if ( ! keep_cnvs ) {
	mf.exclude_cnvs();
  }

}



/////////////////////////////////////////////////////////////////
//
// main()
//
/////////////////////////////////////////////////////////////////

int main (int argc, char** argv) {

  
  string manifest_file, inputfile, listing, output_snps_file, final_samples_file;
  bool using_wide_format = false;

  b2g beatty;

  beatty.process_command_line(argc, argv, manifest_file, inputfile, listing, output_snps_file, final_samples_file, using_wide_format);

  if ( using_wide_format ) {
	cout << "\nUsing wide format probe file.\n" << flush;
  }


  beatty.mf.open(manifest_file, using_wide_format);
  beatty.mf.order_by_position();

  beatty.write_SNPs_file (output_snps_file);


  //
  // Create samples file:
  //
  string tempsamplefile;
  FILE* pFile = beatty.create_samples_file(tempsamplefile);


  //
  // Get list of file(s) which we are to process.
  //
  vector<string> filenames;

  if ( ! inputfile.empty() ) { // Single file to process.
	filenames.push_back(inputfile);
  }
  else { // Process list of files.
	ifstream filestr;
	filestr.open (listing.c_str(), fstream::in);
	if ( filestr.fail() ) {
      cout << "\n\nError in opening input list file " 
		   << listing << endl << flush;
      exit(1);
    }
	string line_string;
	while ( getline(filestr, line_string) ) {

	  // Remove any whitespace characters before, in, or after filename:
	  string no_ws_line_string;
	  for (unsigned int charpos = 0; charpos < line_string.length(); charpos++)
      {
        char candidate = line_string[charpos];
        if ( ! isspace(candidate) ) {
          no_ws_line_string.append(1, candidate);
        }
      }
	  
	  filenames.push_back(no_ws_line_string);
	}

  }


  //
  // OK, now process each file in our vector:
  //
  int num = filenames.size();
  for ( int i = 0; i < num; i++ ) {
	string fname = filenames[i];
	cout << "\nAbout to process " << fname << endl << flush;
	beatty.process_input_file (fname, pFile);
  }

  fclose (pFile);

  // Now copy temp file across:
   if ( ! beatty.old_copy_file(tempsamplefile.c_str(), final_samples_file.c_str()) ) {
	 cout << "\nERROR in copying temp file " << tempsamplefile << " to "
		  << final_samples_file << "\n" << flush;
	 exit(1);
   }
  

  // Do integrity checks:

  


  return 0;


}


/*
 chr1:109457160
 chr1:109457233
 chr1:109457614
 chr1:109457618
 chr1:109457943
 chr1:109458224
 chr1:109458469
 chr1:109458772
 chr1:109459424
 chr1:109460430

200003	58805_E09_PROMIS375810	A	A	0.7593	0.032	1.224	1.164	0.059	15568	1652*
OK, the input file is delimited by 	.

 at 6
 at 29
 at 31
 at 33
 at 40
 at 46
 at 52
 at 58
 at 64
 at 70


  // Testing
  // manifest_file = "test.bpm.csv"; //"/nfs/wtccc/data5/genotyping/chip_descriptions/Illumina_Infinium/Human_Cardio-Metabo/Cardio-Metabo_Chip_11395247_A.bpm.csv";
  // input = fcr.test
  //output_snps_file = "test.snps";
  //final_samples_file = "test.samples";

./b2g -m test.bpm.csv -i  fcr.test -s test.snps -d test.samples

 ./b2g -i /nfs/new_illumina_geno05/FINISH_HAP_MAP/FINISH_HAP_MAP_FinalReport_LRR_BAF.txt -s fins.snps -d fins.samples -m /nfs/wtccc/data5/genotyping/chip_descriptions/Illumina_Infinium/Human_1-2M-2/Human1-2M-DuoCustom_v1_A.bpm.csv

*/
