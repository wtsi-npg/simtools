// b2i.cpp
//
// $Id: b2i.cpp 1380 2010-11-25 15:43:17Z mg10 $
//
// Convert BeadStudio output files into Illuminus input intensity files.
// 
// b2i class implementation.
//
// Matthew Gillman, WTSI, 16th March 2009


#include "b2i.h"


/*

TODO

      - delete temp files at exit unless emergency
	  
	  - add atexit() to clear out temp files.

*/



using namespace std;


// Set all required values to -1.
void b2i::initialise_colpos () {

  colpos.clear();

  colpos[SNP_name] = -1;
  colpos[SampleID] = -1;
  colpos[Xnorm] = -1;
  colpos[Ynorm] = -1;

}



/////////////////////////////////////////////////////////////////
//
// void b2i::process_command_line (int argc, char** argv, std::string& input,
//                       std::string& output, std::string& snpfile)
//
/////////////////////////////////////////////////////////////////

void b2i::process_command_line (int argc, char** argv,
						   std::string& input,
						   std::string& listing,
						   std::string& output,
						   std::string& snp_file)
  
{
  bool gotP = false;
  bool gotD = false;
  bool gotC = false;
  bool gotF = false;
  char c;
  
  while ((c = getopt (argc, argv, "s:i:l:o:ghnpdc:m:f:")) != -1) {
    switch (c) {
	  case 's':	
		snp_file = optarg; 
		break;
	  case 'i':	// Processing a single BeadStudio file as input.
		input = optarg; 
		break;
	  case 'l': // Processing a file which lists a number of BS input files
		listing = optarg;
		break;
      case 'o':
		output = optarg;
		break;
      case 'h':
		usage();
		break;
      case 'p':
        SORT_BY_POS = true;
		gotP = true;
        break;
      case 'd':
        SORT_BY_POS = false;
		gotD = true;
        break;
	  case 'c':
	  {
		string mychrom = optarg;
	
		if ( 0 != mychrom.compare("all") ) {
		  
		  // Special case for mitochondrials.
		  // Whether user has specified "M", "MT" or "Mt", store as "MT".
		  if ( ( 0 == mychrom.compare ("M") )
			   ||
			   ( 0 == mychrom.compare ("Mt") ) )
		  {
			  SELECTED_CHROM = "MT";
		  }
		  else 
		  {
			SELECTED_CHROM = optarg;
		  }
		  
		  // Check for non-existent chromosomes!
		  char* chkstr;
          try {
            chkstr = new char [SELECTED_CHROM.size() + 1];
		  }
          catch (bad_alloc& ba) {
            cerr << "\nERROR: in process_command_line(), chkstr, caught bad_alloc: " 
                 << ba.what() << "\n" << flush;
            exit(1);
          }
		  strcpy (chkstr, SELECTED_CHROM.c_str());
		  int mycheck = atoi (chkstr);
		  delete[] chkstr;
		  bool ok = false;
		  if ( 0 != mycheck ) // It's an integer
		  { 
			if ( ( mycheck >= 1 ) && ( mycheck <= 22 ) ) 
			{
			  ok = true;
			}
		  }
		  else // Not an integer
		  { 
			if ( ( 0 == SELECTED_CHROM.compare("MT") ) ||
				 ( 0 == SELECTED_CHROM.compare("X") ) ||
				 ( 0 == SELECTED_CHROM.compare("XY") ) ||
				 ( 0 == SELECTED_CHROM.compare("Y") ) )  
			{
			  ok = true;
			}
		  }
		  
		  if ( false == ok ) {
			cout << "\nError. Chromosomes should be 1-22, X, Y, XY or M/MT/Mt (or \"all\")." << flush;
			exit(1);
		  }
		  
		} // End (if)
		gotC = true;
		break;
	  } // End (case 'c')
		  
	case 'f': // force b2i to use this chromosome (for future compatability)
	{
	  SELECTED_CHROM = optarg;
	  gotF = true;
	  break;
	}
	  
	case 'm': // memory estimation
	  string quantities = optarg; //format L:N (snps:samples)
	  size_t colon_pos = quantities.find_first_of(":");
	  string L_string = quantities.substr(0, colon_pos);
	  string N_string = quantities.substr(colon_pos + 1, quantities.size() - (colon_pos + 1)); 
	  cout << "\n\nEstimating memory usage with L = " << L_string << " probes and N = " 
		   << N_string << " samples.";

	  char* L_chars;
	  try {
		L_chars = new char [L_string.size() + 1];
	  }
	  catch (bad_alloc& ba) {
	    cerr << "\nERROR: in process_command_line(), L_chars, caught bad_alloc: " 
             << ba.what() << "\n" << flush;
        exit(1);
      }
	  strcpy (L_chars, L_string.c_str());
	  long L = atol(L_chars);
	  delete[] L_chars;
	  
	  char* N_chars;
	  try {
		N_chars = new char [N_string.size() + 1];
	  }
	  catch (bad_alloc& ba) {
		cerr << "\nERROR: in process-command_line(), N_chars, caught bad_alloc: " 
             << ba.what() << "\n" << flush;
		exit(1);
	  }
	  strcpy (N_chars, N_string.c_str());
	  long N = atol(N_chars);
	  delete[] N_chars;
	  
	  float /*long*/ mem = (188 * L) + ( 96 * L * N) + (52 * N);
	  cout << "\nI think you'll need about " << mem << " bytes, or about ";
	  float mega = mem / (1024*1024);
	  cout << mega << " MB.";
	  
	  float cache_factor = 25; // Additional percent required
	  cout << "\nAdding an extra " << cache_factor << "% for cache memory brings this total to ";
	  mem = mem * ( 1 + (cache_factor/100) );
	  mega = mem / (1024*1024);
	  cout << mem << " bytes, or about " << mega << " MB.\n\n";
	  
	  exit (1);
	  break;
	
	}
	
  } // End (while)


  cout << "\n\n*** PRODUCING OUTPUT FOR ILLUMINUS INPUT ***\n" << flush;
  

  if ( gotC && gotF ) {
	cout << "\n Do not use -c and -f options together."
		 << "\n(Use -h option to see usage.)\n";
    exit(1);
  }
  

  if ( gotP && gotD ) {
	cout << "\nPlease decide: do you want to sort by probe position or not?"
		 << "\n(Use -h option to see usage.)\n";
	exit(1);
  }

  if ( 0 == snp_file.size() ) {
    cout << "\n\n ERROR: No SNP (probe) file specified" << endl;
    usage();
  }

  if ( ( 0 == input.size() ) && ( 0 == listing.size() ) ) {
    cout << "\n\n ERROR: Either a single input file, or a file listing multiple such, must be specified." << endl;
    usage();
  }

  if ( ( 0 != input.size() ) && ( 0 != listing.size() ) ) {
    cout << "\n\n ERROR: Don't use -i and -l options together." << endl;
    usage();
  }

  if ( 0 == output.size() ) {
    cout << "\n\n ERROR: No output file specified." << endl;
    usage();
  }
  
  cout << "\n\nExecuting with:"
       << "\nSNP (probe) file: " << snp_file << flush;

  if ( input.size() > 0 ) {
    cout << "\nInput file: " << input << flush;
  }
  else if ( listing.size() > 0 ) {
    cout << "\nList of input files: " << listing << flush;
  }
  
  cout << "\nOutput file: " << output
       << "\nProcessing: " << flush;


  if ( SELECTED_CHROM.size() > 0 ) {
    cout << "chromosome " << SELECTED_CHROM << " only." << flush;
  }
  else {
    cout << "all chromosomes." << flush;
  }

  if ( SORT_BY_POS ) {
    cout << "\nOutput file will be sorted in ascending probe position order." << flush;
  }

  cout << "\nA temporary file will be created to try and solve I/O problems." << flush;
  
  //cout << "\n";
  

} // End of b2i::process_command_line()



/////////////////////////////////////////////////////////////////
//
// b2i::usage()
//
// Prints usage information, then exits.
//
/////////////////////////////////////////////////////////////////

void b2i::usage () {

  cout << "\nb2i: produce files from BeadStudio format in Illuminus input format.";

  cout << "\n\nUsage:"
       << "\n b2i -s snpfile [-i input_file | -l list_of_files] -o output_file \\"
       << "\n    [-c chrom | -f chrom] [-p|-d]"
	   << "\n OR   b2i -h    OR   b2i -m L:N"
       << "\n* snpfile is file containing SNP (& other probe) positions."
       << "\n* input_file is a single BeadStudio file to be processed."
       << "\n* list_of_files is a file listing multiple BeadStudio files. This file must consist"
	   << "\n    only of the absolute path to one such file per line. There can be whitespace"
       << "\n    around each such file path, but not IN any of the file paths."
       << "\n* output_file is name of Illuminus output file produced."
       << "\n* the -h option prints this usage message."
       << "\n* chrom is either a chromosome or \"all\" for all chromosomes."
       << "\n* OPTIONAL the -p option means that the output file will be sorted in ascending"
	   << "\n    probe position order. (This is now the default and is retained for"
	   << "\n    backwards-compatibility.)"
       << "\n* OPTIONAL the -d option means DON'T sort by ascending probe position"
	   << "\n    (i.e. sort probes by name)."
       << "\n\nIf the -c option \"all\" is used (rare), all calls across all chromosomes"
	   << "\n    will be processed. Otherwise (usual case), only the single"
	   << "\n    chromosome specified will be processed. Valid values are the following:"
       << "\n    1-22, X, Y, XY. Mitochondrials can be specified with  M, MT or Mt.\n";


  cout << "\nThe -f (force b2i to use this) option can be used to process individual chromosomes"
	   << "\nonly, as per the non-\"all\" -c options. However, no checking is performed of the input."
	   << "\nThis is to future-proof b2i; if a snpfile is encountered with chromosome names not"
	   << "\nmatching the options enforced by -c, the -f option should be used.\n\n"; 

  cout << "The -m option estimates the *theoretical minimum* memory requirement,"
       << "\ngiven L probes and N samples. If you have specified a single chromosome (normal"
       << "\ncase), L is the number of probes for this chromosome in the SNP file. It's"
       << "\nbest to use AT LEAST the estimated memory + 25% for caching, etc.\n\n";  

  //  cout << "\n\n NB snpfile's format must be compatible with that required "
  //       << "\n by the Manifest class.\n\n";

  exit(1);

} // End of b2i::usage().



/////////////////////////////////////////////////////////////////
//
// void read_input_file (...)
//
// Extracts data from the file to be converted
//
/////////////////////////////////////////////////////////////////

void b2i::read_input_file (std::string input, 
                      map < std::string, hash_map <std::string, std::string> >& readings,
                      map <std::string, int>& samples,
                      Manifest& mf) {

  // cout << "\n DEBUG: read_input_file() processing file " << input << "\n";

  // Reset stored column numbers in case we are processing multiple input files
  // (they could vary between files).
  initialise_colpos();

  ifstream filestr;
  char* in_string;
  try {
	in_string = new char[input.size() + 1];
  }
  catch (bad_alloc& ba) {
    cerr << "\nERROR: in read_input_file(), in_string, caught bad_alloc: " << ba.what() << "\n" << flush;
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
	cout << "\nERROR: I/O file error in read_input_file() (loop top).\n" << flush;
	exit(1);
  }

  while ( string::npos == (size_t) str.find ("SNP Name") ) { // Equal until found.
    //cout << "\nRead:" << str << "\n";
    getline (filestr, str);
	if ( ! filestr.good() ) {
	  cout << "\nERROR: I/O file error in read_input_file() (loop body).\n" << flush;
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
      cout << "\n\nERROR I haven't found all required columns.\n"
		   << "Either they are not all present or the input file has\n"
		   << "a delimiter which I don't know about. Goodbye...\n\n" << flush;
      exit(1);
    }

  }
  
  char delim = delimiter.at(0); // Note that this assumes the delimiter is a single char.
  cout << "\nOK, the input file is delimited by " << delimiter << "." << flush;


  //
  // OK, we can now cycle through the input file and read in the data values...
  //

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
		  cout << "\nERROR: I'm in read_input_file(), but in an impossible place!\n" << flush;
		  exit(1);
		}

	  }

    } // End if ( ! filestr.good() )

	/*
	// Final getline() will not have entered loop above but probably set an error flag.
	if ( ! filestr.good() ) {
	if ( filestr.fail() ) {
	cout << "\nERROR mmm: badbit and/or failbit were set\n" << flush;
	if ( filestr.bad() ) {
	cout << "\nmmmThe bad bit was definitely set. I don't know about the failbit. \n" << flush;
	}
	exit(1);		
	}
	else {
	if ( filestr.eof() ) {
	cout << "\nmmmReached end of file." << flush;
	
	// Clear error state
	filestr.clear();
	}
	else {
	cout << "\nmmmERROR: I'm in read_input_file(), but in an impossible place!\n" << flush;
	exit(1);
	}
	}
	}
	*/


    // Remove any nasty Windows carriage return characters at end of string:
    size_t where = str.find_last_of ('\r');
    if ( where != str.npos ) {
      str.erase (where, 1);
    }

    vector<string> line_things;
    char *line = strdup(str.c_str());
	if ( NULL == line ) {
	  cout << "\nERROR: unable to strdup() memory in read_input_file()\n" << flush;
	  exit(1);
	}
     
    int length = strlen(line);
    // cout <<"\n\nlength = " << length << "\n";
    int delimnum = 0;

    // For this line, locate and store the positions of each delimiter:
    unsigned int* positions;
	try {
	  positions = new unsigned int [LAST_COL_NUM + 1];
	}
   	catch (bad_alloc& ba) {
      cerr << "\nERROR: in read_input_file(), positions, caught bad_alloc: " 
		   << ba.what() << "\n" << flush;
	  exit(1);
	}
    for (int pos = 0; pos < length; pos++ ) {      
      if ( line[pos] == delim ) {
        // cout << "\n at " << pos;
        *(positions + delimnum) = pos;
        delimnum++;
      }      
    }
	// Now delimnum stores the number of entries in positions array.

    // Now pull the info out of "line":
    
    std::string snpname, sample, normX, normY, rawX, rawY; 

    extract_datum_from_string (positions, delimnum, line, length, snpname, colpos[SNP_name]);


    // Filter out non-selected chromosomes, if required: 
    if ( SELECTED_CHROM.size() > 0 ) {
      if ( 0 != SELECTED_CHROM.compare(mf.get_chromosome_for_SNP(snpname)) ) {
		delete[] positions;
		free(line);
		continue;
      }
    }

    extract_datum_from_string (positions, delimnum, line, length, sample, colpos[SampleID]);
	extract_datum_from_string (positions, delimnum, line, length, normX, colpos[Xnorm]);
	extract_datum_from_string (positions, delimnum, line, length, normY, colpos[Ynorm]);

    delete[] positions;
    free(line);


    //cout << "\nSNP name is " << snpname << "\n";
    //cout << "\n value of norm X is " << normX << endl;

    // Update data structures
    samples[sample] = 1;
	//cout << "\nNo. of instances of this sample = " <<   samples[sample] << "\n";
  
	std::string XYints = normX + "," + normY;

    readings[snpname][sample] = XYints;

  } // End while ( getline(filestr, str) )


  filestr.close();

  /*
  // Check closure was clean.
  if ( filestr.fail() ) {
  cout << "\nERROR in read_input_file(): failbit and/or badbit set when file was closed\n" << flush;
  if ( filestr.bad() ) {
  cout << "\nbadbit was definitely set. I don't know about failbit\n" << flush;
  }
  exit(1);
  }
  */

} // End of read_input_file()









/////////////////////////////////////////////////////////////////
//
// char* b2i::write_file_for_illuminus (...)
//
// Write out the new Illuminus input-format file.
//
// Generate name of output file written. Return this; caller must delete[].
//
/////////////////////////////////////////////////////////////////

char* b2i::write_file_for_illuminus 
         ( std::string outfile,
           map < std::string, hash_map <std::string, std::string> >& data,
           map < std::string, int >& samples,
           Manifest& mf ) {

  
  FILE* pFile = NULL;
  char* buf;

  // Create uniquely-named temporary file:   // TODO umask check
  string tempname = "/tmp/b2i_illuminus_temp_XXXXXX";
  buf = create_tempfile (string("write_file_for_illuminus"), tempname);
	
  // Associate pFile with this file; and, if file is not already open, do so:
  pFile = fopen (buf, "w");
	

  if ( NULL == pFile) {
    cout << "\n\nError in opening output file. \n\n" << flush;
    exit(1);
  }

  // Set up sample iterator:
  map < std::string, int >::iterator sample_it;


  //
  // Header.
  //
  fprintf(pFile, "SNP\tCoor\tAlleles"); ////fout << "SNP\tCoor\tAlleles";
  for ( sample_it = samples.begin(); sample_it != samples.end(); sample_it++ ) {
    std::string a_sample = (*sample_it).first;
	char* samp;
	try {
	  samp = new char[a_sample.size() + 1];
	}
	catch (bad_alloc& ba) {
	  cerr << "\nERROR: in write_file_for_illuminus(), populating header of O/P file, caught bad_alloc: " 
           << ba.what() << "\n" << flush;
	  exit(1);
	}
    strcpy(samp, a_sample.c_str());
	fprintf (pFile, "\t%sA\t%sB", samp, samp);
    delete[] samp;
  }


  //
  // Data.
  //
  // Outer loop = probes; inner loop = samples.
  // Probe loop will proceed in the same order as Manifest's probe vector.
  //

  int numprobes = mf.snps.size();

  for ( int i = 0; i < numprobes; i++ ) {

	snpClass s = mf.snps[i];
	std::string name = s.name;

	// We don't bother outputting anything if we have no data for this probe:
	map <std::string, hash_map <std::string, std::string> >::iterator it = data.find(name);
	if ( it == data.end() ) {
	  continue; // Not found in data.
	}

	// Now get the properties for this probe.

	if ( ! s.position || ! s.snp ) {
	  cout << "\nWARNING: insufficient data found in list for SNP " << name << "." << flush;
	  continue; // Bad data - go to next in file.
	}


	// OK, we have data; let's output it. Probe details first:
	fprintf (pFile, "\n%s\t%li\t%c%c", name.c_str(), s.position, s.snp[0], s.snp[1]);

	// Now add norm X (or A) and norm Y (or B) readings for this probe for each sample:
	for ( sample_it = samples.begin(); sample_it != samples.end(); sample_it++ ) {
	  
	  std::string a_sample = (*sample_it).first;
	  std::string readings = data[name][a_sample];  // readings in form "Xnorm,Ynorm".

	  char* rdgs = (char*)readings.c_str();
	  
	  char* p = strtok (rdgs, ","); //TODO replace strtok with strstr
	  char* first;
	  char* second;
	  if ( ! p ) {
	 	first = "NaN";
		second = "NaN";
	  }
	  else {
		first = p;
		p = strtok (NULL, "\t"); //TODO replace strtok with strstr
		if ( p ) {
		  second = p;
		}
		else {
		  second = "NaN";
		}
		
	  }

	  fprintf (pFile, "\t%s\t%s", first, second);
 

	} // End of inner loop (samples)

  } // End of outer loop (probes)

  // Add newline to final line:
  fprintf (pFile, "\n");

  fclose (pFile);

  return buf;

} // End of b2i::write_file_for_illuminus()



//////////////////////////////////////////////////////////////////////
//
// bool b2i::check_file_integrity (...)
//
// Inputs: name of file to check; ref to intensity readings map;
//         ref to samples map.
//
// Checks that file:
// - has expected number of lines, given number of probes
// - has expected number of columns on each line, given no. of unique sampleIDs
//
// Returns true only if all tests passed, otherwise false.
//
//////////////////////////////////////////////////////////////////////

bool b2i::check_file_integrity ( const char* const source,
								 map < std::string, hash_map <std::string, std::string> > &readings,
								 map <std::string, int> &samples ) {

  const unsigned char delim = '\t';
  const unsigned char newline = '\n';

  unsigned int expected_num_cols = 3 + ( 2 * samples.size() );
  unsigned int expected_num_lines = 1 + readings.size();

  FILE* instream = fopen (source, "r");
  
  int saved_errno;
  
  if ( NULL == instream ) {
	saved_errno = errno;
	cout << "\nERROR in b2i::check_file_integrity(): null input stream; errno = " << saved_errno;
	//	 << errmsg << flush;
	return false;
	
  }
  
  char c;
  unsigned int numlines = 0;
  unsigned int num_words_read = 0;
  bool in_a_word = false;
  
  while ( (c = fgetc(instream)) != EOF ) {
	
	switch ( c ) {

	case delim:
	  if ( in_a_word ) {
		in_a_word = false;
		num_words_read++; // We've now just left (processed) a word.
	  }
	  break;

   
	case newline:
	  numlines++;
	  if ( in_a_word ) {
		in_a_word = false;
		num_words_read++; // as we've just left a word
	  }
	  // cout << "\nDBG: read " << num_words_read;
	  
	  if ( expected_num_cols != num_words_read ) {
		cout << "\nERROR in b2i::check_file_integrity(): line " << numlines << " of file " 
			 << source << " has " << num_words_read << " words, whereas " 
			 << expected_num_cols << " were expected.\n" << flush;
		return false;
	  }

	  num_words_read = 0;
	  break;


	default:
	  if ( isalnum(c) ) {
		in_a_word = true;
	  }
	  else {
		// Could be: colons, decimal points, underscores,
		// spaces...I'd leave this section blank if I were you!
	  }

	}
 
  }


  // EOF returned.
  // fgetc() retrurns EOF whether it's really eof or an error. Which is it?
  if ( ferror(instream) ) { 
	cout << "\nERROR (unknown) in b2i::check_file_integrity()\n";
	return false;
  }

  if ( ! feof(instream) ) {
	cout << "\nERROR (unknown) in b2i::check_file_integrity(); not ferror() and not feof() either!\n";
	return false;
  }
	
  // OK, it's eof.


  // The last line should have a newline character. Deal with situation where this is
  // not so, i.e. we have content after the last newline character in the file.

  if ( num_words_read ) {  // We have content on this line to be accounted for.
	// Last char on last line may well be an alphanumeric char

	// If we are _in_ a word when we get right to the end of the file, we need to
	// add 1 to the number of words we have processed on this line. We'll also
	// add 1 to the number of lines processed.

	if ( in_a_word ) {
	  //cout << "\nDBG: in last word of last line (no newline char)\n" << flush;
	  num_words_read++;
	  if ( expected_num_cols != num_words_read ) {
		cout << "\nERROR in b2i::check_file_integrity(): line " << numlines << " of file "
			 << source << " has " << num_words_read << " words, whereas "
			 << expected_num_cols << " were expected.\n" << flush;
		return false;
	  }
	  numlines++;
	}

  }


  // cout << "\n\nDBG: numlines: " << numlines << "\n";
  if ( expected_num_lines != numlines ) {
	cout << "\nERROR in b2i::check_file_integrity(): file " << source << " has "
		 << numlines << " lines, whereas " << expected_num_lines 
		 << " were expected\n" << flush;
	return false;
  }

  cout << "\nINTEGRITY: File " << source << " has PASSED its integrity check." << flush;
  return true;


} // End of b2i::check_file_integrity()



/////////////////////////////////////////////////////////////////
//
// main()
//
/////////////////////////////////////////////////////////////////

int main (int argc, char** argv) {

  b2i beady;

  std::string input_file, list_file, output_file, snp_file, selected_chrom, 
	output_snps_file, output_samples_file;


  // Key: probe name. value = hash_map where key = sampleID, value = string containing
  // both A and B intensities, in the form "A_intensity,B_intensity":
  map < std::string, hash_map <std::string, std::string> > intensities;

  // Use a map to keep track of every unique sampleID encountered, so
  // that they can be output in the same order as the sampleIDs in
  // the inner map above. Value is a meaningless int:
  map < std::string, int > sampleIDs;


  beady.process_command_line(argc, argv, input_file, list_file, output_file, snp_file);

  // Read in probe positions
  cout << " \nNow reading in probe positions." << flush;
  Manifest mani; 
  
  if ( beady.SELECTED_CHROM.size() > 0 ) {
    mani.open(snp_file, beady.SELECTED_CHROM);
  }
  else {
    mani.open(snp_file);  // TODO can we check for exceptions?
  }

  cout << "\nProbes were read in successfully." << flush;

  if ( beady.SORT_BY_POS ) {
	mani.order_by_position();
  }


  // Read in file(s) to be converted
  if ( input_file.size() > 0 ) {
    cout << "\nAbout to read input file..." << input_file << flush;
    beady.read_input_file (input_file, intensities, sampleIDs, mani);
  }
  else if ( list_file.size() > 0 ) {
    
    // Open the list of files, and process each one:
    char* in_string;
	try {
	  in_string = new char[list_file.size() + 1];
	}
	catch (bad_alloc& ba) {
	  cerr << "\nERROR: in main(), in_string, caught bad_alloc: " << ba.what() << "\n" << flush;
	  exit(1);
	}
    strcpy(in_string, list_file.c_str());
    //FILE* pFile = fopen (in_string, "r");
    ifstream filestr;
    filestr.open (in_string, fstream::in);
    delete[] in_string;

    if ( filestr.fail() ) {   //( NULL == pFile) {
      cout << "\n\nError in opening input list file " << list_file << ".\n\n" << flush;
      exit(1);
    }
    string line_string;
    while ( getline(filestr, line_string) ) {
      
      // Remove any whitespace characters before, in, or after filename:
      string no_ws_line_string;
      unsigned int charpos;
      for (charpos = 0; charpos < line_string.length(); charpos++)
      {
        char candidate = line_string[charpos];
        if ( ! isspace(candidate) ) {
          no_ws_line_string.append(1, candidate);
        }
      }

      // Now process the file listed in the string line_string:
      beady.read_input_file (no_ws_line_string, intensities, sampleIDs, mani);    

    } // End while.

  }
  else {
    cout << "\nNeither a single input file nor multiple files specified. Goodbye.\n";
    exit(1);
  }

  // Write out new file. 
  cout << "\nAbout to write o/p file...\n" << flush;
 
  // 'opfile' will be name of temp file created
  
  char* opfile = beady.write_file_for_illuminus (output_file, intensities, sampleIDs, mani);
  //TODO check result.
  
	
  // As we've used a temporary file, we need to copy this across to the desired
  // op file location, and check they are both the same before deleting tempfile
  beady.old_copy_file(opfile, output_file.c_str());
	
  delete [] opfile; // Deletes the string containing its name, not the file itself!

	
  // Finally we do a couple of integrity checks on the output file
  bool file_ok 
	= beady.check_file_integrity ((const char* const)output_file.c_str(), intensities, sampleIDs);
  
  if ( ! file_ok ) {
	cout << "\nINTEGRITY CHECK FAILED\n" << flush;
	cout << "\nPlease check the file manually: expected no. of lines and cols,"
		 << "\nand that no. of cols in line 1 = no. in last line.\n" << flush;
	exit (1);
  }

  
  cout << "\n\nb2i has finished.\n\n" << flush;
  

  return 0;

} // End of main().




 /*
  cout << "\nsizeof(int) is " << sizeof(int)
       << "\nsizeof(long) is " << sizeof(long)
       << "\nsizeof (long int) is " << sizeof (long int)
       << "\nsize of float is " << sizeof(float)
       << "\nsize of char is " << sizeof(char)
       << "\nsizeof snpClass is " << sizeof (snpClass)
       << "\nsizeof string is " << sizeof(string)
       << "\nsizeof Manifest is " << sizeof(Manifest)
       << "\nsize of vector is " << sizeof(vector<snpClass> )
       << "\nsizeof hash_map<std::string, std::string> is "<< sizeof(hash_map<std::string, std::string>)
       << "\n size of map<std::string, hash_map <std::string, std::string> > is "
       << sizeof ( map < std::string, hash_map <std::string, std::string> > )
       << "\nsizeof map < std::string, int > is " << sizeof( map < std::string, int >)
       << "\n\n";

  exit(1);
  */


// rm 9.txt; ./b2i -s /nfs/new_illumina_geno04/call/Human670-QuadCustom_v1_A.bpm.csv -i input.test -o 9.txt -c 9



/*

./b2i -s /nfs/wtccc/data5/genotyping/chip_descriptions/Illumina_Infinium/Human_Cardio-Metabo/Cardio-Metabo_Chip_11395247_A.bpm.csv -i av.test -o av.out -c 1

*/

// EOF
