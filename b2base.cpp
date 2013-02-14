// b2base.cpp
//
// Matthew Gillman, WTSI, 13th October 2010.
//
// $Id: b2base.cpp 1380 2010-11-25 15:43:17Z mg10 $
//
// b2base class implementation.
//

#include "b2base.h"



////////////////////////////////////////////////////////////////////////
//
// bool b2base::check_colpos_complete()
//
// Checks that colpos has values for all required fields.
//
////////////////////////////////////////////////////////////////////////

bool b2base::check_colpos_complete() {

  bool ok = true;
  map <string, int>::iterator it;

  for ( it = colpos.begin(); it != colpos.end(); it++ ) {
	string field = (*it).first;
	int value = (*it).second;
	//cout << "\nDBG: check_colpos_complete(): colpos[" << field << "] = " 
	//	 << value << flush;
	if ( -1 == value ) {
	  ok = false;
	}
  }
  
  return ok;

}



////////////////////////////////////////////////////////////////////////
//
// void b2base::extract_datum_from_string()
//
// Extracts item of data lying between two delimiters in a string.
// Places this in 'datum'.
//
// INPUTS:
//
// unsigned int* delimlocs : ptr to array of delimiter locations in linestring
// unsigned int num_delimlocs : the number of entries in the delimlocs array.
// const char* linestring  : the string containing the data
// unsigned int length     : the length of linestring
// std::string& datum      : the string to be initialized
// unsigned int col        : the column number of this item of data
// 
///////////////////////////////////////////////////////////////////////

void b2base::extract_datum_from_string 
                  (unsigned int* delimlocs,
				   unsigned int num_delimlocs,
                   const char* linestring,
				   unsigned int length, 
				   std::string& datum,
				   unsigned int col) {

  // Delimiter (D) positions: zeroth position is at delimlocs[0], 1st at delimlocs[1], etc.

  // If first char in linestring is a delimiter (unlikely), the nth column will lie between 
  // (and excluding) the nth and n+1th delimiter positions in the string:
  // D---D-------D-----D--etc
  // 0   1       2     3     delimiters
  //   0     1      2     3  columns  
 
  // But generally this will be true:
  // --------D----D---------D-- etc
  //         0    1         2    delimiters
  //    0       1      2      3  columns
  // in which case the 0th column lies from linestring[0] and upto but excluding
  // the 0th delimiter position; the 1st between (but excluding) the 0th and 1st positions; etc.

  // NB the last column in the string, number LAST_COL_NUM, will be just after the (final) delimiter:
  // D--- at the end of linestring. (If the final item in the string is itself a delimiter, the final
  // column will still exist but have no content(zero length)). 


  // The highest delim number (in positions array) is one less 
  // than the number of delimiters in the array.
  unsigned int final_delim = num_delimlocs - 1;

  //cout << "\nDBG in extract_datum_from_string(); num_delimlocs = "
  //<< num_delimlocs << "; length of linestring = " << length << endl << flush;


  // DEBUG
  /*
     cout << "\nDEBUG; linestring = *" << linestring << "*\n" << flush;
     for ( unsigned int i = 0; i < num_delimlocs; i++ ) {
	   unsigned int apos = *(delimlocs + i);
	   cout << "\nDEBUG: when i = " << i << ", position is " << apos << endl << flush;
     }
  */  


  // Deal with first (unlikely) situation above:
  if ( 0 == *delimlocs ) {
	unsigned int startindex =  *(delimlocs + col) + 1; // Start char is 1 beyond delimloc number 'col'.
     
     // linestring will stretch from [0] to [x], so its length is (x+1). The biggest possible
 	 // valid value for startindex is x. If we have gone beyond this we are in danger!
     if ( startindex >= length ) {
	   cout << "\nERROR in extract_datum_from_string(). Data are:"
			<< "\nlinestring: *" << linestring << "*"
			<< "\nlength: *" << length << "*"
			<< "\ncol: *" << col << "*\n" << flush;
	   exit(1);
     }
     
     if ( LAST_COL_NUM == col ) {  // 3 in top figure above.
	   // cout << "\nDEBUG: startindex = " << startindex 
	   //      << ", length = " << length << endl << flush;
       datum = string(linestring + startindex, length - startindex); // e.g. linestring is [0] to [9]
     }                                                               // startindex = 6, length = 10,
                                                                     // len-start = 4 (6,7,8,9). OK. 
     else {
       unsigned int len = *(delimlocs + col + 1) - startindex; 
       datum = string(linestring + startindex, len); 
     }
     
  }
  else { // Normal (second) situation:
    
    if ( col == 0 ) { // Deal with situation where datum occurs before 0th delimiter:
      unsigned int len = *delimlocs;   // Position of the zeroth delim = 1 beyond last char
      datum = string(linestring, len); // String runs from start (0) of linestring for len chars
                                       // so last one is 1 before zeroth delim
    }
    else {
      unsigned int low = *(delimlocs + col - 1) + 1;  // Pos. of 1st char 
                                                      // needed. (e.g. 7)
	  unsigned int high; // 1 beyond last char needed (e.g. 21).

	  // If we have passed the final delimiter, we take account of this fact:
	  
	  if ( col > final_delim ) { 
		high = length;  
		// cout << "\nDEBUG: gone past : col = " << col 
		// << ", final_delim = " << final_delim << endl << flush;
	  }
	  else {
		high = *(delimlocs + col); 
		// cout << "\nDEBUG: still in : col = " << col 
	    // << ", final_delim = " << final_delim << endl << flush;
	  }
	  // cout << "\nDEBUG: low = " << low << ", high = " << high << endl << flush;

      datum = string(linestring + low, high - low); // length = 21-7 = 14
                                                    // ( = 7,8,...19,20)
	  //cout << "\nDEBUG: low = " << low << ", high = " << high 
	  //      << ", datum =*" << datum << "*\n" << flush;
	  
    }

  }


  //cout << "\nDEBUG: col = " << col << ", datum = '" << datum << "' " << endl << flush;

}

 

//////////////////////////////////////////////////////////////////////
//
// const char* get_file_hash (const char* const filename,
//                            const char* digest = "sha1");)
//
// Get hash (default: SHA1) of file 'filename'.
//
// You can use other digests, e.g. "md5"
//
// Typical SHA1 file hash: e8f754b02379248931739e550fa50bb54679ea71
//
// Caller must delete[] char array returned.
//
// http://www.cplusplus.com/forum/unices/3929/
//
//////////////////////////////////////////////////////////////////////

 
const char* b2base::get_file_hash (const char* const filename,
								   const char* digest) {

  ifstream filestr;
  filestr.open ( filename, fstream::in );
  if ( filestr.fail() ) {
    cout << "\n\nError in  b2base::get_file_hash(): opening input file\n\n" 
		 << flush;
    exit(1);
  }

  EVP_MD_CTX mdctx;
  const EVP_MD *md;
  unsigned char md_value[EVP_MAX_MD_SIZE];
  unsigned int md_len, i;
  
  OpenSSL_add_all_digests();
  
  md = EVP_get_digestbyname(digest);
  
  if ( ! md ) {
	cout << "\nERROR in b2base::get_file_hash():"
		 << " Unknown message digest " << digest << "\n" << flush;
	exit(1);
  }
  
  EVP_MD_CTX_init(&mdctx);
  if ( ! EVP_DigestInit_ex(&mdctx, md, NULL) ) {
	cout << "\nERROR in b2base::get_file_hash(): EVP_DigestInit_ex()\n"
		 << flush;
	exit(1);
  }
  
  string str;
  int numread = 0;
  char c;
  char* buffer = new char[ EVP_MAX_MD_SIZE ];

  // Loop through the contents of the file.
  // DO NOT use getline() as it strips off the '\n' at
  // the end of the line, but we need it for computing the file's hash.

  while ( filestr.get(c) ) {

	buffer[numread] = c;
	numread++;
	
	if ( EVP_MAX_MD_SIZE == numread ) {
	  
	  if ( ! EVP_DigestUpdate(&mdctx, buffer, EVP_MAX_MD_SIZE) ) {
		cout << "\nERROR in b2base::get_file_hash(): EVP_DigestUpdate()\n"
			 << flush;
		exit(1);
	  }
	  
	  numread = 0;
	  delete[] buffer;
	  buffer = new char[ EVP_MAX_MD_SIZE ];		
  
	}

  }
	
  // Deal with any remainder
  if ( numread > 0 ) {
	
	if ( ! EVP_DigestUpdate(&mdctx, buffer, numread) ) {
	  cout << "\nERROR in b2base::get_file_hash() (remainder): "
		   << "EVP_DigestUpdate()\n"
		   << flush;
	  exit(1);
	}
	delete[] buffer;
	
  }
  

  if ( ! EVP_DigestFinal_ex(&mdctx, md_value, &md_len) ) {
	cout << "\nERROR in b2base::get_file_hash(): EVP_DigestFinal_ex()\n"
		 << flush;
	exit(1);	  
  }


  EVP_MD_CTX_cleanup(&mdctx);
	
  //printf("\nDigest is: ");
  // for(i = 0; i < md_len; i++) printf("%02x", md_value[i]);
  //printf("\n");
  
  char* result = new char[2*md_len + 1];

  for ( i = 0; i < md_len; i++ ) {
	
	char temp[3]; // 2 chars + terminating null.
	sprintf(temp, "%02x", md_value[i]);
	// cout << " *" << temp << "*";
	result[2*i] = temp[0];
	result[2*i + 1] = temp[1];

  }

  result[2 * md_len] = '\0';

  cout << "\nFile " << filename << " has hash value " << result << "." << flush;

  return result; // Caller must delete[]

}



/////////////////////////////////////////////////////////////////
//
// unsigned int b2base::check_substring ( const char* const content, int col )
//
// Given string 'content' (e.g. the entire header line, or part of it),
// see if the start of it matches any of the required column names.
// If it does, assign the value of the relevant column to 'col'.
//
// Returns length of substring matched.
//
/////////////////////////////////////////////////////////////////

unsigned int b2base::check_substring ( const char* const content, int col ) {

  // cout << "\nDBG: check_substring(), content = " << content << endl << flush;

  
  // Does the content match any of our required fields?
  map <string, int>::iterator it;
  it = colpos.find(content);

  if ( colpos.end() == it ) {
	cout << "\n(Not required) " <<  content << " col. no. is " << col << "." << flush;
	return 0;
  }
  else {
	colpos[content] = col;
	cout << "\n" << content << " col. no. is " << col << "." << flush;
	return strlen(content);
  }


}



////////////////////////////////////////////////////////////////////////
//
// b2base::find_column_numbers() - find column numbers in the BeadStudio data file.
//
// INPUTS:
// std::string& mystring  : header line immediately above values in the file.
// std::string& delimiter : delimiter to use when breaking up line.
// 
// WARNING: at present this cannot cope with a space, or spaces, as delimiter.
//
// RETURNS:
// bool value - true if success, false otherwise
//
// This is a helper function for read_input_file().
//
// TODO add ellipsis so you can choose which columns it searches for...
// TODO if it doesn't find all the required column numbers, report an error
//
///////////////////////////////////////////////////////////////////////

bool b2base::find_column_numbers (std::string &mystring, std::string &delimiter) {

  // cout << "\nDEBUG: in find_column_numbers(), mystring = *" << mystring << "*\n" << flush;


  char* header;
  try {
	header = new char[mystring.size() + 1];
  }
  catch (bad_alloc& ba) {
    cerr << "\nERROR: in find_column_numbers(), header, caught bad_alloc: " 
         << ba.what() << "\n" << flush;
	exit(1);
  }
  strcpy(header, mystring.c_str());

  char* the_delim;
  try {
	the_delim = new char[delimiter.size() + 1];
  }
  catch (bad_alloc& ba) {
    cerr << "\nERROR: in find_column_numbers(), the_delim, caught bad_alloc: " 
         << ba.what() << "\n" << flush;
	exit(1);
  }
  strcpy(the_delim, delimiter.c_str());
 
  char* p = strtok(header, the_delim); //TODO replace strtok with strstr ??

  int col = 0; // This must start at 0 to be consistent with vector indexing.
 
  while ( p != NULL) {
	
	if ( p != NULL ) {
      //cout << "\nDBG: find_column_numbers(), col = " << col << ", p = *" << p << "*\n" << flush;
    }
	//	else { cout << "\nNULL\n" << flush; }

    LAST_COL_NUM = col; // Final value from this loop will be the proper value for this.
	//cout << "\nDEBUG: find_column_numbers(), LAST_COL_NUM = " << LAST_COL_NUM << endl << flush;

	check_substring(p, col);

    p = strtok(NULL, the_delim); //TODO replace strtok with strstr
	//if ( p != NULL ) {cout << "\nDBG: find_column_numbers(), col = " << col << ", p = *" << p << "*\n" << flush;}

    col++;
  }

  delete[] header;
  delete[] the_delim;

  bool try_again = false;

  if ( ! check_colpos_complete() ) {
	try_again = true;
  }

  if ( try_again ) {
    cout << "\n\nOops: columns not found in input file header...wrong format...try again...\n\n" << flush;
    return false;
  }

  return true;

} // End of b2base::find_column_numbers()



/////////////////////////////////////////////////////////////////
//
// char* b2base::create_tempfile ( std::string func, std::string basis )
//
// Creates a tempfile using the string 'basis', called
// from function 'func'. 
//
// Closes the file, then returns its name. Caller must delete[] it.
//
// We use mkstemp() to avoid potential race conditions with other
// processes.
//
// The last 6 characters of 'basis' must be "XXXXXX". 'basis' should
// include the path to the file (unless desired to create it in 
// current directory). Best for basis to be "/tmp/something_XXXXXX".
//
/////////////////////////////////////////////////////////////////

char* b2base::create_tempfile ( std::string func, std::string basis ) {

  char* tn = (char*) basis.c_str();
  char* buf;

  try {
	buf = new char[strlen(tn) + 1];
  }
  catch (bad_alloc& ba) {
	cerr << "\nERROR: in " << func << "(), buf alloc, caught bad_alloc: " 
		 << ba.what() << "\n" << flush;
	exit(1);
  }

  strcpy(buf, tn);

  int fd = mkstemp(buf);
  int stored_errno = errno;
  
  if ( fd == -1 ) {
	cout << "\nError in create_tempfile(), called from " << func << ": " << flush;

	switch ( stored_errno ) {

	case EEXIST:
	  cout  << "\nCould not create a unique temporary filename.\n" << flush;
	  break;
	case EINVAL:
	  cout << "\nThe last 6 characters of 'basis' were not XXXXXX.\n" << flush;
	  break;
	default:
	  cout << "\nUnknown error, with errno = " << stored_errno << ".\n" << flush;
	}

	exit(1);
  }

  // Now can use fopen to associate FILE* stream with the filename,
  // even though it is already open. (Do this in client code.)

  return buf;

} // End of b2base::create_tempfile().



/////////////////////////////////////////////////////////////////
//
// bool b2base::copy_file (const char* const source, const char* const dest)
//
// Copies file "source" to new file "dest".
//
/////////////////////////////////////////////////////////////////

bool b2base::copy_file (const char* const source, const char* const dest) {

  char* errmsg = "\nFile copy aborted (use the temp file instead).\n";

  int saved_errno = 0;

  size_t num_read, num_written;
  size_t size = 1;
  size_t nmemb = 1024;
  
  char buffer [size * nmemb];

  FILE* instream = fopen (source, "r");
  if ( NULL == instream ) {
	saved_errno = errno;
	cout << "\nERROR in b2base::copy_file: null input stream; errno = " << saved_errno
		 << errmsg << flush;
	return false;

  }

  FILE* outstream = fopen (dest, "w+");
  if ( NULL == outstream ) {
	saved_errno = errno;
	cout << "\nERROR in b2base::copy_file: null output stream; errno = " << saved_errno
		 << errmsg << flush;
	return false;
  }
  
  // Get file descriptor for the output stream (neneded with fsycn to force flushing;
  // we need to do this to ensure the entire new file has been written before we
  // find its checksum, otherwise this will be wrong.
  int fd = fileno(outstream);
  if ( -1 == fd ) {
	saved_errno = errno;
	cout << "\nERROR in b2base::copy_file: unable to get file descriptor for output stream; errno = " 
		 << saved_errno << errmsg << flush;
	return false;
  }

  clearerr(instream);
  clearerr(outstream);

  while ( (num_read = fread (buffer, size, nmemb, instream)) ) {

	//cout << "\nDBG: numread: " << num_read;

	if (num_read != size * nmemb) break;

	num_written = fwrite (buffer, size, nmemb, outstream);

	if ( 0 != fflush(outstream) ) {
	  saved_errno = errno;
	  cout << "\nERROR in b2base::copy_file: fwrite failed; errno = " 
		   << saved_errno << errmsg << flush;
	  return false;
	}

	if ( -1 == fsync(fd) ) {
	  saved_errno = errno;
	  cout << "\nERROR in b2base::copy_file: fsync failed; errno = " 
		   << saved_errno << errmsg << flush;
	  return false;
	}

	//cout << "\nDBG: numwritten: " << num_written;

	if ( num_written != num_read ) { // Error!
	  saved_errno = errno;
	  cout << "\nERROR in  b2base::copy_file(): did not write all bytes;"
		   << " errno = " << saved_errno << "\n" << flush;
	  return false;
	}

  }  

  // Jumping to here means we have either finished reading the file (EOF),
  // or we've reached an error condition. NB There may (or may not) be a
  // few bytes left to write (but not if file size is an exact multiple
  // of size * nmemb).

  if ( feof(instream) ) {
	if ( num_read ) {
	  num_written = fwrite (buffer, 1, num_read, outstream);

	  if ( 0 != fflush(outstream) ) {
		saved_errno = errno;
		cout << "\nERROR in b2base::copy_file: fwrite failed (last block); errno = " 
			 << saved_errno << errmsg << flush;
		return false;
	  }

	  if ( -1 == fsync(fd) ) {
		saved_errno = errno;
		cout << "\nERROR in b2base::copy_file: fsync failed (last block); errno = " 
			 << saved_errno << errmsg << flush;
		return false;
	  }

	  //cout << "\nDBG:numwritten: " << num_written;
	  if ( num_written != num_read ) {
		saved_errno = errno;
		cout << "\n(Jumped) ERROR in b2base::copy_file(): "
			 << "did not write all bytes; errno = " 
			 << saved_errno << "\n" << flush;
		return false;
	  }
	}

  }
  else {
	if ( ferror(instream) ) {
	  saved_errno = errno;
	  cout << "\nERROR in b2base::copy_file() in reading file; errno = " 
		   << saved_errno << "\n" << flush;
	  return false;
	}
  }


  // OK so far. Just need to check the file hashes are the same.

  const char* firstHash = get_file_hash((const char* const)source);
  const char* secondHash = get_file_hash((const char* const)dest);
  
  bool goodbye = false;
  if ( strcmp (firstHash, secondHash) != 0 ) {
	goodbye= true;
  }
  delete[] firstHash;
  delete[] secondHash;
  
  if ( goodbye ) {
	cout << "\nDANGER the temp file and final file have different hashes,"
		 << "\nsuggesting the copy failed. Please check the temp file created:"
		 << "\ndoes it have the expected number of lines and columns, and"
		 << "\ndoes it have the same number of columns in first and last lines?"
		 << "\nIf it seems OK then copy it across manually.\n" << flush;
	exit(1);
  }
  else {
	cout << "\nPhew! The temp file was copied OK, so I'm going to delete it." 
		 << flush;
	
	// Temp file should be 'opfile' at this stage, but just in case something
	// has gone wrong we'll double check this is so before deletion.
	if ( strcmp (source, dest) == 0 ) {
	  cout << "\nHmm, very strange: the temp file and final file appear to be"
		   << "\n  the *same* file (i.e. not just same content!). I won't delete it."
		   << flush;
	}
	else {
	  unlink(source);
	}
  }



  return true;

}


//
// Old version, which used execute_process().
//
 
bool b2base::old_copy_file (const char* const source, const char* const dest) {

  return copy_file (source, dest);

}



//////////////////////////////////////////////////////////////////////
//
// const char* execute_process (const char* const cmd)
//
// Opens a process and execute the command given in 'cmd'.
// Returns the result of this comand.
//
// Caller must check to see if result is NULL. For commands with no
// output, of course, this is expected; otherwise, this indicates an
// error.
//
// Also, caller must delete[] the result.
//
//////////////////////////////////////////////////////////////////////
/*
const char* b2base::execute_process (const char* cmd) {
  
cout << "\nDBG:Attempting to execute cmd: " << cmd << flush;
  
int bufSize = 2000; // Should be plenty big enough.
char* data;
try {
data = new char[bufSize];
}
 catch (bad_alloc& ba) {
 cerr << "\nERROR: in execute_process(), caught bad_alloc: " << ba.what() << "\n" << flush;
exit(1);
}
fflush(NULL); // flush all open output streams

FILE* processStream = popen(cmd, "r");

int tmp = errno;
if ( tmp != 0 ) {
cout << "\n value of errno " << tmp << "\n" << flush;
}

int saved_errno = 0;

if ( processStream == NULL ) { // Error
saved_errno = errno;
cout << "\nERROR: NULL popen result in execute_process()" << flush;
if (saved_errno) {
cout << "\nValue of errno = " << saved_errno << endl << flush;
}
else {
cout << "\nNo errno set, so probably a memory allocation problem\n" 
<< flush;
}
exit (1);

}
  
char* result = fgets(data, bufSize, processStream);
if ( ferror(processStream) ) {
cerr << "\nERROR: in execute_process(), got ferror on stream\n" << flush;
}


if ( result == NULL ) { 
pclose(processStream);
delete[] data;
return NULL;
}


// Replace any space or newline character in the string containing
// the output data from the command with a terminating null.
int length = strlen(data);
for (int i = 0; i < length; i++) {
if ( (data[i] == ' ') || (data[i] == '\n') ) {
data[i] = '\0';
}
}

int pclose_ret = pclose(processStream);
  
if ( -1 == pclose_ret ) {
saved_errno = errno;
cout << "\nERROR:Got error in execute_process(), pclose() step..." 
<< "\nerrno = " << saved_errno << "\n" << flush;
}

// cout << "\nDBG: Result from cmd is: " << data << ".\n";
return (const char*) data; 

} // End of b2base::execute_process()

*/
