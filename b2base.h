// b2base.h
//
// Matthew Gillman, WTSI, 13th October 2010.
//
// $Id: b2base.h 1380 2010-11-25 15:43:17Z mg10 $
//
// The b2base class is the base class for the b2i and b2g classes. It attempts
// to host the functions and variables which they share.
//


#ifndef _B2_BASE_H
#define _B2_BASE_H


#include "Manifest.h"

#include <iostream>
#include <fstream> 
#include <unistd.h>
#include <cstdio>
#include <ctime>
#include <locale>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <openssl/evp.h>

using namespace std;


/* 
 
Input will be a BeadStudio file of format similar to:
 
 [Header]
 BSGT Version    3.3.4
 Processing Date 2/22/2009 1:58 PM
 Content         Human670-QuadCustom_v1_A.bpm
 Num SNPs        660447
 Total SNPs      660447
 Num Samples     683
 Total Samples   683
 [Data]
 SNP Name  Sample ID  Allele1 - Top  Allele2 - Top  GC Score  Theta  R  X  Y  X Raw  Y Raw  ...
 200003  88036_A01_WTCCCT492730  A   A  0.9299  0.034   1.046   0.994   0.052   9788    997 
 200006  88036_A01_WTCCCT492730  G   G  0.7877  0.982   1.487   0.040   1.447   848     16488 
 
 
 NB (1) It may not have all the columns above. In fact this program should work if the
        only columns in the [Data] section are "SNP Name" and "Sample ID",
        plus "X" and "Y" for Illuminus or "X Raw" and "Y Raw" for GenoSNP.. 
    (2) It may be comma-delimited instead of tab-delimited.


 For Illuminus, output of this prog must be a single file of the format...
 
SNP         Coor  Alleles 74230_A08_sample1A 74230_A08_sample1B 74230_B08_sample2A  ...   mmmmm_nnn_sampleNB
rs2461547   60864   TC               0.151              1.858           0.140        ...        1.483
cnvi0070164 62810   NA               1.635              0.063           1.297        ...        0.081
 
NB A == X and B == Y.


NOTE: column numbers in this program are numbered 0 onwards. The presence of two adjacent delimiters
implies there is an empty column between them, and this still has a valid column number.
           
*/


class b2base {

 public:

  b2base () {
  
	//reset_columns();

	SORT_BY_POS = true;

	SNP_name = "SNP Name";

	SampleID = "Sample ID";

  }

  virtual ~b2base() {}

  virtual void usage() = 0;
  
  void extract_datum_from_string (unsigned int* delimlocs,
								  unsigned int num_delimlocs, 
								  const char* linestring,
								  unsigned int length,
								  std::string& datum,
								  unsigned int col);
  
   const char* get_file_hash (const char* const filename, const char* digest = "sha1");
  
  unsigned int check_substring ( const char* const p, int col );
  
  char* create_tempfile ( std::string func, std::string basis );
  
  bool copy_file (const char* const source, const char* const dest);

  bool old_copy_file (const char* const source, const char* const dest);

  bool SORT_BY_POS;


 protected:

  bool check_colpos_complete();

  // }{const char* execute_process (const char* cmd);

  /*
	virtual void read_input_file (std::string input, 
	map < std::string, hash_map <std::string, std::string> >& readings,
	map <std::string, int>& samples,
	Manifest& mf) = 0;
  */

  bool find_column_numbers (std::string &mystring, std::string &delimiter);

  // Key = name of a column (e.g. "X Raw"); value = column number (0 upwards).
  // NB important - it should _only_ store the particular fields required for the child class.
  map <string, int> colpos; 

  virtual void initialise_colpos () = 0;

  unsigned int LAST_COL_NUM; // Columns are numbered 0,1,2,...,LAST_COL_NUM

  string SNP_name, SampleID;

};



#endif // _B2_BASE_H
