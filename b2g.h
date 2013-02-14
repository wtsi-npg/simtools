// b2g.h
//
// Matthew Gillman, WTSI, 13th October 2010.
//
// $Id: b2g.h 1331 2010-10-25 16:07:24Z mg10 $
//
// The b2g class is to convert BeadStudio Final Call Report files
// into files required for GenoSNP input.
//

#ifndef _B2G_H
#define _B2G_H

#include "b2base.h"

class b2g : public b2base {

 public:

  b2g () : b2base() {

	XRaw = "X Raw";

	YRaw = "Y Raw";
  
  }

	
  ~b2g () {}



  void process_command_line (int argc, char** argv,
							 string& manifest_file,
							 string& input,
							 string& listing,
							 string& snps,
							 string& data,
							 bool& wf);


  void usage();


  void initialise_colpos ();

  
  bool check_file_integrity 
	( const char* const filename,
	  map < std::string, hash_map <std::string, std::string> > &readings,
	  map <std::string, int> &samples );


  void process_input_file ( std::string input, 
							FILE* outFile );


  void write_SNPs_file ( std::string& final );
  

  FILE* create_samples_file ( string& tempname );


  void write_line_to_samples_file ( FILE* outFile,
                                    map<string, string>& mymap,
									string& sampleID );

  string XRaw, YRaw;

  Manifest mf;

};


#endif // _B2G_H
