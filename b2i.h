// b2i.h
//
// Matthew Gillman, WTSI, 13th October 2010.
//
// $Id: b2i.h 1380 2010-11-25 15:43:17Z mg10 $
//
// The b2i class exists to convert BeadStudio Final Call Report files
// into files required for Illuminus input. This format is also
// used with CNV pipelines.
//

#ifndef _B2I_H
#define _B2I_H

#include "b2base.h"

class b2i : public b2base {

 public:

  b2i () : b2base() {

	Xnorm = "X";

	Ynorm = "Y";	 

  }
	
  ~b2i () {}

  void process_command_line (int argc, char** argv,
						   std::string& input,
						   std::string& listing,
						   std::string& output,
						   std::string& snp_file);

  void usage();


  char* write_file_for_illuminus 
         ( std::string outfile,
           map < std::string, hash_map <std::string, std::string> >& data,
           map < std::string, int >& samples,
           Manifest& mf );


  bool check_file_integrity
	( const char* const source,
	  map < std::string, hash_map <std::string, std::string> > &readings,
	  map <std::string, int> &samples );


  void read_input_file (std::string input, 
						map < std::string, hash_map <std::string, std::string> >& readings,
						map <std::string, int>& samples,
						Manifest& mf);

  
  string SELECTED_CHROM; // Used if we are only processing one chrom (usual case).

  string Xnorm, Ynorm;

  void initialise_colpos ();

};


#endif // _B2I_H
