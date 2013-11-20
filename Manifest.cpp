//
// Manifest.cpp
//
// Author: Jennifer Liddle (js10)
//
// $Id: Manifest.cpp 1354 2010-11-11 16:20:09Z js10 $
//
#include "Manifest.h"
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

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstdio>


using namespace std;


///////////////////////////
//
// Constructor
//
//////////////////////////

Manifest::Manifest(void) 
{
  filename = "";

  EXCLUDE_CNVS = false;

#ifdef _DEBUG      // compile with g++ -D_DEBUG Manifest.cpp ...
  test_convert();
#endif

}



//////////////////////////////////////////
//
// void Manifest::open (...)
//
// Overloaded function. This version will
// only process SNPs from the specified 
// chromosome, rather than every chromosome.
//
// INPUTS:
// string filename: SNP file to be parsed.
// string chromosome: only SNPs in this chromosome
//   will be added. e.g. "1", "23", "X"; format
//   must match that used in the SNP file.
//
/////////////////////////////////////////

void Manifest::open (string filename, string chromosome, bool wide) {

  cout << "\nDEBUG: in chrom-specific open()..\n\n";

  this->selectedChromosome = chromosome;
  open (filename, wide);

}



//////////////////////////////////////////
//
// void Manifest::open (string filename, bool wide = false;) 
//
// Open file and parse it one line at a time.
//
/////////////////////////////////////////

void Manifest::open(string filename, bool wide) 
{
	string s;
	ifstream file;
	snpClass *snp = new snpClass();

	map<string, int> widecols; // Only used if opening a wide-format file.
                               // Key = col. name; value = col. number (0 onwards).

	file.open(filename.c_str());
	if (!file) {
	  cout << "Can't open file: " << filename << endl << flush;
		exit(1);
	}

	// Deal with header line(s).
	if ( wide ) {
	  bool gotcols = false;
	  string findme = "IlmnID";
	  while ( getline(file,s) ) {
		if ( s.find(findme) != string::npos ) {
		  cout << "\nDBG: got column header line.";
		  gotcols = find_wide_columns(s, widecols);
		  break;
		}
	  }
	  if ( ! gotcols ) {
		cout << "\nError: couldn't find a line containing " 
			 << findme << "\n" << flush;
		exit(1);
	  }
	}
	else {
	  getline(file,s);
	}


	// Set up column numbers. Numbers will be -1 if not found/not available.
	int INDEX_COL, NAME_COL, CHROM_COL, POS_COL, SCORE_COL, SNP_COL, 
	  I_COL, C_COL, NORMID_COL, BEADSETID_COL;

	if ( wide ) {

	  INDEX_COL     = get_map_value(widecols, "Index");          // Probably -1
	  NAME_COL      = get_map_value(widecols, "Name");
	  CHROM_COL     = get_map_value(widecols, "Chr");
	  POS_COL       = get_map_value(widecols, "MapInfo");
	  SCORE_COL     = get_map_value(widecols, "GenTrain Score"); // Probably -1
      SNP_COL       = get_map_value(widecols, "SNP");
	  I_COL         = get_map_value(widecols, "IlmnStrand");
	  C_COL         = get_map_value(widecols, "SourceStrand");
	  NORMID_COL    = get_map_value(widecols, "NormID");         // Probably -1 
	  BEADSETID_COL = get_map_value(widecols, "BeadSetID");

	}
	else {
	  // Index,Name,Chromosome,Position,GenTrain Score,SNP,ILMN Strand,Customer Strand,NormID
	  //   0     1       2        3             4       5      6               7          8
      if (s.find("Index") != 0) {
		cerr << "Manifest::open(" << filename << ")  :  invalid file header (" << s.find("Index") << ")" << endl;
	    exit(1);
	  }

	  INDEX_COL = 0;
	  NAME_COL = 1;
	  CHROM_COL = 2;
	  POS_COL = 3;
	  SCORE_COL = 4;
	  SNP_COL = 5;
	  I_COL = 6;
	  C_COL = 7;
	  NORMID_COL = 8; 
	  BEADSETID_COL = -1; 
	}

	//
	// OK, now ready to acquire data
	// If it's wide. make sure we stop at end of data section.
	int numsnps = 1; // Use if INDEX_COL missing; only increment on success.
	while ( getline(file,s) ) {	// read line by line

	  //cout << "\nDEBUG: " << s << flush;


	  if ( wide ) {
		size_t found = s.find("[Controls]");
		if ( found != string::npos ) {
		  break;
		}
	  }

	  // Store each token (substring) from s, inclusing some which may be
	  // empty, in same order as they are encountered in s.
	  vector<string> a; 

	  string DELIM = ",";
	  size_t start = 0;
	  size_t found = s.find(DELIM, start);

	  while ( found != string::npos ) {

		string datum = s.substr(start, found-start);
		a.push_back(datum);
		start = found + DELIM.length();
		found = s.find(DELIM, start);
	  }

	  // Get last (or only!) candidate:
	  string dat = s.substr(start, s.length() - start);
	  a.push_back(dat);

	  if ( (a.size() != 9) && (! wide) ) { 
		string err="line too short or too long"; 
		err += "\n";
		err += s;
		throw err; 
	  }

	  // DBG
	  /*
		cout << "\nsize = " << a.size() << endl << flush;
		for ( int i = 0; i < a.size(); i++ ) {
		string thing = a[i];
		cout << "\n" << i << ": ";
		if ( thing.empty() ) { cout << " NULL" << flush; }
		else { cout << thing << flush; }
		}
		cout << "\nNow size = " << a.size() << endl << flush;
		exit(1);
	  */

		
		if ( -1 == INDEX_COL ) {
		  snp->index = numsnps;
		}
		else {
		  snp->index = atoi(a[INDEX_COL].c_str());
		}

		// We always need the probe name.
		if ( -1 == NAME_COL ) {
		  cout <<"\nSadly, the name column was not found in the manifest file. Goodbye.\n" << flush;
		  exit(1);
		}
		snp->name = a[NAME_COL];
		if ( EXCLUDE_CNVS && ( 0 == strncmp((snp->name).c_str(), "cnv", 3) ) ) {
		  continue;
		}

		if ( -1 == CHROM_COL ) {
		  snp->chromosome = "??";
		}
		else {
		  snp->chromosome = a[CHROM_COL];
		}

		// Always store mitochondrials as "MT"!
		if ( ( 0 == snp->chromosome.compare ("M") )
			 ||
			 ( 0 == snp->chromosome.compare ("Mt") ) )
		{
		  snp->chromosome = "MT";
		}
		
		// Filter on selected chromosome if necessary:
		if ( this->selectedChromosome.length() > 0 ) {
		  if ( 0 != this->selectedChromosome.compare(snp->chromosome) ) {
		    continue;
		  }
		  else {
		    //    cout << "\nDEBUG: adding in SNP " << snp->name;
		  }
		}

		if ( -1 == POS_COL ) {
		  cout << "\nPosition column not found in Manifest file. Cheerio.\n" << flush;
		  exit(1);
		}
		else {
		  snp->position = atol(a[POS_COL].c_str()); 
		  // cout << "\nDBG; pos is " << a[POS_COL].c_str() << flush;
		}

		if ( -1 == SCORE_COL ) {
		  snp->score = -1;
		}
		else {
		  snp->score = atof(a[SCORE_COL].c_str());
		}

		if ( -1 == I_COL ) {
		  snp->iStrand = '?';
		}
		else {
		  snp->iStrand = a[I_COL].at(0);
		}
		
		if ( -1 == C_COL ) {
		  snp->cStrand = '?';
		}
		else {
		  snp->cStrand = a[C_COL].at(0);
		}

		if ( -1 == NORMID_COL ) {
		  snp->normId = -1;
		}
		else {
		  snp->normId = atoi(a[NORMID_COL].c_str());
		}

		if ( -1 == SNP_COL ) {
		  cout << "\nOops! SNP column not found in Manifest file!\n" << flush;
		  exit(1);
		}
		else {
		  convert (snp, a[SNP_COL]);
		}

		if ( -1 == BEADSETID_COL ) {
		  snp->BeadSetID = -1;
		}
		else {
		  snp->BeadSetID = atoi(a[BEADSETID_COL].c_str());

		}
		
		snps.push_back(*snp);

		/*
		  cout << "\nDBG: added in probe " << snp->name 
		  << "; pos = " << snp->position << ", snp: " 
		  << snp->snp[0] << snp->snp[1] << flush;
		*/

		normIdMap[snp->normId] = 1;

		numsnps++; // If we got this far, must be OK to increment!
	}

	file.close();

	map<int,int>::iterator i;
	int n;
      	for (n=0, i = normIdMap.begin(); i != normIdMap.end(); n++, i++) {
		normIdMap[i->first] = n;
	}

	
	populate_hashmap();


} // End of Manifest::open()



void Manifest::open(char *filename, bool wide)
{
	string f = filename;
	open(f, wide);
}



//////////////////////////////////////////
//
// void Manifest::populate_hashmap ()
//
// Populate the snpNames hash_map for quick look-ups of SNPs.
//
//////////////////////////////////////////

void Manifest::populate_hashmap () {

	//vector<snpClass>::iterator it;
	//for (it = snps.begin(); it != snps.end(); it++) {
	int num_snps = snps.size();
	for (int index = 0; index < num_snps; index++) {
	      snpClass s = snps[index];
	      string the_name = s.name;
	      snpNames[the_name] = index;

	}

}



///////////////////////////////////////////////
//
// void Manifest::convert (snpClass* snip, 
//                         std::string input_snp)
//
// Converts input_snp (the SNP as read in, e.g. 
// "[A/C]") to shorter format (e.g. "AC"), 
// having first changed the SNP - if it is on the
// BOT strand - to the corresponding TOP strand
// SNP. The final format is stored as the
// snp->snp
//
// In the case of any confusion, at least one
// of the two chars in the resulting format
// will be a '?'.
//
/////////////////////////////////////////////

void Manifest::convert (snpClass* snip, std::string input_snp) {

 
  // input_snp will be something like "[A/C]"; we want to store it as "AC"
  // Illumina method aims to designate A as Allele A on TOP, and the T as Alelle A
  // on BOT.
  // 
  // Read as [A/B]    Strand      Store as Allele A, Allele B for TOP
  
  // [C/A]             BOT (B)            GT    
  // [G/A]             BOT                CT
  // [G/C]             BOT                CG
  // [T/G]             BOT                AC 
  // [T/C]             BOT                AG
  // [T/A]             BOT                AT
  
  // So BOT to TOP is C->G, A->T, G->C, T->A

  // [A/C]             TOP (T)            AC
  // [A/G]             TOP                AG
  // [A/T]             TOP                AT
  // [C/G]             TOP                CG
  // [C/T]             TOP                CT
  // [G/T]	           TOP                GT
  
  // weird cases
  // [D/I]             M                  DI
  // [D/I]             P                  DI
  // [I/D]             M                  ID
  // [I/D]             P                  ID
  // [N/A]             P                  NA


  switch ( snip->iStrand ) {

    case 'T': // Already a TOP, or one of the weird cases
    case 'M':
    case 'P':
      snip->snp[0] = input_snp.at(1);
      snip->snp[1] = input_snp.at(3);
      break;

    case 'B':  // BOT
	{
		char first = input_snp.at(1);
		char second = input_snp.at(3);
      
		// This relies on the above table (and conversions) being correct and complete!
		switch ( first ) {
		case 'C':
		  snip->snp[0] = 'G';
		  break;
	  
		case 'G':
		  snip->snp[0] = 'C';
		  break;
	
		case 'T':
		  snip->snp[0] = 'A';
		  break;

		case 'A':
		  snip->snp[0] = 'T';
		  break;

		default:
		  snip->snp[0] = '?';
		}

		switch ( second ) {
		case 'T':
		  snip->snp[1] = 'A';
		  break;
		  
		case 'A':
		  snip->snp[1] = 'T';
		  break;
	  
		case 'C':
		  snip->snp[1] = 'G';
		  break;
		  
		case 'G':
		  snip->snp[1] = 'C';
		  break;
		  
		default:
		  snip->snp[1] = '?';
		}
		snip->converted = true;
	}
	break;
   
  default:   // Unknown
    snip->snp[0] = '?';
    snip->snp[1] = '?';
	
  }

#ifdef _DEBUG
  cout << flush << "\n" << input_snp << " with strand " << snip->iStrand << " produces: " 
       << snip->snp[0] << snip->snp[1] << endl;  
#endif

} // End of Manifest::convert()



///////////////////////////////////////////////
//
// void Manifest::test_convert ()
//
// Use when _DEBUG defined.
// Simple test utility to check convert() works
// OK.
//
///////////////////////////////////////////////

void Manifest::test_convert () {

  cout << "\nstarting test_convert()\n";

  snpClass *mysnp = new snpClass();

  // TOP cases
  mysnp->iStrand = 'T';
  convert (mysnp, "[A/C]");
  convert (mysnp, "[A/G]");
  convert (mysnp, "[A/T]");
  convert (mysnp, "[C/G]");
  convert (mysnp, "[C/T]");
  convert (mysnp, "[G/T]");

  // BOT cases
  mysnp->iStrand = 'B';
  convert (mysnp, "[C/A]");
  convert (mysnp, "[G/A]");
  convert (mysnp, "[G/C]");
  convert (mysnp, "[T/G]");
  convert (mysnp, "[T/C]");
  convert (mysnp, "[T/A]");  
  convert (mysnp, "[A/T]"); // illegal - produces "??"
  convert (mysnp, "[T/T]"); // illegal - produces "A?"

  // weird cases
  mysnp->iStrand = 'M';
  convert (mysnp, "[D/I]");
  convert (mysnp, "[I/D]");
 
  mysnp->iStrand = 'P';
  convert (mysnp, "[D/I]");
  convert (mysnp, "[I/D]");
  convert (mysnp, "[N/A]");



  delete mysnp;

  cout << "\nfinished test_convert()\n";

} // End of Manifest::test_convert()



///////////////////////////////////////////////
//
// string Manifest::get_chromosome_for_SNP (string snpname)
//
// Returns the chromosome which this SNP is on,
// or an empty string if not found.
//
///////////////////////////////////////////////

string Manifest::get_chromosome_for_SNP (string snpname) {


  int index = snpNames[snpname];

  if ( index != 0 ) { // Normal case.

    return snps[index].chromosome;

  }

  else {
    // EITHER SNP is in index 0 of our vector, OR it hasn't
    // been found in our hash_map. So we need to do a special 
    // check for this case, by checking to see if the name of 
    // the SNP in index 0 is the same as the snpname string passed in:
    
    snpClass s = snps[0];
    if ( 0 == s.name.compare(snpname) ) { // It's the same
      return s.chromosome;
    }
    else {
      return string();
    }

  }


}



///////////////////////////////////////////////
//
// snpClass* Manifest::lookup_SNP_by_name (string snpname) 
//
// Given the SNP name snpname, this function
// returns a pointer to a new instance of the 
// snpClass object with this name, or NULL if 
// this isn't found in the Manifest.
//
// Assumes unique SNP names, of course.
//
// Note that modifications to the object whose
// pointer is returned by this function won't affect
// the original object in the Manifest's vector.
//
// Note also that calling code must delete the new
// object when it has finished with it.
//
///////////////////////////////////////////////

snpClass* Manifest::lookup_SNP_by_name (string snpname) {

  int index = snpNames[snpname];
  snpClass s = snps[index];


  //vector<snpClass>::iterator it;
  //for (it = snps.begin(); it != snps.end(); it++) {

    //snpClass s = *it;

    if ( 0 == s.name.compare(snpname) ) {

      snpClass* snippy = new snpClass(s);
      return snippy;

    }

  //}

  return NULL;

}



///////////////////////////////////////////////
//
// snpClass Manifest::lookup_SNP_by_position (long pos)
//
// Given the SNP position pos, this function
// returns a pointer to the snpClass object with this
// position, or NULL if this isn't found in the Manifest.
//
// Assumes unique SNP positions, of course. Note
// that there could be problems using this function
// with CNV regions, as they can overlap.
//
// UPDATE 7th October 2010: The Manifest files do sometimes
// contain more than 1 probe at the same position. 
// Therefore I am removing this function. MG 
//
///////////////////////////////////////////////
/*
  snpClass* Manifest::lookup_SNP_by_position (long pos) {
  
    vector<snpClass>::iterator it;
    for (it = snps.begin(); it != snps.end(); it++) {
  
      snpClass s = *it;
  
      if ( pos == s.position ) {

        snpClass* snippy = new snpClass(s);
        return snippy;

      }

    }

    return NULL;

  }

*/



///////////////////////////////////////////////
//
// void Manifest::dump ()
//
// Produces some output; probably better for 
// client code to write their own versions as
// required.
//
///////////////////////////////////////////////

void Manifest::dump(void)
{
  /*
	map<string,snpClass>::iterator i;
	for (i = snps.begin(); i != snps.end(); i++) {
		snpClass snp = i->second;
		cout << snp.index << ',' 
		     << snp.name << ',' 
			 << snp.chromosome << ',' 
			 << snp.position << ',' 
			 << snp.score << endl;
	}
  */
	map<int,int>::iterator i;
	for (i = normIdMap.begin(); i != normIdMap.end(); i++) {
		cout << i->first << '\t' << i->second << endl;
	}

} 


///////////////////////////////////////////////
//
// void Manifest::write (string outpath)
//
// Write normalized .bpm.csv output to given path
// Use to create a normalized .bpm.csv file for input to genotype callers
//
///////////////////////////////////////////////

void Manifest::write(string outPath) {

  ofstream outFile;
  outFile.open(outPath.c_str());


  outFile << "BPM output goes here.\n";

  // Write SNP data to file
  int snpTotal = snps.size(); // Number of elements.

  // TODO ensure SNPs are output in sorted order?
  for (int i = 0; i < snpTotal; i++) {
	snpClass s = snps[i];
	outFile << s.toString().c_str() << endl;

  }


  outFile.close();


}



///////////////////////////////////////////////
//
// void Manifest::order_by_position()
//
// Re-order the vector of snpClass objects such that 
// they are stored in ascending position order.
//
// Must be called after Manifest::open().
//
///////////////////////////////////////////////

void Manifest::order_by_position() {

  // Multimap allows duplicate entries having the same key. So it
  // is fine if the Manifest file contains multiple probes at the
  // same position.
  multimap<long, snpClass> posmap; // Key = position.

  // Populate the map.
  int siz = snps.size(); // Number of elements.
  for (int i = 0; i < siz; i++) {
	snpClass s = snps[i];
	long pos = s.position;
	//posmap[pos] = s;
	posmap.insert( pair<long,snpClass>(pos, s) );
  }

  // Clear out old contents of vector:
  snps.clear();

  // Now use the map, which stores elements in key (pos) order,
  // to order vector by position
  int mapsize = posmap.size();
  if ( mapsize != siz ) {
	cout << "\nERROR in Manifest::order_vector_by_position(): mapsize = " << mapsize
         << " but siz = " << siz << "\n" << flush;
	exit(1);
  }
  
  multimap<long, snpClass>::iterator posit;
  for ( posit = posmap.begin(); posit != posmap.end(); posit++ ) {
	snpClass s = (*posit).second;
	snps.push_back(s);
  }


  populate_hashmap();

}



///////////////////////////////////////////////
//
// void Manifest::exclude_cnvs()
//
// Sets EXCLUE_CNVS to false; this means that
// CNV probes will not be stored by the Manifest
// instance.
//
// If required, must be called before Manifest::open().
//
///////////////////////////////////////////////

void Manifest::exclude_cnvs() {

  EXCLUDE_CNVS = true;

}



///////////////////////////////////////////////
//
// bool Manifest::find_wide_columns(...)
//
// Given the data header line 'linestr' of a "wide-format" (.csv) Manifest
// file, discover which header corresponds to which column
// number of the header and store in map 'dict'.
//
// Header line will be something like this (wrapped for clarity):
//
// IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,
// AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,
// SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,
// Intensity_Only,Exp_Clusters,CNV_Probe
//
// Hence dict[0] = "IlmnID", dict[1] = "Name", etc.
//
///////////////////////////////////////////////

bool Manifest::find_wide_columns(string& linestr, map<string,int>& dict) {

  //cout << "\nDBG: find_wide_columns(): start\n" << flush;

  // Remove any newline or \r characters etc at end.
  size_t lastpos = linestr.find_last_of("\r\n");
  //cout << " folder: " << str.substr(0,found) << endl;
  string line = linestr.substr(0, lastpos);

  string delim = ",";
  
  int len = line.length();
  size_t start = 0;
  size_t found = line.find(delim, start);
  int colnum = 0;

  while ( found != string::npos ) {
	
	string candidate = line.substr(start, found-start);
	if ( ! candidate.empty() ) {
	  dict[candidate] = colnum;
	}
	start = found + delim.length();
	found = line.find(delim, start);		
	colnum++;
  }
  
  // Get last (or only!) candidate:
  string cand = line.substr(start, len - start);
  if ( ! cand.empty() ) {
	  dict[cand] = colnum;
	}


  /*
    // Debug - show what we have:
    map<string, int>::iterator it;
    for ( it=dict.begin(); it != dict.end(); it++ ) {
	  cout << "\n" << (*it).first << " => " << (*it).second << flush;
    }
    cout << "\n";
    exit(1);
  */

  return true;
}




///////////////////////////////////////////////
//
// int Manifest::get_map_value (...)
//
// Given a <string,int> map 'mymap', and a look-up 
// string 'treasure', find the int valule stored
// as the value corresponding to thsi string,
// or -1 if not present.
//
///////////////////////////////////////////////

int Manifest::get_map_value (map<string, int>& mymap, const char* const treasure) {
  
  // We can't look for the item in mymap using [] syntax as that creates
  // a new entry!

  int answer;

  map<string, int>::iterator it;

  it = mymap.find(treasure);
  if ( it == mymap.end() ) {
	answer = -1;
  }
  else {
	answer = (*it).second;
  }


  // cout << "\nDBG: get_wide_col() for " << treasure << ": " << answer << flush;

  return answer;


}


///////////////////////////////////////////////
//
// string snpClass::toString()
//
// method of snpClass objects to convert data to a comma-delimited string
// use for .bpm.csv output in Manifest::write
//
///////////////////////////////////////////////

string snpClass::toString() {

  // placeholders for char array to C++ string conversion
  char buffer[100];
  int n;

  // index
  n = sprintf(buffer, "%d", index);
  string i = string(buffer, n);

  // position
  n = sprintf(buffer, "%ld", position);
  string p = string(buffer, n);

  // score
  n = sprintf(buffer, "%f", score);
  string s = string(buffer, n);


  // generate the allele string
  string a = string(1, snp[0]); // this constructor makes 1 copy of snp[0]
  string b = string(1, snp[1]);
  string alleles;
  if (a=="?" or b=="?") {
    alleles = "[N/N]";
  } else {
    alleles = "["+a+"/"+b+"]";
  }
  
  // strands after normalization
  string iStrandNorm, cStrandNorm; 
  iStrandNorm = strandToString(iStrand, converted);
  cStrandNorm = strandToString(cStrand, converted);

  // norm ID
  n = sprintf(buffer, "%d", normId);
  string nid = string(buffer, n);

  string snpString = i+","+name+","+chromosome+","+s+","+alleles+\
    ","+p+","+iStrandNorm+","+cStrandNorm+","+nid;

  /*
    TODO work out why compiler complains:  no member named ‘BeadSetId’

  if (this->BeadSetID != -1) {
     n = sprintf(buffer, "%d", this->BeadSetId);
     string bsid = string(buffer, n);
     snpString = snpString+","+bsid;
  }

  */

  cout << snpString << endl;
  exit(1);

  return snpString;

}


///////////////////////////////////////////////
//
// string snpClass::strandToString(char strand, bool converted)
//
// Convert char strand ID to original long-form string
// use in snpClass::toString
//
///////////////////////////////////////////////


string snpClass::strandToString(char strand, bool converted) {
  // generate string to represent strand
  // may be normalized to Illumina top strand
  string normStrand;
  if (converted && strand=='B') { strand='T'; }
  if (strand=='T') { normStrand = "TOP"; }
  else if (strand=='B') { normStrand = "BOT"; }
  else if (strand=='M') { normStrand = "MINUS"; }
  else if (strand=='P') { normStrand = "PLUS"; }
  else { normStrand = "UNKNOWN"; }
  return normStrand;
}


// EOF
