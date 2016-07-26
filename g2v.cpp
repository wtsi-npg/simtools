//
// g2v.cpp
//
// Program to convert GTC format genotypes to VCF.
//
// Copyright (C) 2016 Genome Research Ltd.
//
// Author: Keith James
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer
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

#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>

#include "Gtc.h"
#include "Manifest.h"

using namespace std;

string VCF_VERSION = "4.2";

int CLI_OPTIONS_ERR = 2;
int    MANIFEST_ERR = 3;
int       GTC_ERROR = 4;

bool verbose = false;

string gtc_file;
string manifest_file;
string reference_name;

int GENCALL_THRESHOLD = 15; // 0.15 * 100

void print_usage() {
  cout << "g2v "
       << "-g <GTC file> "
       << "[-h] "
       << "-m <manifest file> "
       << "[-r <reference name>] "
       << "[> <VCF file>]" << endl;
}

char *timestamp(void) {
  static char buffer[64];
  time_t t = time(NULL);
  strftime(buffer, 64, "%F", localtime(&t));

  return buffer;
}

// Remove indel probes and anything with a failed GenCall score
void filter_fails_and_indels(Manifest & manifest, Gtc & gtc,
                             vector<snpClass> & called_snps) {
  for (auto si = manifest.snps.begin(); si != manifest.snps.end(); si++) {
    if (si->snp[0] != 'D' && si->snp[0] != 'I') {

      int i = si->index - 1;
      float score = gtc.scores[i];
      if (score * 100 > GENCALL_THRESHOLD) {
        called_snps.push_back(*si);
      }
    }
  }
  sort(called_snps.begin(), called_snps.end());
}

void print_boilerplate_fields() {
  cout << "##fileformat=VCFv" << VCF_VERSION
       << endl;

  cout << "##fileDate=" << timestamp()
       << endl;

  if (reference_name.size() > 0) {
    cout << "##reference=" << reference_name
         << endl;
  }

  cout << "##INFO=<ID=NS,"
       << "Number=1,Type=Integer,"
       << "Description=\"Number of samples with data\">"
       << endl;

  cout << "##INFO=<ID=AC,"
       << "Number=.,Type=Integer,Description=\"Allele count\">"
       << endl;

  cout << "##INFO=<ID=AN,"
       << "Number=1,Type=Integer,"
       << "Description=\"Number of alleles with data\">"
       << endl;

  cout << "##FILTER=<ID=gencall,Description=\"Illumina GenCall score\">"
       << endl;

  cout << "##FORMAT=<ID=GT,"
       << "Number=1,Type=String,"
       << "Description=\"Genotype\">"
       << endl;

  return;
}

void collect_chr_names(Manifest & manifest, vector<string> & chr_names) {
  set<string> cnames;
  for (auto si = manifest.snps.begin(); si != manifest.snps.end(); si++) {
    cnames.insert(si->chromosome);
  }

  chr_names.assign(cnames.begin(), cnames.end());
  sort(chr_names.begin(), chr_names.end());

  return;
}

void collect_unique_alts(Gtc & gtc, vector<snpClass> & snps,
                         vector<char> & alts) {
  for (auto si = snps.begin(); si != snps.end(); si++) {
    int i = si->index - 1;
    int gt_code = gtc.genotypes[i];
    if (gt_code == 2 or gt_code == 3) {
      alts.push_back(si->snp[1]);
    }
  }
  sort(alts.begin(), alts.end());
  alts.erase(unique(alts.begin(), alts.end()), alts.end());
}

void collect_unique_identifiers(vector<snpClass> & snps,
                                vector<string> & ids) {
  for (auto si = snps.begin(); si != snps.end(); si++) {
    ids.push_back(si->name);
  }

  sort(ids.begin(), ids.end());
  ids.erase(unique(ids.begin(), ids.end()), ids.end());
}

void format_identifiers(vector<string> & identifiers, string & field) {
  for (auto ni = identifiers.begin(); ni != identifiers.end(); ni++) {
    field += *ni;

    if (next(ni) != identifiers.end()) {
      field += ";";
    }
  }
}

void format_alts(vector<char> & alts, string & field) {
  if (alts.empty()) {
    field = ".";
  }
  else {
    for (auto ai = alts.begin(); ai != alts.end(); ai++) {
      field += *ai;

      if (next(ai) != alts.end()) {
        field += ",";
      }
    }
  }
}

void format_info(int num_called, vector<char> & alts, string & field) {
  field += "NS=1";
  field += ";AC=" + to_string(alts.size());
  field += ";AN=" + to_string(num_called);
}

void format_genotypes(vector<string> & genotypes,
                      string & format_field, string & genotypes_field) {
  for (auto gi = genotypes.begin(); gi != genotypes.end(); gi++) {
    format_field    += "GT";
    genotypes_field += *gi;

    if (next(gi) != genotypes.end()) {
      format_field    += ":";
      genotypes_field += ":";
    }
  }
}

void print_chr_name_fields(vector<string> & chr_names) {
  for (auto ni = chr_names.begin(); ni != chr_names.end(); ni++) {
    cout << "##contig=<ID=" << *ni << ">" << endl;
  }

  return;
}

void print_data_header(string sample_name) {
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
       << sample_name
       << endl;

  return;
}

void print_locus_records(Gtc & gtc, vector<snpClass> & snps) {
  vector<string> identifiers;
  collect_unique_identifiers(snps, identifiers);

  vector<char> alts;
  collect_unique_alts(gtc, snps, alts);

  // Map sorted alt bases to an alt index number
  map<char, int> alt_nums;
  for (auto ai = alts.begin(); ai != alts.end(); ai++) {
    alt_nums[*ai] = alt_nums.size();
  }

  vector<string> genotypes;
  int n = 0; // Total number of called alleles

  // Calculate data fields 
  for (auto si = snps.begin(); si != snps.end(); si++) {
    int i = si->index - 1;
    int gt_code = gtc.genotypes[i];

    string gt_ref = "0";
    string gt_alt;
    if (gt_code == 2 or gt_code == 3) {
      gt_alt = to_string(alt_nums[si->snp[1]]);
    }

    string genotype;
    switch (gt_code) {
        case 0: genotype += ".";                           break;
        case 1: genotype += gt_ref + "/" + gt_ref; n += 2; break;
        case 2: genotype += gt_ref + "/" + gt_alt; n += 2; break;
        case 3: genotype += gt_alt + "/" + gt_alt; n += 2; break;

        default:
          cerr << "Invalid genotype code " << gt_code
               << " for index " << i << endl;
          exit(GTC_ERROR);
    }

    genotypes.push_back(genotype);
  }

  // VCF spec section 1.4.1.3
  string identifiers_field;
  format_identifiers(identifiers, identifiers_field);

  // VCF spec section 1.4.1.5
  string alts_field;
  format_alts(alts, alts_field);

  // VCF spec section 1.4.1.8
  string info_field;
  format_info(n, alts, info_field);

  string format_field;
  string genotypes_field;
  format_genotypes(genotypes, format_field, genotypes_field);

  cout << snps.at(0).chromosome << "\t"
       << snps.at(0).position   << "\t"
       << identifiers_field     << "\t"
       << snps.at(0).snp[0]     << "\t"
       << alts_field            << "\t"
       << "."                   << "\t"
       << "PASS"                << "\t"
       << info_field            << "\t"
       << format_field          << "\t"
       << genotypes_field       << endl;
}

int main(int argc, char *argv[]) {

  char c;
  while ((c = getopt (argc, argv, "g:hm:r:v")) != -1) {
    switch (c) {
        case 'g': gtc_file       = optarg; break;
        case 'h': print_usage(); exit(0);  break;
        case 'm': manifest_file  = optarg; break;
        case 'r': reference_name = optarg; break;
        case 'v': verbose = true;          break;
    }
  }

  if (gtc_file.size() == 0) {
    print_usage();
    cerr << "No GTC file specified" << endl;
    exit(CLI_OPTIONS_ERR);
  }
  if (manifest_file.size() == 0) {
    print_usage();
    cerr << "No manifest file specified" << endl;
    exit(CLI_OPTIONS_ERR);
  }

  Manifest *manifest = new Manifest();
  try {
    if (verbose) {
      cout << "Using manifest file " << manifest_file << endl;
    }
    manifest->open(manifest_file);
    manifest->order_by_locus();
  }
  catch (string s) {
    cerr << s << endl;
    exit(MANIFEST_ERR);
  }

  Gtc gtc;
  gtc.open(gtc_file, Gtc::GENOTYPES | Gtc::BASECALLS | Gtc::SCORES);

  vector<string> chromosomes;
  collect_chr_names(*manifest, chromosomes);

  vector<snpClass> called_snps;
  filter_fails_and_indels(*manifest, gtc, called_snps);

  print_boilerplate_fields();
  print_chr_name_fields(chromosomes);
  print_data_header(gtc.sampleName);

  // There are multiple probes at the same locus to be combined into a
  // single VCF record
  vector<snpClass> accum;
  for (auto si = called_snps.begin(); si != called_snps.end(); si++) {
    // Collect first result
    if (si == called_snps.begin()) {
      accum.push_back(*si);
      continue;
    }

    if (si->chromosome.compare(prev(si, 1)->chromosome) == 0 &&
        si->position == prev(si, 1)->position) {
      // Collect current locus
      accum.push_back(*si);
    }
    else {
      // Process current locus and clear
      print_locus_records(gtc, accum);
      accum.clear();

      // Collect new current locus
      accum.push_back(*si);
    }
  }

  if (!accum.empty()) {
    // Print final locus
    print_locus_records(gtc, accum);
  }

  // The destructor causes a segfault by attempting to free something
  // it should not.

  // delete manifest;

  return 0;
}
