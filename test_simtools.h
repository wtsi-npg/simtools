// Author: Iain Bancarz <ib5@sanger.ac.uk>
//
// Copyright (c) 2013 Genome Research Ltd. All rights reserved.
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


// Test classes for simtools in cxxtest framework

#include <cerrno>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cxxtest/TestSuite.h>
#include "commands.h"
#include "Manifest.h"
#include "unistd.h"
#include "win2unix.h"

using namespace std;


class TestBase :  public CxxTest::TestSuite
{
  // base class with shared methods

 public:

  static const int sim_size = 1491; // size of data/example.raw.sim 
  string tempdir;
  string sim_raw;
  string manfile;
  bool verbose;

  void setUp() 
  {
    // called before every test
    char *_template = new char[20]; 
    strcpy(_template,  "temp_XXXXXX"); // template string for directory name
    tempdir = string(mkdtemp(_template)); // returns name of new directory
    TS_TRACE("Created tempdir:"+tempdir);
    string cmd = "/bin/cp data/example.raw.sim data/example.bpm.csv "+tempdir;
    int status = system(cmd.c_str());
    if (status!=0) { // test assertions are not allowed in setup
      cerr << "Failed to copy test data to temporary directory: ";
      cerr << tempdir << endl;
      throw 1;
    }
    sim_raw = tempdir+"/example.raw.sim";
    manfile = tempdir+"/example.bpm.csv";
    verbose = false;    
  }

  void tearDown()
  {
    // this is called after every test (successful or not)
    string cmd = "/bin/rm -Rf " + tempdir;
    int status = system(cmd.c_str());
    if (status!=0) { // test assertions are not allowed in teardown
      cerr << "Failed to delete temporary directory: " << tempdir << endl;
      throw 1;
    }
    TS_TRACE("Removed tempdir: "+tempdir);
  }

  void assertFilesIdentical(string path1, string path2, int size)
  {
    // check if given files contain same data; size = expected size in bytes
    char *buf1 = new char[size];   
    char *buf2 = new char[size];   
    ifstream stream1;
    ifstream stream2;
    stream1.open(path1.c_str());
    stream1.read(buf1, size);
    stream1.close();
    stream2.open(path2.c_str());
    stream2.read(buf2, size);
    stream2.close();
    TS_ASSERT_SAME_DATA(buf1, buf2, size);
    delete buf1;
    delete buf2;
  }

  void assertFileSize(string path, int size)
  {
    ifstream testStream;
    testStream.open(path.c_str());
    testStream.seekg(0, testStream.end);
    int testSize = testStream.tellg();
    testStream.close();
    TS_ASSERT_EQUALS(size, testSize);
  }

  int stdoutRedirect(string tempfile) {

    fflush(stdout);
    int sout = dup(fileno(stdout));
    if (sout == -1) {
      cerr << "Failed to dup stdout: " << string(strerror(errno)) << endl;
      throw 1;
    }

    // redirect standard output to given path, which may be /dev/null
    FILE *ptr;
    ptr = freopen(tempfile.c_str(), "w", stdout);
    if (!ptr) { // null pointer
      cerr << "Failed to redirect standard output to file: ";
      cerr << tempfile << endl;
      throw 1;
    }

    return sout;
  }

  void stdoutRestore(int sout) {
    // restore stdout to standard terminal

    int fd = dup2(sout, fileno(stdout));
    if (fd == -1) {
      //cerr << "Failed to restore stdout: " << string(strerror(errno)) << endl;
      throw 1;
    }

    close(sout);
    clearerr(stdout);
  }
};


class NormalizeTest : public TestBase
{
 public:

  void testManifest(void)
  {
    // test creation of Manifest objects
    string infile = "data/mock_1000.bpm.csv";
    string outfile = tempdir+"/normalized.bpm.csv";
    Manifest *manifest;
    TS_TRACE("Starting manifest test");
    manifest = new Manifest();
    TS_ASSERT_THROWS_NOTHING(manifest->open(infile, "1", false)); // chr1 only
    TS_ASSERT_EQUALS(manifest->snps.size(), 500);
    TS_TRACE("Read manifest for chromosome 1 only");
    delete manifest;
    manifest = new Manifest();
    TS_ASSERT_THROWS_NOTHING(manifest->open(infile));
    TS_ASSERT_EQUALS(manifest->snps.size(), 1000);
    TS_ASSERT_THROWS_NOTHING(manifest->write(outfile));
    TS_TRACE("Read manifest and write normalized for all chromosomes");
    delete manifest;
    TS_TRACE("Finished manifest test");
  }

  void testNormalize(void)
  {
    // compare normalized output with reference file
    string infile = "data/mock.bpm.csv";
    string normfile = "data/mock_normalized.bpm.csv";
    string outfile = tempdir+"/mock_normalized.bpm.csv";
    Manifest *manifest = new Manifest();
    TS_TRACE("Normalizing "+infile);
    TS_ASSERT_THROWS_NOTHING(manifest->open(infile));
    TS_ASSERT_THROWS_NOTHING(manifest->write(outfile));

    int size = 1275; // expected file size
    assertFileSize(outfile, size);
    TS_TRACE("Normalized .csv file is of correct length");

    // compare output data
    assertFilesIdentical(normfile, outfile, size);
    TS_TRACE("Normalized .csv file is identical to master");

    delete manifest;
  }

};

class SimTest : public TestBase {

 public:

  void testExecutable(void) {
    TS_TRACE("Testing the 'sim' executable");
    // suppress standard output form command
    string cmd = "./sim "+sim_raw+" 10 >/dev/null";
    TS_ASSERT_EQUALS(system(cmd.c_str()), 0);
  }

};


class SimtoolsTest : public TestBase
{

 public:

  void testCreate(void) {
    TS_TRACE("Testing .sim create command");
    string infile = "data/example.json";
    string expected = "data/example.raw.sim";
    string outfile = tempdir+"/test.sim";
    bool normalize = false;
    Commander *commander = new Commander();
    TS_ASSERT_THROWS_NOTHING(commander->commandCreate(infile, outfile, normalize, manfile, verbose));
    delete commander;
    TS_TRACE("SIM file successfully created from GTC");
    assertFileSize(outfile, sim_size);
    TS_TRACE("SIM file created from GTC is of expected length");
    assertFilesIdentical(outfile, sim_raw, sim_size);
    TS_TRACE("SIM file created from GTC is identical to master");
  }

  void testGenoSNP(void) {
    // Tests of GenoSNP mode:
    // 1. Input from file, output all SNPs
    // 2. Input from stdin, output all SNPs
    // 3. Input from file, output subset of SNPs
    int size_all = 430;
    int size_single = 168;
    Commander *commander = new Commander();
    int start_pos = 0;
    int end_pos = -1;
    string outfile1 = tempdir+"/genosnp01.gsn";
    string expected = "data/example_all.gsn";
    TS_TRACE("Testing GenoSNP command with .sim input from file");
    TS_ASSERT_THROWS_NOTHING(commander->commandGenoSNP(sim_raw, outfile1, manfile, start_pos, end_pos, verbose));
    assertFileSize(outfile1, size_all);
    assertFilesIdentical(outfile1, expected, size_all);

    // now try with standard input; TODO do this without a system call
    string outfile2 = tempdir+"/genosnp02.gsn";
    char cmd[1024];
    snprintf(cmd, sizeof cmd,
             "./simtools genosnp --infile - --outfile %s --man_file %s < %s",
             outfile2.c_str(), manfile.c_str(), sim_raw.c_str());

    TS_TRACE("Testing GenoSNP command with .sim input from STDIN");
    TS_ASSERT_EQUALS(system(cmd), 0);
    assertFileSize(outfile2, size_all);
    assertFilesIdentical(outfile2, expected, size_all);

    // test output of a single sample
    start_pos = 2;
    end_pos = 3;
    string outfile3 = tempdir+"/genosnp03.gsn";
    expected = "data/example_single.gsn";
    TS_TRACE("Testing GenoSNP command with output of a single sample");
    TS_ASSERT_THROWS_NOTHING(commander->commandGenoSNP(sim_raw, outfile3, manfile, start_pos, end_pos, verbose));
    delete commander;
    assertFileSize(outfile3, size_single);
    assertFilesIdentical(outfile3, expected, size_single);
  }

  void testIlluminus(void) {
    // Tests of Illuminus mode:
    // 1. Input from file, output all SNPs
    // 2. Input from stdin, output all SNPs
    // 3. Input from file, output subset of SNPs
    int size_all = 1268; // output size for all SNPs
    int size_single = 349; // SNP 3 only
    Commander *commander = new Commander();
    int start_pos = 0;
    int end_pos = -1;
    string outfile1 = tempdir+"/illuminus01.iln";
    string expected = "data/example_all.iln";
    TS_TRACE("Testing Illuminus command with .sim input from file");
    TS_ASSERT_THROWS_NOTHING(commander->commandIlluminus(sim_raw, outfile1, manfile, start_pos, end_pos, verbose));
    assertFileSize(outfile1, size_all);
    assertFilesIdentical(outfile1, expected, size_all);

    // now try with standard input
    string outfile2 = tempdir+"/illuminus02.iln";
    // TODO Would be cleaner to run test without the system call
    TS_TRACE("Testing Illuminus command with .sim input from STDIN");
    char cmd[1024];
    snprintf(cmd, sizeof cmd,
             "./simtools illuminus --infile - --outfile %s --man_file %s < %s",
             outfile2.c_str(), manfile.c_str(), sim_raw.c_str());
    TS_ASSERT_EQUALS(system(cmd), 0);
    assertFileSize(outfile2, size_all);
    assertFilesIdentical(outfile2, expected, size_all);

    // input a subset of SNPs
    TS_TRACE("Testing Illuminus command with output of a single SNP");
    start_pos = 3;
    end_pos = 3;
    string outfile3 = tempdir+"/illuminus03.iln";
    expected = "data/example_single.iln";
    TS_ASSERT_THROWS_NOTHING(commander->commandIlluminus(sim_raw, outfile3, manfile, start_pos, end_pos, verbose));
    delete commander;
    assertFileSize(outfile3, size_single);
    assertFilesIdentical(outfile3, expected, size_single);
  }

  void testQC(void) {
    TS_TRACE("Testing QC command");
    string infile = "data/example.larger_with_inf.sim";
    string mag = tempdir+"/mag.txt";
    string xyd = tempdir+"/xyd.txt";
    string mag_expected = "data/mag.txt";
    string xyd_expected = "data/xyd.txt";
    bool verbose = false;
    Commander *commander = new Commander();
    TS_ASSERT_THROWS_NOTHING(commander->commandQC(infile, mag, xyd, verbose));
    delete commander;
    TS_TRACE("QC files successfully created from SIM");
    int mag_size = 4500;
    assertFileSize(mag, mag_size);
    TS_TRACE("QC magnitude output is of expected size");
    assertFilesIdentical(mag, mag_expected, mag_size);
    TS_TRACE("QC magnitude output is identical to master");
    int xyd_size = 4599;
    assertFileSize(xyd, xyd_size);
    TS_TRACE("QC xydiff output is of expected size");
    assertFilesIdentical(xyd, xyd_expected, xyd_size);
    TS_TRACE("QC xydiff output is identical to master");
  }

  void testView(void) {
    TS_TRACE("Testing .sim view command");
    string simfile = "./data/example.raw.sim";
    string tempfile = tempdir+"/simview.txt";
    string viewfile = "data/simview.txt"; // expected view output
    int viewsize = 355;
    bool verbose = false;

    int sout = stdoutRedirect(tempfile);

    Commander *commander = new Commander();
    TS_ASSERT_THROWS_NOTHING(commander->commandView(simfile, verbose));
    delete commander;

    assertFileSize(tempfile, viewsize);
    TS_TRACE("Sim view output is of expected size");
    assertFilesIdentical(tempfile, viewfile, viewsize);
    TS_TRACE("Sim view output is identical to master");

    stdoutRestore(sout);
  }
};


class Win2UnixTest : public TestBase
{
 public:

  void testWin2Unix(void) {

    string filepath = "\\\\fastnfs\\illumina_geno4\\call\\20130605\\9298751015_R01C01.gtc";
    string expected = "/nfs/new_illumina_geno04/call/20130605/9298751015_R01C01.gtc";
    filepath = win2unix(filepath);
    TS_ASSERT_EQUALS(filepath, expected);
   
  }

};

// Putting TestSuite classes in separate files appears not to work
// See Cxx manual section 4.4; possible issue with compiler options
