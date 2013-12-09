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

#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <cxxtest/TestSuite.h>

// include statements from simtools
#include "commands.h"
#include "Sim.h"
#include "Gtc.h"
#include "QC.h"
#include "Manifest.h"
#include "json/json.h"

using namespace std;


class TestBase :  public CxxTest::TestSuite
{
  // base class with shared methods

 public:

  string tempdir;

  void setUp() 
  {
    // called before every test
    cerr << endl;
    char *_template = new char[20]; 
    strcpy(_template,  "temp_XXXXXX"); // template string for directory name
    tempdir = string(mkdtemp(_template)); // returns name of new directory
    TS_TRACE("TEMPDIR:"+tempdir);
    
  }

  void tearDown()
  {
    // this is called after every test (successful or not)
    string cmd = "rm -Rf " + tempdir;
    //int result = 0;
    int result = system(cmd.c_str());
    cerr << cmd << " " << result << endl;
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
    TS_TRACE(infile);
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


class SimtoolsTest : public TestBase
{


 public:


  void testCreate(void) {
    // duplicate ../simtools create --infile example.json --outfile test.sim --man_file example.bpm.csv
    // TODO revise example.json to have absolute, or correct relative, paths


    TS_TRACE("Testing .sim create command");

    string infile = "data/example_with_dir.json";
    string expected = "data/example.raw.sim";
    string outfile = tempdir+"/test.sim";
    bool normalize = false;
    string manfile = "data/example.bpm.csv";
    bool verbose = false;

    Commander *commander = new Commander();
    TS_ASSERT_THROWS_NOTHING(commander->commandCreate(infile, outfile, normalize, manfile, verbose));
    TS_TRACE("SIM file successfully created from GTC");
    int size = 1491; // expected file size
    assertFileSize(outfile, size);
    TS_TRACE("SIM file created from GTC is of expected length");
    assertFilesIdentical(outfile, expected, size);
    TS_TRACE("SIM file created from GTC is identical to master");
    delete commander;

  }

  void testView(void) {
    TS_TRACE("Testing .sim view command");
    Commander *commander = new Commander();
    string simfile = "./data/example.raw.sim";
    bool verbose = false;
    TS_ASSERT_THROWS_NOTHING(commander->commandView(simfile, verbose));
    delete commander;
    TS_TRACE("View command test finished");

  }



};


// Putting TestSuite classes in separate files appears not to work
// See Cxx manual section 4.4; possible issue with compiler options
