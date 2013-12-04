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

#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <cxxtest/TestSuite.h>
#include "Manifest.h"

using namespace std;


class NormalizeTest : public CxxTest::TestSuite
{
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

  }
  
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

    // read test output
    ifstream testStream;
    TS_ASSERT_THROWS_NOTHING(testStream.open(outfile.c_str()));
    int size = 1275; // expected file size
    testStream.seekg(0, testStream.end);
    int testSize = testStream.tellg();
    testStream.seekg(0, testStream.beg);
    TS_ASSERT_EQUALS(size, testSize);
    TS_TRACE("Normalized .csv file is of correct length");
    char *testBuf = new char[size];   
    testStream.read(testBuf, size);
    testStream.close();

    // read master normalized file
    ifstream normStream;
    normStream.open(normfile.c_str());
    char *normBuf = new char[size];   
    normStream.read(normBuf, size);
    TS_ASSERT_SAME_DATA(testBuf, normBuf, size);
    TS_TRACE("Normalized .csv file is identical to master");
    normStream.close();

    delete manifest;
  }

};


class SimtoolsTest : public CxxTest::TestSuite
{


 public:

  void testTrace(void) {
    // simple "Hello world" test of tests
    TS_TRACE("Testing simtools");
    int foo = 0;
    TS_ASSERT_EQUALS(foo, 0);

  }



};


// Putting TestSuite classes in separate files appears not to work
// See Cxx manual section 4.4
// May be an issue with compiler options for simtools?
