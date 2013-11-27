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

#include <stdio.h>
#include <stdlib.h>

#include <cxxtest/TestSuite.h>
#include "Manifest.h"

using namespace std;


class NormalizeTest : public CxxTest::TestSuite
{
 public:

  char *tempdir;

  void setUp() 
  {
    // called before every test
    cerr << endl;
    tempdir = new char[20]; 
    strcpy(tempdir,  "test_XXXXXX"); // template string for directory name
    tempdir = mkdtemp(tempdir); // returns name of new directory
    cerr << tempdir << endl;
    
  }


  void tearDown()
  {
    // this is called after every test (successful or not)
    string cmd = "rm -Rf " + string(tempdir);
    int result = system(cmd.c_str());
    cerr << cmd << " " << result << endl;

  }
  
  // example tests taken from Cxxtest docs

  void testAddition(void)
  {
    TS_ASSERT(1 + 1 > 1);
    TS_ASSERT_EQUALS(1 + 1, 2);
  }


  void testMultiplication(void)
  {
    TS_TRACE("Starting multiplication test");
    TS_ASSERT_EQUALS(2 * 2, 5);
    TS_TRACE("Finishing multiplication test");
  }
};
