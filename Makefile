
#
# Makefile for Gtc classes
#
# Author: Jennifer Liddle (js10)
#
# $Id: Makefile 1512 2013-02-11 17:23:41Z js10 $
#

# Author: Jennifer Liddle <js10@sanger.ac.uk, jennifer@jsquared.co.uk>
#

# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
# 3. Neither the name of Genome Research Ltd nor the names of the 
# contributors may be used to endorse or promote products derived from 
# software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR 
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL GENOME RESEARCH LTD. BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



.PHONY: test # always run, regardless of timestamps

PREFIX=/usr/local
INSTALL_INC=$(PREFIX)/include
INSTALL_LIB=$(PREFIX)/lib
INSTALL_BIN=$(PREFIX)/bin

# other useful paths
STLPORT_INC=/software/solexa/pkg/STLport/current/stlport
STLPORT_LIB=/software/solexa/pkg/STLport/current/build/lib/obj/gcc/so

EXECUTABLES=gtc g2i gtc_process sim simtools normalize_manifest
INCLUDES=Sim.h Gtc.h Manifest.h win2unix.h 
LIBS=libsimtools.so libsimtools.a
PERL_MODULES = Gtc.pm Sim.pm
TARGETS=$(EXECUTABLES) $(LIBS)
PERL_TARGETS=$(PERL_MODULES)

PERL_CC_OPTS=$(shell perl -MExtUtils::Embed -e ccopts)
PERL_LD_OPTS=$(shell perl -MExtUtils::Embed -e ldopts)

# To compile with debug information added, invoke as (say): 
# make DEBUG='y'
# to build all targets with debug info. 
# For just one target, say:
# make DEBUG='y' simtools

# do NOT use -ffast-math, as it causes errors in infinity/NaN handling
ifeq ($(DEBUG),y)
	CXXFLAGS+=-g -O0
else
	CXXFLAGS+=-O3
endif

CXXFLAGS+=-Wall -fPIC -I$(STLPORT_INC) -std=c++0x

# Set runpath instead of relying on LD_LIBRARY_PATH
LDFLAGS=-L./ -L$(STLPORT_LIB) -Wl,-rpath -Wl,$(STLPORT_LIB)


default: all

usage:
	@echo -e "Usage: make install DIR=<destination directory>\nOther targets: make all, make test, make install_g2i"

clean:
	rm -f *.o *.so Gtc_wrap.cxx Gtc.pm Sim_wrap.cxx Sim.pm runner.cpp runner $(TARGETS)

test: Sim.o Gtc.o Manifest.o QC.o win2unix.o json/json_reader.o json/json_writer.o json/json_value.o commands.o runner.o
	$(CXX) $(CXXFLAGS) -Wno-deprecated $(LDFLAGS) -o runner $^ -lstlport
	LD_LIBRARY_PATH=. ./runner # run "./runner -v" to print trace information

runner.cpp: test_simtools.h
	cxxtestgen --error-printer -o runner.cpp test_simtools.h

runner.o: runner.cpp all
	$(CXX) -c $(CXXFLAGS) -Wno-deprecated $(LDFLAGS) -o $@ runner.cpp

install: all
	@echo "Installing to "$(PREFIX)
	install -d $(INSTALL_INC) $(INSTALL_LIB) $(INSTALL_BIN)
	install $(INCLUDES) $(INSTALL_INC)
	install $(LIBS) $(INSTALL_LIB)
	install $(PERL_MODULES:pm=so) $(INSTALL_LIB)
	install $(PERL_MODULES) $(INSTALL_LIB)
	install $(EXECUTABLES) $(INSTALL_BIN)

all: $(TARGETS) $(PERL_MODULES)

perl: $(PERL_MODULES)

gtc: gtc.o libsimtools.a
	$(CXX) $< $(LDFLAGS) -o $@ -lm -lstlport -Wl,-Bstatic -lsimtools -Wl,-Bdynamic

normalize_manifest: normalize_manifest.o libsimtools.a
	$(CXX) $< $(LDFLAGS) -o $@ -lm -lstlport -Wl,-Bstatic -lsimtools -Wl,-Bdynamic

sim: sim.o libsimtools.a
	$(CXX) $< $(LDFLAGS) -o $@ -lm -lstlport -Wl,-Bstatic -lsimtools -Wl,-Bdynamic

simtools: simtools.o commands.o libsimtools.a
	$(CXX) simtools.o commands.o $(LDFLAGS) -o $@ -lm -lstlport -Wl,-Bstatic -lsimtools -Wl,-Bdynamic

g2i: g2i.o libsimtools.a
	$(CXX) $< $(LDFLAGS) -o $@ -lm -lstlport -Wl,-Bstatic -lsimtools -Wl,-Bdynamic

gtc_process: gtc_process.o libsimtools.a
	$(CXX) $< $(LDFLAGS) -o $@ -lm -lstlport -Wl,-Bstatic -lsimtools -Wl,-Bdynamic

gtc_process.o: gtc_process.cpp
	$(CXX) -c -DTEST $(CXXFLAGS) -o $@ $<

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

%.swig.o: %.cpp
	$(CXX) -c -DSWIG $(CXXFLAGS) $(PERL_CC_OPTS) -o $@ $<

%.swig.o: %.cxx
	$(CXX) -c -DSWIG $(CXXFLAGS) $(PERL_CC_OPTS) -o $@ $<

Gtc_wrap.cxx: Gtc.i
	swig -perl -c++ -shadow -Wall Gtc.i

Sim_wrap.cxx: Sim.i
	swig -perl -c++ -shadow -Wall Sim.i

Gtc.so Gtc.pm: Gtc_wrap.swig.o Gtc.swig.o Manifest.swig.o gtc_process.swig.o win2unix.swig.o
	$(CXX) -shared $(PERL_LD_OPTS) -o $@ $^

Sim.so Sim.pm: Sim_wrap.swig.o Sim.swig.o
	$(CXX) -shared $(PERL_LD_OPTS) -o $@ $^

libsimtools.so: Sim.o Gtc.o Manifest.o QC.o json/json_reader.o json/json_writer.o json/json_value.o utilities.o plink_binary.o gtc_process.o win2unix.o
	$(CXX) -shared $(LDFLAGS) -o $@ $^

libsimtools.a: Sim.o Gtc.o Manifest.o QC.o json/json_reader.o json/json_writer.o json/json_value.o utilities.o plink_binary.o gtc_process.o win2unix.o
	$(AR) rcs $@ $^
