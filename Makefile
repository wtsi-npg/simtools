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

DIR = INSTALL_DIRECTORY # placeholder for main installation directory
BIN=$(DIR)/bin
# default directories for g2i 
INSTALLLIB=/software/varinf/lib 
INSTALLBIN=/software/varinf/bin
# other useful paths
STLPORT_INC=/software/solexa/pkg/STLport/current/stlport
STLPORT_LIB=/software/solexa/pkg/STLport/current/build/lib/obj/gcc/so
#STLPORT_INC=/software/varinf/lib/STLport/include/stlport
#STLPORT_LIB=/software/varinf/lib/STLport/lib
PERL_CORE=/usr/lib/perl/5.8.8/CORE
PERL_CORE=/software/perl-5.8.8/lib/5.8.8/x86_64-linux-thread-multi/CORE


TARGETS=libplinkbin.so gtc g2i gtc_process sim simtools normalize_manifest
PERL_TARGETS=Gtc.so Sim.so
LIBS=Gtc.o win2unix.o Sim.o

# To compile with debug information added, invoke as (say): 
# make DEBUG='y'
# to build all targets with debug info. 
# For just one target, say:
# make DEBUG='y' simtools

CC=/usr/bin/g++

# do NOT use -ffast-math, as it causes errors in infinity/NaN handling
ifeq ($(DEBUG),y)
	CFLAGS=-g -Wall -fPIC -O0 -I$(STLPORT_INC) -std=c++0x
else
	CFLAGS=-Wall -fPIC -O3 -I$(STLPORT_INC) -std=c++0x
endif
# Set runpath instead of relying on LD_LIBRARY_PATH
LDFLAGS=-Wl,-rpath -Wl,$(STLPORT_LIB) -L$(STLPORT_LIB) -lstlport -lm 
CXXFLAGS=-Wno-deprecated -I/software/gapi/pkg/cxxtest/4.2.1/

usage:
	@echo -e "Usage: make install DIR=<destination directory>\nOther targets: make all, make test, make install_g2i"


clean:
	rm -f *.o Gtc_wrap.cxx Gtc.pm Sim_wrap.cxx Sim.pm runner.cpp runner $(TARGETS)

test: Sim.o Gtc.o Manifest.o QC.o win2unix.o json/json_reader.o json/json_writer.o json/json_value.o commands.o runner.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(CXXFLAGS) -o runner $^ -lstlport
	./runner # run "./runner -v" to print trace information

runner.o: all
	cxxtestgen --error-printer -o runner.cpp test_simtools.h
	$(CC) -c $(CFLAGS) $(LDFLAGS) $(CXXFLAGS)  -o $@ runner.cpp

install: all
	install -d $(BIN) 
	install g2i gtc normalize_manifest sim simtools $(BIN)
	@echo -e "Simtools successfully installed."

install_g2i: all
	cp Gtc.pm $(INSTALLLIB)
	cp Gtc.so $(INSTALLLIB)
	cp Sim.pm $(INSTALLLIB)
	cp Sim.so $(INSTALLLIB)
	cp g2i $(INSTALLBIN)

all: $(TARGETS)

perl: $(PERL_TARGETS)

gtc: gtc.o Gtc.o Manifest.o
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

normalize_manifest: normalize_manifest.o Manifest.o
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

sim: sim.o Sim.o
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

simtools: simtools.o commands.o Sim.o Gtc.o Manifest.o QC.o json/json_reader.o json/json_writer.o json/json_value.o
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

commands.o: commands.cpp
	$(CC) -c  $(CFLAGS) $(CPPFLAGS) -o $@ $<

manifest: manifest.o Manifest.o
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

g2i: g2i.o Gtc.o Manifest.o win2unix.o Sim.o json/json_reader.o json/json_writer.o json/json_value.o utilities.o plink_binary.o
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

gtc_process: Gtc.o Manifest.o gtc_process.o 
	$(CC) $(LDFLAGS) -o $@ $^ -lstlport

gtc_process.o: gtc_process.cpp
	$(CC) -c -DTEST $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

Gtc.so: Gtc.cpp Gtc.h Manifest.cpp gtc_process.cpp win2unix.cpp Gtc.i
	swig -perl -c++ -shadow -Wall Gtc.i
	$(CC) -c -DSWIG -I$(PERL_CORE) -fPIC Gtc.cpp Gtc_wrap.cxx Manifest.cpp gtc_process.cpp win2unix.cpp `perl -MExtUtils::Embed -e ccopt`
	$(CC) -shared Gtc.o Gtc_wrap.o Manifest.o gtc_process.o win2unix.o -o Gtc.so
	perl -e 'use Gtc; $$g=new Gtc::Gtc(); $$g->open("/nfs/new_illumina_geno04/call/20090407/4439467315_R01C01.gtc",$$Gtc::Gtc::BASECALLS);$$m=new Gtc::Manifest();'
	rm Gtc.o Manifest.o gtc_process.o win2unix.o

Sim.so: Sim.cpp Sim.h Sim.i
	swig -perl -c++ -shadow -Wall Sim.i
	$(CC) -c -DSWIG -I$(PERL_CORE) -fPIC Sim.cpp Sim_wrap.cxx `perl -MExtUtils::Embed -e ccopt`
	$(CC) -shared Sim.o Sim_wrap.o -o Sim.so
	rm Sim.o

libplinkbin.so: utilities.o plink_binary.o
	$(CXX) -shared utilities.o plink_binary.o -o $@

Gtc.o: Gtc.cpp Gtc.h
Sim.o: Sim.cpp Sim.h
Manifest.o: Manifest.cpp Manifest.h
gtc.o: Gtc.h Manifest.h
sim.o: Sim.h
simtools.o: Sim.h
g2i.o: Gtc.h Manifest.h plink_binary.h
win2unix.o: win2unix.cpp win2unix.h
plink_binary.o: plink_binary.cpp plink_binary.h

