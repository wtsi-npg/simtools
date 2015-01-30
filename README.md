Simtools
========

Program to create and process SIM files for the <a href="https://github.com/wtsi-npg/genotyping">WTSI genotyping pipeline</a>

Executables
-----------

* simtools
  * Create SIM from GTC
  * View SIM contents
  * Generate input for Illuminus or GenoSNP from SIM
  * Compute QC metrics
  * Create Illumina Final Call Report (FCR) from GTC
* normalize_manifest
  * Normalize a .bpm.csv manifest to Illumina TOP strand
* sim
  * Demonstration program to parse SIM files
* g2i
  * Older program to convert GTC directly to Illuminus input
  * No longer supported
* gtc
  * Process GTC files for g2i
  * No longer supported

Run any executable with `--help` for more information.

Installation and testing
------------------------

* To compile: `make all`
* To run tests: `make test`. Requires `cxxtestgen` on the PATH. In the WTSI environment, can use `module load cxxtest/$(VERSION)`.
* To install simtools and normalize_manifest: `make install DIR=$(DESTINATION)` will install executables to the `bin` subdirectory of `$DESTINATION`.
* To install g2i to /software/varinf: `make install_g2i`




