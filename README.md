PaperBLAST is a tool to find papers about homologs of a protein of interest. For an example see

http://fit.genomics.lbl.gov/cgi-bin/papers/cgi/litSearch.cgi?vimss=14484

# System Requirements

These scripts should work on any Linux system, and would probably work
on other Unix or MacOS as well. All of the code is written in perl.

# Installation

Please put the blast executables (formatdb,
blastall, and fastacmd) in the bin/blast/ subdirectory.

Please put the usearch executable in the bin/ subdirectory. I have only
tested usearch version 8.0.

The cgi/ subdirectory contains the common gateway interface script,
litSearch.cgi, which is the interface for using PaperBLAST.

For litSearch.cgi to work, the data/ subdirectory must include the
sqlite3 database and the protein BLAST database. These are both built
by bin/run_final.pl (which invokes bin/buildLitDb.pl).

Create a tmp/ subdirectory and set it to writable by apache (or
everyone) so that litSearch.cgi can write to ../tmp/ (relative to the
directory that it is invoked from, which is usually the cgi directory)

# Dependencies

sqlite3 is required

To identify search terms from MicrobesOnline, the MicrobesOnline code
base must be in ~/Genomics. (This is required for building the
database; it is not required by the CGI scripts.) See
http://www.microbesonline.org/programmers.html#source

The swissknife perl library must be in the SWISS subdirectory (so that
the *.pm files are in SWISS/lib/SWISS/). It is available at

http://swissknife.sourceforge.net/docs/
https://sourceforge.net/projects/swissknife/files/latest/download

Finally, PaperBLAST uses these standard Perl libraries:

Getopt::Long
FindBin
File::Which
IO::Handle
XML::LibXML
JSON
POSIX
Time::HiRes
LWP::Simple
URI::Escape
CGI
DBI

# Building the Database

The database can be built with the following series of scripts. As of
March 2017, the process requires approximately 324 GB of free disk
space and takes around a week to run. The majority of the time is for
querying EuropePMC.

Download information from RefSeq, UniProt, EcoCyc, EuropePMC, and
PubMed:

	mkdir indata
	nice bin/download.pl -dir indata >& download.log

Choose which queries to run against EuropePMC:

	mkdir work
	nice bin/run_terms.pl -in indata -work work >& run_terms.log

Run all the queries against EuropePMC:

	nice bin/run_search.pl -in indata -work work >& run_search.log

Extract snippets:

	nice bin/run_snippets.pl -in indata -work work >& run_snippets.log

Before using run_snippets, you should set up a key for the elsevier API
and a token for the CrossRef API. These should go in the cache/ directory, in files named elsevier key and crossref_token.

The cache/ directory is also where the articles that are downloaded
from the Elsevier API or from Crossref are stored. They are stored as
cache_nnnn.pdf or .xml, where the number is the pubmed id. If you use
other tools such as pubMunch2 to obtain articles, you need to rename
the PDFs to this format and put them in the cache/ directory.

Build the database:

	nice bin/run_final.pl -in indata -work work >& run_final.log

This will create the directory work/data, which will include a BLAST
database of unique protein sequences and a sqlite3 relational
database. See bin/litsearch.sql for the sqlite3 database schema. These
contents can be copied or symbolically linked to the data/ directory,
which is where the CGI scripts expect the database to be. (More
precisely, the data should be in ../data/ relative to the directory
that the CGI script is invoked from.)

-- Morgan Price, Arkin group, Lawrence Berkeley National Lab, February 2017
