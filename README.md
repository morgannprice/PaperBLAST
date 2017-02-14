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
by bin/buildLitDb.pl

Create a tmp/ subdirectory and set it to writable by apache (or
everyone) so that litSearch.cgi can write to ../tmp/

# Dependencies

sqlite3 is required

To identify search terms from MicrobesOnline, the MicrobesOnline code
base must be in ~/Genomics. (This is not required by the CGI scripts.) See
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

To prepare the input files for buildLitDb.pl, you need to download
RefSeq, UniProt, the open-access part of EuropePMC, and EcoCyc, and
you need to run many queries with the API of EuropePMC. The steps are:

mkdir oa # directory for open access part of EuropePMC

mkdir oa.out # directory for intermediate results

wget -O - --no-check-certificate https://europepmc.org/ftp/oa > oa/listing

perl -ane 'print $1."\n" if m/"(PMC[a-zA-Z90-9_.-]+[.]gz)"/' < oa/listing > oa/files

### Find potential locus tags
(for i in `cat oa/files | sed -e s/.xml.gz//`; do echo $i; gunzip -c oa/$i.xml.gz | bin/words.pl > oa.out/$i.words; done) >& words.log

### Find locus tags corresponding to those words, using MicrobesOnline
(cat oa.out/PMC*.words | bin/moIds.pl > PMC.oa.links) >& moIds.log

### Build queries from those locus tags
(bin/oaquery.pl < PMC.oa.links > queryprot.oa) >& queryprot.oa.log &

### Build queries from associations found using an earlier version of PaperBLAST
bin/hitsToTerms2.pl < pubcache > queryprot.pubcache

### Build queries for the most popular genomes (so that locus tags can be found even if the paper is secret)
cat queryprot.pubcache queryprot.oa | sort -u | ~/Genomics/util/tabulate.pl -skip 0 -column 1 > cnt.taxnames

(head -301 cnt.taxnames | tail -300 | nice bin/queryProtByOrg.pl > queryprot.pop) >& queryprot.pop.log &

### download refseq (all the complete*.genomic.gbff.gz files). Note release79 would change.
mkdir refseq

cd refseq

mkdir complete

cd complete

(for i in `grep complete ../release79.files.installed  | cut -f 2 | grep genomic.gbff.gz | sort`; do echo Fetching $i; wget -nv ftp://ftp.ncbi.nih.gov/refseq/release/complete/$i; done) >& fetch.log

### Build queries from RefSeq
cut -f 2 oa.out/PMC*.words | sort -u > oa.words

(for i in `(cd /usr2/people/gtl/data/refseq/release79; ls complete.*.genomic.gbff.gz | sed -e 's/.genomic.gbff.gz//')`; do echo "zcat /usr2/people/gtl/data/refseq/release79/$i.genomic.gbff.gz | bin/findRefSeqQueries.pl oa.words > refseq/$i.loci"; done) > refseq/cmds

nice bin/submitter.pl refseq/cmds >& refseq/cmds.log

(cat refseq/*.loci | bin/convertRefSeqQueries.pl > queryprot.refseq) >& refseq.convert.log

(bin/removeDupQueries.pl queryprot.oa2 queryprot.pop < queryprot.refseq > queryprot.refseqnew) >& queryprot.refseq.duplog

### Run all the queries against EuropePMC

mkdir parts # staging area for temporary files

cat queryprot.oa queryprot.pubcache | sort -u > queryprot.oa2

###### queryEuropePMCBatch.pl submits queries in parallel, throttled to 5/second
###### Running all of these will probably take over a week
bin/queryEuropePMCBatch.pl -dir parts -in queryprot.oa2 -out hits.oa2 >& hits.oa2.log

bin/queryEuropePMCBatch.pl -dir parts -in queryprot.pop -out pop.oa2 >& hits.pop.log

bin/queryEuropePMCBatch.pl -dir parts -in queryprot.refseqnew -out refseq.hits >& refseq.hits.log

### Combine the results, find snippets, and build the database
bin/parseEuropePMCHits.pl -in queryprot.oa2 -hits hits.oa2 -out epmc_oa2 >& epmc_oa2.log

bin/parseEuropePMCHits.pl -in queryprot.pop -hits pop.oa2 -out epmc_pop >& epmc_pop.log

bin/parseEuropePMCHits.pl -in queryprot2.pop -hits pop2.oa2 -out epmc_pop2 >& epmc_pop2.log

bin/parseEuropePMCHits.pl -in queryprot.refseqnew -hits tmp.hits -out epmc_refseq >& epmc_refseq.log

(for i in `cat oa/files | sed -e s/.xml.gz//`; do echo "./buildSnippets.pl -in oa/$i.xml.gz -out oa.out/$i.snippets epmc_oa2.papers epmc_pop.papers epmc_pop2.papers epmc_refseq.papers >& oa.out/$i.snippets.log"; done) > oa.out/snippets.cmds

bin/submitter.pl oa.out/snippets.cmds >& oa.out/snippets.cmds.log

bin/buildLitDb.pl -snippets snippets -dir data -sprot sprot.char.tab epmc_oa2 epmc_pop epmc_pop2 epmc_refseq -ecocyc ecocyc/data

# The Databases

PaperBLAST uses a SQL database (currently implemented with sqlite3)
with information about proteins and links to papers as well as a
protein BLAST database.

See bin/litsearch.sql for the sqlite3 database schema. Note that genes
curated functions in UniProt and from EcoCyc do not link to papers in
the same way that the other genes (from MicrobesOnline or RefSeq)
do. UniProt.comment contains the information about papers for UniProt
genes. EcoCycToPubMed contains the links to papers from EcoCyc. The
other links are in GenePaper. The sqlite3 database must be in
../data/litsearch.db (this path is relative to the CGI scripts).  The
full set of proteins are in ../data/litsearch.faa. The proteins with
unique sequences is in ../data/uniq.faa, which is also a BLASTp
database.

-- Morgan Price
Arkin group
Lawrence Berkeley National Lab
February 2017
