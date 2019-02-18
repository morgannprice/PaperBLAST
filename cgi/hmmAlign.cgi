#!/usr/bin/perl -w
#######################################################
## hmmAlign.cgi
##
## Copyright (c) 2018 University of California
##
## Authors: Morgan Price
#######################################################
#
# Required CGI parameters:
# hmmId -- the model to use
# acc -- the accession to compare the model to

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use lib "../lib";
use pbweb; # for FetchFasta()

my $cgi = CGI->new;
my $hmmId = $cgi->param('hmmId');
my $acc = $cgi->param('acc');
die "Must specify hmmId and acc" unless $hmmId && $acc;

my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $tmpDir = "../tmp";
my $sqldb = "$base/litsearch.db";
die "No such file: $sqldb" unless -e $sqldb;
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $hmmfile = HmmToFile($hmmId);
die "No hmm file for $hmmId" unless defined $hmmfile;
my $isUploaded = $hmmId =~ m/^hex[.]/;

my $hmmName = `egrep '^NAME' $hmmfile`;
$hmmName =~ s/[\r\n].*//;
$hmmName =~ s/^NAME +//;

my $showId = $isUploaded ? "uploaded HMM " . escapeHTML($hmmName) : $hmmId;
# Quote accession to avoid cross-site scripting attacks
my $title = "Align " . escapeHTML($acc) . " to $showId";
$title .= " (" . escapeHTML($hmmName) . ")" unless $isUploaded || $hmmName eq $hmmId;
start_page('title' => $title);

my $db = "$base/uniq.faa";
die "No such file: $db" unless -e $db;

my $protseq = FetchFasta($dbh, $db, $acc);
print p($acc, "has", length($protseq), "amino acids");

my $tmppre = "$tmpDir/hmmAlign.$$";
my $faafile = "$tmppre.faa";
open(FAA, ">", $faafile) || die "Cannot write to $faafile";
print FAA ">$acc\n", $protseq;
close(FAA) || die "Error writing to $faafile";
my $hmmsearch = "../bin/hmmsearch";
die "No such executable: $hmmsearch" unless -x $hmmsearch;
my $alnfile = "$tmppre.aln";
system("$hmmsearch", "-E", 0.01, "-o", $alnfile, $hmmfile, $faafile) == 0
  || die "Error running hmmsearch: $!";
unlink($faafile);
my @lines = ();
open(ALN, "<", $alnfile) || die "Cannot read $alnfile";
while(my $line = <ALN>) {
  push @lines, $line;
}
close(ALN) || die "Error reading $alnfile";
unlink($alnfile);
my @out = ();
my $writing = 0;
foreach my $line (@lines) {
  $writing = 1  if $line =~ m/^Query:/;
  $writing = 0 if $line =~ m/^Internal pipeline/;
  push @out, $line if $writing;
}


print pre(escapeHTML(join("",@out)));
my $newline = "%0A";
my $CDDURL = "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${acc}$newline$protseq";
print p("Or compare $acc to",
        a({ -href => $CDDURL, -title => 'Conserved Domains Database' }, "CDD"),
        "or",
        a({ -href => "litSearch.cgi?query=$acc" }, "PaperBLAST"));
print end_html;
