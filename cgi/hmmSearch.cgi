#!/usr/bin/perl -w
#######################################################
## hmmSearch.cgi
##
## Copyright (c) 2018 University of California
##
## Authors: Morgan Price
#######################################################
#
# Optional CGI parameters:
# hmmId -- what hmm to search for
# curated -- show curated hits only

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils;
use pbweb;

my $cgi = CGI->new;
my $hmmId = $cgi->param('hmmId');
my $curatedOnly = $cgi->param('curated');

my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $tmpDir = "../tmp";
my $sqldb = "$base/litsearch.db";
die "No such file: $sqldb" unless -e $sqldb;
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $hmmfile = HmmToFile($hmmId);

if (!defined $hmmfile) {
  # print a form
  my $title = "Family Search (HMMer) vs. Papers";
  print
    header(-charset => 'utf-8'),
    start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -title => $title),
    h2($title);
  print p("Did not find an HMM matching $hmmId") if $hmmId;
  print
    GetMotd(),
    start_form(-name => 'select', '-method' => "GET", -action => 'hmmSearch.cgi'),
    p(small("Examples:",
            a({-href => "hmmSearch.cgi?hmmId=TIGR00001&curated=1",
               -title => "ribosomal protein L35 from TIGRFam"},
              "TIGR00001"),
            a({-href => "hmmSearch.cgi?hmmId=PF05853",
               -title => "&beta;-keto acid cleavage enzyme family from Pfam"},
              "PF05853"))),
    p("HMM Accession:", textfield(-name => "hmmId", -value => '', -size => 12)),
    p(checkbox(-name => 'curated', -checked => 0, -label => "Show hits to curated sequences only")),
    p(submit('Search')),
    end_form,
    a({-href => "litSearch.cgi"}, "Or search by sequence"),
    end_html;
  exit(0);
}
# else have $hmmfile
my $title = "Family Search for $hmmId";
  print
    header(-charset => 'utf-8'),
    start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -title => $title),
    qq{<SCRIPT src="../static/pb.js"></SCRIPT>\n},
    h2($title),
    GetMotd();
autoflush STDOUT 1; # show preliminary results
print "\n";

# See if the search results exist already
my $resultsFile = "$tmpDir/$hmmId.hmmer.dom";
unless (-e $resultsFile
        && NewerThan($resultsFile, $hmmfile)
        && NewerThan($resultsFile,  $blastdb)) {
  my $hmmsearch = "../bin/hmmsearch";
  die "No such executable: $hmmsearch" unless -x $hmmsearch;
  my $tmpResultsFile = "$resultsFile.$$.tmp";
  print p("Running HMMer for $hmmId") . "\n";
  system($hmmsearch, "--cut_tc", "-o", "/dev/null", "--domtblout", $tmpResultsFile, $hmmfile, $blastdb) == 0
    || die "Error running hmmsearch: $!";
  rename($tmpResultsFile, $resultsFile)
    || die "Error renaming to $resultsFile";
}

# Parse the results. For each unique sequence id, a hash
# of hmmAcc, hmmlen, evalue, score, and hits, which is a list of hashes,
# each containing domscore, hmmfrom, hmmto, envfrom, envto [envelope coordinates]
my %hits = (); 

open(IN, "<", $resultsFile) || die "Cannot read $resultsFile";
while(my $line = <IN>) {
  next if $line =~ m/^#/;
  my ($uniqId, undef, $seqlen, undef, $hmmAcc, $hmmlen,
      $evalue, $score, undef,
      undef, undef, undef, undef, $domscore, undef, $hmmfrom, $hmmto, undef, undef,
      $envfrom, $envto) = split / +/, $line;
  die "Cannot parse line $line"
    unless defined $envto
      && $seqlen =~ m/^\d+$/ && $hmmlen =~ m/^\d+$/
      && $hmmfrom =~ m/^\d+$/ && $hmmto =~ m/^\d+$/
      && $envfrom =~ m/^\d+$/ && $envto =~ m/^\d+$/;
  $hits{$uniqId} = { evalue => $evalue, score => $score, hmmAcc => $hmmAcc, hmmlen => $hmmlen, hits => [] }
    unless exists $hits{$uniqId};
  push @{ $hits{$uniqId}{hits} },
    { domscore => $domscore, hmmfrom => $hmmfrom, hmmto => $hmmto, envfrom => $envfrom, envto => $envto };
}
close(IN) || die "Error reading $resultsFile";

if (keys %hits > 0) {
  my $first = (values %hits)[0];
  my $URL = undef;
  $URL = "http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=$hmmId" if $hmmId =~ m/^TIGR\d+$/;
  $URL = "https://pfam.xfam.org/family/$hmmId" if $hmmId =~ m/^PF\d+$/;
  my $acc = $first->{hmmAcc};
  $acc = a({ -href => $URL }, $acc) if defined $URL;
  my @parts = ($acc, "hits", scalar(keys %hits), "sequences in PaperBLAST's database.");
  if ($curatedOnly) {
    push @parts, "Showing hits to curated sequences only.",
      "Or see",
        a({-href => "hmmSearch.cgi?hmmId=$hmmId"}, "all hits");
  } else {
    push @parts, "Showing all hits.",
      "Or show",
      a({-href => "hmmSearch.cgi?hmmId=$hmmId&curated=1"}, "only hits to curated sequences");
  }
  push @parts, "or try", a({-href => "hmmSearch.cgi"}, "another family.");
  print p(@parts);
} else {
  print p("Sorry, no hits for $hmmId. Try",
          a({-href => "hmmmSearch.cgi"}, "another family."));
}
print "\n";

# For speed, sort by overall score
my @uniq_sorted = sort { $hits{$b}{score} <=> $hits{$a}{score} } (keys %hits);

my $maxPapers = 8;
my $nShown = 0;
my $maxShow = 250;
foreach my $uniqId (@uniq_sorted) {
  my @genes = UniqToGenes($dbh, $uniqId);
  # 5 is fitness browser, the lowest curated setting
  @genes = grep { $_->{priority} <= 5 } @genes
    if $curatedOnly;
  next if @genes == 0;
  $nShown++;
  if ($nShown > $maxShow) {
    print p("Additional hits are not shown.");
    last;
    next;
  }
  my $alns = $hits{$uniqId}{hits};
  # minFrom and maxTo are in the hit protein's coordinates
  my $seqlen = $genes[0]{protein_length};
  my $minFrom = $seqlen;
  my $maxTo = 1;
  my $seqCover = 0; # total coverage of query, either summing aligments or total range
  my $maxHMMCoverage = 0;
  foreach my $aln (@$alns) {
    $minFrom = $aln->{envfrom} if $aln->{envfrom} < $minFrom;
    $maxTo = $aln->{envto} if $aln->{envto} > $maxTo;
    my $cov = ($aln->{hmmto} - $aln->{hmmfrom} + 1) / $hits{$uniqId}{hmmlen};
    $maxHMMCoverage = $cov if $cov > $maxHMMCoverage;
    $seqCover += ($aln->{envto} - $aln->{envfrom} + 1) / $seqlen;
  }
  my $seqCover2 = ($maxTo-$minFrom+1) / $seqlen;
  $seqCover = $seqCover2 if $seqCover2 < $seqCover;
  my $alnstring = scalar(@$alns) == 1 ? "Aligns from $minFrom:$maxTo/$seqlen"
    :  scalar(@$alns) . " alignments in $minFrom:$maxTo/$seqlen";
  my $coverage_html = br()
    . join(" ",
           a({ -href => "hmmAlign.cgi?hmmId=$hmmId&acc=$uniqId" }, $alnstring),
           "(" . sprintf("%.1f%%", $seqCover*100) . "),",
           scalar(@$alns) == 1 ? "covers" : "covering up to",
           sprintf("%.1f%%", $maxHMMCoverage*100.0),
           "of $hmmId,",
           $hits{$uniqId}{score}, "bits");
  print GenesToHtml($dbh, $uniqId, \@genes, $coverage_html, $maxPapers);
  print "\n";
}

print end_html;

