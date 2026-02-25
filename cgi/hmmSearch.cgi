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
# hmmfile -- an uploaded HMM
# byeval -- ignore trusted cutoff and use (per-domain) E < 0.01 instead

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils;
use pbweb;
use Digest::MD5;

my $maxMB = 10; # maximum size of uploaded HMM
$CGI::POST_MAX = $maxMB*1024*1024;

my $cgi = CGI->new;
my $hmmId = $cgi->param('hmmId');
my $curatedOnly = $cgi->param('curated') ? 1 : 0;

my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $tmpDir = "../tmp";
my $sqldb = "$base/litsearch.db";
die "No such file: $sqldb" unless -e $sqldb;
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $up = $cgi->upload('hmmfile');
if ($up) {
  # Check the uploaded file, save it, set up an MD5-based hmmId, and redirect
  my $fh = $up->handle;
  my @lines = <$fh>;
  my @errors = ();
  push @errors, "Uploaded file does not start with HMMER3"
    unless ($lines[0] =~ m/^HMMER3/i);
  push @errors, "Uploaded file does not end with a // line"
    unless (@lines > 0 && $lines[-1] =~ m!^//\r?\n?!);
  my @endlines = grep m!^//!, @lines;
  push @errors, "Uploaded file has more than one record-ending // line"
    unless @endlines == 1;
  my @namelines = grep m/^NAME /, @lines;
  push @errors, "Uploaded file has no NAME field"
    unless @namelines > 0;
  push @errors, "Uploaded file has more than one NAME"
    unless @namelines <= 1;
  if (@errors > 0) {
    print header,
      start_html(-title => "HMM Upload failed"),
      h2("Not a valid HMM file"),
      p(join(". ", @errors)),
      a({ -href => "hmmSearch.cgi"}, "Try another search");
    finish_page();
  }
  my $name = $namelines[0];
  chomp $name;
  $name =~ s/^NAME +//;

  my $hex = Digest::MD5::md5_hex(@lines);
  my $hmmId = "hex.$hex";
  my $file = "../tmp/$hmmId.hmm";
  unless (-e $file) {
    open(my $sfh, ">", $file) || die "Cannot write to $file";
    print $sfh @lines;
    close($sfh) || die "Error writing to $file";
  }
  print $cgi->redirect("hmmSearch.cgi?hmmId=$hmmId");
  exit(0);
}

$hmmId =~ s/ +$// if defined $hmmId;
$hmmId =~ s/^ +// if defined $hmmId;
my $hmmfile = HmmToFile($hmmId);
my $isUploaded = defined $hmmfile && $hmmId =~ m/^hex[.]/;

if (!defined $hmmfile) {
  # print a form
  my $title = "Family Search (HMMer) vs. Papers";
  start_page('title' => $title);

  if (defined $hmmId && $hmmId ne "") {
    if ($hmmId =~ m/^[a-zA-Z0-9._-]+$/) {
      print p("Did not find an HMM matching $hmmId");
    } else {
      print p("Invalid characters in hmm id");
    }
  }
  print
    GetMotd(),
    start_form(-name => 'select', '-method' => "GET", -action => 'hmmSearch.cgi'),
    p(small("Examples:",
            a({-href => "hmmSearch.cgi?hmmId=TIGR00001&curated=1",
               -title => "ribosomal protein L35 from TIGRFam"},
              "TIGR00001"),
            a({-href => "hmmSearch.cgi?hmmId=PF05853",
               -title => "beta-keto acid cleavage enzyme family from Pfam"},
              "PF05853"),
            a({-href => "hmmSearch.cgi?hmmId=NF001994",
               -title => "FliE (alpha-Proteobacteria type) from NCBI"},
              "NF001994"))),
    p("HMM Accession:", textfield(-name => "hmmId", -value => '', -size => 12)),
    p(checkbox(-name => 'curated', -checked => 0, -label => "Show hits to curated sequences only")),
    p(checkbox(-name => 'byeval', -checked => 0, -label => "Use E < 0.01 instead of trusted cutoff")),
    p(submit('Search')),
    end_form,
    start_form(-name => 'upload', -method => 'POST', -action => 'hmmSearch.cgi'),
    p("Or upload an HMM:",
      filefield(-name => 'hmmfile', -size => 50),
      submit('Go')),
    end_form,
    a({-href => "litSearch.cgi"}, "Or search by sequence");
  finish_page();
}
# else have $hmmfile

my $hmmName = `egrep '^NAME' $hmmfile`;
$hmmName =~ s/[\r\n].*//;
$hmmName =~ s/^NAME +//;
my $showId = $isUploaded ? "uploaded HMM " . escapeHTML($hmmName) : $hmmId;
my $title = "Family Search for $showId";
$title .= " (" . escapeHTML($hmmName) . ")" unless $isUploaded || $hmmName eq $hmmId;

start_page('title' => $title);
print GetMotd();
autoflush STDOUT 1; # show preliminary results
print "\n";

my $byE = $cgi->param('byeval') ? 1 : 0; # ignore trusted cutoff if set

my @tc;
open(my $fh, "<", $hmmfile) || die "Cannot read $hmmfile";
while(my $line = <$fh>) {
  my @parts = split /\s+/, $line;
  if (@parts == 3 && $parts[0] eq "TC") {
    shift @parts;
    @tc = @parts;
    last;
  }
}
close($fh) || die "Error reading $hmmfile";

$byE = 1 if @tc == 0; # Always use by E if no trusted cutoff
my $resultsFile = "$tmpDir/$hmmId.hmmer.dom";
$resultsFile .= ".byE" if $byE;

# See if the search results exist already
unless (-e $resultsFile
        && NewerThan($resultsFile, $hmmfile)
        && NewerThan($resultsFile,  $blastdb)) {
  my $hmmsearch = "../bin/hmmsearch";
  die "No such executable: $hmmsearch" unless -x $hmmsearch;
  my $tmpResultsFile = "$resultsFile.$$.tmp";
  print p("Running HMMer for $showId") . "\n";
  # If it has no trusted cutoff line, try running it with -E 0.001
  my @cmd;
  if ($byE) {
    print p("Searching for hits with E &le; 0.01") . "\n";
     @cmd = ($hmmsearch, "-E", 0.01, "-o", "/dev/null", "--domtblout", $tmpResultsFile, $hmmfile, $blastdb);
  } else {
    @cmd = ($hmmsearch, "--cut_tc", "-o", "/dev/null", "--domtblout", $tmpResultsFile, $hmmfile, $blastdb);
  }
  system(@cmd) == 0
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

my $URL = pbweb::HMMToURL($hmmId);
my $hmmLink = $isUploaded ? "Uploaded HMM" : $hmmId;
$hmmLink = a({-href => $URL}, $hmmLink) if $URL;
if (keys %hits > 0) {
  my $first = (values %hits)[0];
  my $evalueSection;
  if ($byE) {
    $evalueSection = "with E &le; 0.01";
    if (@tc > 0) {
      $evalueSection .= " (or use the "
        . a({-href => "hmmSearch.cgi?hmmId=$hmmId&curated=${curatedOnly}"}, "trusted cutoff")
        . " of " . escapeHTML( join(", ", @tc)  ) . " bits instead)";
    }
    $evalueSection .= ".";
  } else {
    $evalueSection = "above the trusted cutoff of "
      . escapeHTML( join(",", @tc) )
      . " bits (or use "
      . a({-href => "hmmSearch.cgi?hmmId=$hmmId&curated=${curatedOnly}&byeval=1"}, "E < 0.01")
      . " instead).";
  }
  my @parts = ($hmmLink, "hits",
               scalar(keys %hits), "sequences in PaperBLAST's database",
               $evalueSection);
  if ($curatedOnly) {
    push @parts, "<BR>Showing hits to curated sequences only,",
      "or see",
        a({-href => "hmmSearch.cgi?hmmId=$hmmId&byeval=$byE"}, "all hits with papers") . ",";
  } else {
    push @parts, "<BR>Showing all hits,",
      "or show",
      a({-href => "hmmSearch.cgi?hmmId=$hmmId&curated=1&byeval=$byE"}, "hits to curated sequences only") . ",";
  }
  push @parts, "or try", a({-href => "hmmSearch.cgi"}, "another search.");
  print p(@parts);
} else {
  print p("Sorry, no hits for $hmmLink.",
          "Try",
          a({-href => "hmmSearch.cgi"},
            "another family."));
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
  my $alnstring = scalar(@$alns) == 1 ? "Aligns to $minFrom:$maxTo / $seqlen"
    :  scalar(@$alns) . " alignments in $minFrom:$maxTo / $seqlen";
  my $coverage_html = span( { -style => "font-family: sans-serif; font-size: smaller;" },
                            a({ -href => "hmmAlign.cgi?hmmId=$hmmId&acc=$uniqId" }, $alnstring),
                            "(" . sprintf("%.1f%%", $seqCover*100) . "),",
                            scalar(@$alns) == 1 ? "covers" : "covering up to",
                            sprintf("%.1f%%", $maxHMMCoverage*100.0),
                            "of $showId,",
                            a({-title => "per-sequence score"},
                              $hits{$uniqId}{score}, "bits"));
  print GenesToHtml($dbh, $uniqId, \@genes, $coverage_html, $maxPapers);
  print "\n";
}

print p("Or search for genetic data about $hmmId in the",
        a({ -href => "http://fit.genomics.lbl.gov/cgi-bin/myFitShow.cgi?gene=$hmmId" },
          "Fitness Browser"))
  if ! $up && $hmmId && $hmmId =~ m/^(PF|TIGR)\d+$/;

finish_page();
