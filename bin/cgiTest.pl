#!/usr/bin/perl -w
use strict;
use DBI;
use FindBin qw{$RealBin};
use Getopt::Long;
sub Test($$);

my $verbose;
die "Usage: cgiTest.pl [ -verbose ]\n"
  unless GetOptions('verbose' => \$verbose)
  && @ARGV == 0;


# Redirect stdout to /dev/null
open(STDOUT, ">", "/dev/null") || die "Cannot redirect STDOUT";

# cd to the cgi directory
chdir("$RealBin/../cgi") || die "Cannot find cgi directory";

# Tests for PaperBLAST
Test("litSearch.cgi", "query=VIMSS3615187");
Test("litSearch.cgi", "query=AO353_07705");
Test("litSearch.cgi", "more=WP_000168720.1");
Test("litSearch.cgi", "more=VIMSS6582571");
Test("litSearch.cgi", "more=Q7DFN2");

# and for family search vs. papers
Test("hmmSearch.cgi", "hmmId=TIGR00001&curated=1");
Test("hmmSearch.cgi", "hmmId=TIGR00001");
Test("hmmSearch.cgi", "hmmId=PF05853&curated=1");

# Tests for curated blast for genomes
Test("genomeSearch.cgi", "gdb=NCBI&gid=GCF_000236665.1&query=perchlorate");

# Tests for GapMind and curated clusters

# A few locus- or genome-specific tests
my $orgId = "FitnessBrowser__Keio";
Test("gapView.cgi", "set=carbon&orgs=orgsFit&orgId=$orgId&locusId=18312");
Test("gapView.cgi", "set=carbon&orgs=orgsFit&orgId=$orgId&locusId=18312&path=citrate&step=fecD");
Test("gapView.cgi", "orgs=orgsFit&set=carbon&orgId=FitnessBrowser__Keio&path=citrate&step=fecD&locusId=18312");
Test("gapView.cgi", "orgs=orgsFit&set=carbon&orgId=FitnessBrowser__Keio&findgene=dicitrate");

# Test every set, path, step
foreach my $set (qw{aa carbon}) {
  Test("gapView.cgi", "set=$set");
  Test("gapView.cgi", "set=$set&orgs=orgsFit");
  Test("gapView.cgi", "set=$set&orgs=orgsFit&gaps=1");
  Test("gapView.cgi", "set=$set&orgs=orgsFit&orgId=$orgId");

  my $stepsDir = "$RealBin/../tmp/path.$set";
  my $dbhS = DBI->connect("dbi:SQLite:dbname=${stepsDir}/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my $pathIds = $dbhS->selectcol_arrayref(qq{SELECT pathwayId FROM Pathway WHERE pathwayId <> "all"});
  foreach my $pathId (@$pathIds) {
    Test("gapView.cgi", "set=$set&orgs=orgsFit&path=$pathId");
    Test("gapView.cgi", "set=$set&orgs=orgsFit&path=$pathId&showdef=1");
    Test("gapView.cgi", "set=$set&orgs=orgsFit&path=$pathId&showdef=1&literal=1");
    Test("gapView.cgi", "set=$set&orgs=orgsFit&path=$pathId&orgId=$orgId");
    Test("curatedClusters.cgi", "set=$set&path=$pathId");
    my $stepIds = $dbhS->selectcol_arrayref(qq{SELECT stepId FROM Step WHERE pathwayId = ?},
                                            {}, $pathId);
    foreach my $stepId (@$stepIds) {
      Test("gapView.cgi", "set=$set&orgs=orgsFit&path=$pathId&step=$stepId&orgId=$orgId");
      Test("curatedClusters.cgi", "set=$set&path=$pathId&step=$stepId");
      Test("curatedClusters.cgi", "set=$set&path=$pathId&step=$stepId&close=1");
    }
  }
  $dbhS->disconnect || die $DBI::errstr;
}
print STDERR "cgiTest.pl done\n";

sub Test($$) {
  my ($cgi, $arg) = @_;
  my @cmd = ("./$cgi", $arg);
  die "No such cgi: $cmd[0]" unless -x $cmd[0];
  print STDERR "Testing: $cgi $arg\n"
    if defined $verbose;
  my $tmpFile = "/tmp/cgiTest.$$.err";
  # redirect STDERR
  open(my $STDERROLD, ">&", STDERR) || die "Cannot save STDERR";
  open (STDERR, ">", $tmpFile) || die "Cannot redirect STDERR";
  my $status = system(@cmd);
  # restore STDERR
  open(STDERR, ">&", $STDERROLD);
  die "Failed running: $cgi $arg\nSee stderr in $tmpFile\n"
    if $status != 0;
  # warn if stderr is not empty
  print STDERR "Nonempty stderr from $cgi $arg\n"
    if -s $tmpFile;
  unlink($tmpFile);
}
