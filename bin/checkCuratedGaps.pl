#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use Steps qw{ReadOrgTable};
use pbutils qw{ReadTable NewerThan};
use DBI;

my $stepBase = "$RealBin/../gaps";
my $resultsDir = "$RealBin/../tmp";
my $set = "aa";
my ($stepDir, $queryDir);

my $usage = <<END
checkCuratedGaps.pl -org orgs35 ... orgsNif

Given a pathway set and one or more org directories with computed
results, validate that the curated/known gaps are consistent with those
results.

Optional arguments:
-set $set -- which pathway set
-stepDir $stepBase/$set -- the directory with steps.db
  and the $set.curated.gaps.tsv table
-queryDir $RealBin/../tmp/path.$set -- the directory with the
  *.query files
-results $resultsDir -- the directory that contains the
  org directories.
  The analysis results should be in
  $resultsDir/orgDir/$set.sum.*
END
;

my @orgDir;
die $usage
  unless GetOptions('set=s' => \$set,
                    'org=s{,}' => \@orgDir,
                    'stepDir=s' => \$stepDir,
                    'results=s' => \$resultsDir)
  && @ARGV == 0;
die "No organism directories specified\n" unless @orgDir > 0;
$stepDir = "$stepBase/$set" if !defined $stepDir;
$queryDir = "$RealBin/../tmp/path.$set" if !defined $queryDir;
foreach my $dir ($stepDir, $queryDir, $resultsDir, map { "$resultsDir/$_" } @orgDir) {
  die "No such directory: $dir\n"
    unless -d $dir;
}

my $dbhS = DBI->connect("dbi:SQLite:dbname=${queryDir}/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $pathways = $dbhS->selectall_arrayref("SELECT * from Pathway",
                                         { Slice => {} });
my %pathways = map { $_->{pathwayId} => $_ } @$pathways;
my $steps = $dbhS->selectall_arrayref("SELECT * from Step",
                                         { Slice => {} });
my %steps;
foreach my $step (@$steps) {
  $steps{ $step->{pathwayId} }{ $step->{stepId} } = $step;
}

my @cgapHeader = qw{genomeName gdb gid pathway step class comment};
my $cgapFile = "$stepDir/$set.curated.gaps.tsv";
my @cgaps;
if (-e $cgapFile) {
  @cgaps = ReadTable($cgapFile,\@cgapHeader);
  my @noclass = grep { $_->{class} eq "" } @cgaps;
  @cgaps = grep { $_->{class} ne "" } @cgaps;
  print STDERR "Read " . scalar(@cgaps) . " non-empty curated gaps from\n$cgapFile\n";
  print STDERR "Skipped " . scalar(@noclass) . " entries with empty class\n"
    if @noclass > 0;
} else {
  print STDERR "No curated gaps for $set\n";
}

my @kgapHeader = qw{genomeName gdb gid pathway step};
my $kgapFile = "$stepDir/$set.known.gaps.tsv";
my @kgaps;
if (-e $kgapFile) {
  @kgaps = ReadTable($kgapFile, \@kgapHeader);
  print STDERR "Read " . scalar(@kgaps) . " known gaps from\n$kgapFile\n";
} else {
  print STDERR "No known gaps for $set\n";
}

# gdb => gid => pathway => step => curated gap object
my %cgap = ();
foreach my $cgap (@cgaps) {
  print STDERR "Duplicate curated gaps for $cgap->{gdb} $cgap->{gid} $cgap->{pathway} $cgap->{step}\n"
    if exists $cgap{ $cgap->{gdb} }{ $cgap->{gid} }{ $cgap->{pathway} }{ $cgap->{step} };
  $cgap{ $cgap->{gdb} }{ $cgap->{gid} }{ $cgap->{pathway} }{ $cgap->{step} } = $cgap;
}
foreach my $kgap (@kgaps) {
  $cgap{ $kgap->{gdb} }{ $kgap->{gid} }{ $kgap->{pathway} }{ $kgap->{step} } = $kgap
    if !exists $cgap{ $kgap->{gdb} }{ $kgap->{gid} }{ $kgap->{pathway} }{ $kgap->{step} };
}

my %orgSeen = (); # gdb => gid => genome name

my $dateFile = "$queryDir/date";
die "No such file: $dateFile -- queries not built\n" unless -e $dateFile;
foreach my $org (@orgDir) {
  my $doneFile = "$resultsDir/$org/$set.sum.db";
  unless (NewerThan($doneFile, $dateFile)) {
    print STDERR "Skipping the organisms in $org -- not up to date\n";
    next;
  }
  my @orgs = ReadOrgTable("$resultsDir/$org/orgs.org");
  print STDERR "Checking $resultsDir/$org with " . scalar(@orgs) . " organisms\n";
  my @sumSteps = ReadTable("$resultsDir/$org/$set.sum.steps",
                           qw{orgId gdb gid pathway step onBestPath score locusId sysName
                              score2 locusId2 sysName2});
  my %sumSteps; # gdb => gid => list
  foreach my $row (@sumSteps) {
    push @{ $sumSteps{ $row->{gdb} }{ $row->{gid} } }, $row;
  }
  my @sumRules = ReadTable("$resultsDir/$org/$set.sum.rules",
                           qw{orgId gdb gid pathway rule score nHi nMed nLo expandedPath});
  my %sumRules = (); # gdb => gid => list
  foreach my $row (@sumRules) {
    push @{ $sumRules{ $row->{gdb} }{ $row->{gid} } }, $row;
  }
  foreach my $org (@orgs) {
    my $gdb = $org->{gdb};
    my $gid = $org->{gid};
    print STDERR "Inconsistent genome names for $gdb $gid in $org vs. earlier directories\n"
      . " $org->{genomeName} vs. $orgSeen{$gdb}{$gid}\n"
        if exists $orgSeen{$gdb}{$gid}
          && $orgSeen{$gdb}{$gid} ne $org->{genomeName};
    $orgSeen{$gdb}{$gid} = $org->{genomeName};
    foreach my $stepObj (@{ $sumSteps{$gdb}{$gid} }) {
      my $pathwayId = $stepObj->{pathway};
      my $step = $stepObj->{step};
      print STDERR "Curated/known gap $pathwayId $step for $gdb $gid ($org->{genomeName}) has score=2\n"
        if $stepObj->{score} eq "2"
          && exists $cgap{$gdb}{$gid}{$pathwayId}{$step};
      print STDERR "Curated gap $pathwayId $step for $gdb $gid ($org->{genomeName}) is not on the best path\n"
        if ! $stepObj->{onBestPath} && exists $cgap{$gdb}{$gid}{$pathwayId}{$step};
    }
    foreach my $ruleObj (@{ $sumRules{$gdb}{$gid} }) {
      my $pathwayId = $ruleObj->{pathway};
      next unless $ruleObj->{rule} eq "all";
      print STDERR "Curated/known gap $pathwayId (no step) for $gdb $gid ($org->{genomeName}) has a high-confidence path\n"
        if $ruleObj->{nLo} == 0 && $ruleObj->{nMed} == 0
          && exists $cgap{$gdb}{$gid}{$pathwayId}{""};
    }
  }
}

# Check each row of cgaps
foreach my $cgap (@cgaps) {
  my $gdb = $cgap->{gdb};
  my $gid = $cgap->{gid};
  if (!exists $orgSeen{$gdb}{$gid}) {
    print STDERR "Organism $gdb $gid ($cgap->{genomeName}) not checked\n";
    $orgSeen{$gdb}{$gid} = "";
  } elsif ($orgSeen{$gdb}{$gid} ne $cgap->{genomeName} && $orgSeen{$gdb}{$gid} ne "") {
    print STDERR "Inconsistent genome name for $gdb $gid\n"
      . "  Curated: $cgap->{genomeName} Input directory: $orgSeen{$gdb}{$gid}\n";
  }
  print STDERR "No comment for $cgap->{pathway} $cgap->{step} $cgap->{genomeName} (gdb $gdb gid $gid)\n"
    unless $cgap->{comment} =~ m/[a-zA-Z]/;
  my $pathwayId = $cgap->{pathway};
  my $step = $cgap->{step};
  if (!exists $pathways{$pathwayId}) {
    print STDERR "Unknown pathway $pathwayId in curated table\n";
  } elsif ($step ne "" && !exists $steps{$pathwayId}{$step}) {
    print STDERR "Unknown step $step for pathway $pathwayId in curated table\n";
  }
}
