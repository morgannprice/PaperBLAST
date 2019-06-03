#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use Steps qw{ReadSteps ReadOrgTable ReadReqs CheckReqs};
use pbutils qw{ReadTable NewerThan};

my $stepBase = "$Bin/../gaps";
my $resultsDir = "$Bin/../tmp";
my $set = "aa";
my ($stepDir, $queryDir);

my $usage = <<END
checkCuratedGaps.pl -org orgsNif > warnings.tab

Given a pathway set and an org directory with computed results,
check against the dependency requirements and output a table of warnings.
END
;

my $orgDir;
die $usage
  unless GetOptions('set=s' => \$set,
                    'org=s' => \$orgDir,
                    'stepDir=s' => \$stepDir,
                    'results=s' => \$resultsDir)
  && @ARGV == 0;
die "No organism directory specified\n" unless defined $orgDir;
$stepDir = "$stepBase/$set" if !defined $stepDir;
$queryDir = "$Bin/../tmp/path.$set" if !defined $queryDir;
foreach my $dir ($stepDir, $queryDir, $resultsDir, "$resultsDir/$orgDir") {
  die "No such directory: $dir\n"
    unless -d $dir;
}

my @pathInfo = ReadTable("$stepDir/$set.table", ["pathwayId","desc"]);
my @pathwayIds = map { $_->{pathwayId} } grep { $_->{pathwayId} ne "all" } @pathInfo;
my %stepObj = (); # pathwayId => object containing steps and rules
foreach my $pathwayId (@pathwayIds) {
  $stepObj{$pathwayId} = ReadSteps("$stepDir/$pathwayId.steps");
}

my $reqs = ReadReqs("$stepDir/requires.tsv", \%stepObj);

my $doneFile = "$resultsDir/$orgDir/$set.sum.done";
my $dateFile = "$queryDir/date";
die "Organisms in $orgDir are not up to date\n"
  unless NewerThan($doneFile, $dateFile);

my @orgs = ReadOrgTable("$resultsDir/$orgDir/orgs.org");
my @sumSteps = ReadTable("$resultsDir/$orgDir/$set.sum.steps",
                         qw{orgId gdb gid pathway step onBestPath score locusId sysName
                            score2 locusId2 sysName2});
my %sumSteps = ();
foreach my $row (@sumSteps) {
  push @{ $sumSteps{$row->{orgId}} }, $row;
}
my @sumRules = ReadTable("$resultsDir/$orgDir/$set.sum.rules",
                         qw{orgId gdb gid pathway rule score nHi nMed nLo expandedPath});
my %sumRules = ();
foreach my $row (@sumRules) {
  push @{ $sumRules{$row->{orgId}} }, $row;
}

print join("\t", qw{orgId pathway rule requiredPath requiredRule requiredStep not genomeName comment})."\n";
foreach my $org (@orgs) {
  my $orgId = $org->{orgId};
  my $warns = CheckReqs($sumRules{$orgId}, $sumSteps{$orgId}, $reqs);
  foreach my $warn (@$warns) {
    print join("\t", $orgId, $warn->{pathway}, $warn->{rule},
               $warn->{requiredPath}, $warn->{requiredRule} || "", $warn->{requiredStep} || "",
               $warn->{not}, $warn->{genomeName} || "", $warn->{comment})."\n";
  }
}
