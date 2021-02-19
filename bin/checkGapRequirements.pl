#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use DBI;
use lib "$Bin/../lib";
use Steps qw{ReadSteps ReadOrgTable ReadReqs CheckReqs};
use pbutils qw{ReadTable NewerThan};

my $resultsDir = "$Bin/../tmp";
my $set = "aa";

my $usage = <<END
checkGapRequirements.pl -org orgsNif > warnings.tab

Given a pathway set and an org directory with computed results,
check against the dependency requirements and output a table of warnings.

Optional arguments:
-set -- default $set
-results -- defaults to $resultsDir
-stepsDb -- defaults to $resultsDir/path.set/steps.db
END
;

my ($orgDir,$stepsDb);
die $usage
  unless GetOptions('set=s' => \$set,
                    'org=s' => \$orgDir,
                    'stepsDb' => \$stepsDb,
                    'results=s' => \$resultsDir)
  && @ARGV == 0;
die "No organism directory specified\n" unless defined $orgDir;
$stepsDb = "$resultsDir/path.$set/steps.db";
foreach my $dir ($resultsDir, "$resultsDir/$orgDir") {
  die "No such directory: $dir\n"
    unless -d $dir;
}

my $dbhS = DBI->connect("dbi:SQLite:dbname=${stepsDb}","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $doneFile = "$resultsDir/$orgDir/$set.sum.done";
die "Organisms in $orgDir are not up to date:\n$doneFile\nshould be newer than\n$stepsDb\n"
  unless NewerThan($doneFile, $stepsDb);

my $reqs = $dbhS->selectall_arrayref("SELECT * FROM Requirement",
                                     { Slice => {} });

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
    print join("\t", $orgId, $warn->{pathwayId}, $warn->{ruleId},
               $warn->{requiredPathwayId}, $warn->{requiredRuleId} || "", $warn->{requiredStepId} || "",
               $warn->{isNot}, "", $warn->{comment})."\n";
    #XXX genomeName sould ultimately be remove
  }
}
