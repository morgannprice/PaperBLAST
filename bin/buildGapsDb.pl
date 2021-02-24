#!/usr/bin/perl -w
# Build the gaps database from the output of gapsummary.pl

use strict;
use Getopt::Long;
use DBI;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadTable SqliteImport};

my $usage = <<END
buildGapsDb.pl -gaps aa.sum -requirements warningsFile -markersim aa.known.sim -out gaps.db
  Given the output files from gapsummary.pl, checkGapRequirements.pl,
  and orgsVsMarkers.pl, build a sqlite3 gaps database.
END
;

{
  my ($gapsPre, $warningsFile, $markerSimFile, $outDb);
  die $usage
    unless GetOptions('gaps=s' => \$gapsPre,
                      'requirements=s' => \$warningsFile,
                      'markersim=s' => \$markerSimFile,
                      'out=s' => \$outDb)
      && @ARGV == 0
      && defined $gapsPre && defined $warningsFile && defined $markerSimFile && defined $outDb;

  my $schemaFile = "$RealBin/../lib/gaps.sql";
  foreach my $file ("$gapsPre.cand", "$gapsPre.steps", "$gapsPre.rules",
                    $warningsFile, $markerSimFile, $schemaFile) {
    die "No such file: $file\n" unless -e $file;
  }

  my @candFields = qw{orgId pathway step score locusId sysName desc locusId2 sysName2 desc2 blastBits curatedIds identity blastCoverage blastScore curatedDesc hmmBits hmmId hmmCoverage hmmScore hmmDesc otherIds otherBits otherIdentity otherCoverage};
  my @cands = ReadTable("$gapsPre.cand", \@candFields);
  my @candidateOutFields = qw{orgId pathway step locusId sysName desc locusId2 sysName2 desc2 score blastBits curatedIds identity blastCoverage blastScore hmmId hmmDesc hmmBits hmmCoverage hmmScore otherIds otherBits otherIdentity otherCoverage};
  my @candidateOut = ();
  foreach my $row (@cands) {
    my @rowOut = map $row->{$_}, @candidateOutFields;
    push @candidateOut, \@rowOut;
  }

  my @stepFields = qw{orgId pathway step onBestPath score locusId sysName score2 locusId2 sysName2};
  my @steps = ReadTable("$gapsPre.steps", \@stepFields);
  my @stepOutFields = qw{orgId pathway step onBestPath score locusId sysName score2 locusId2 sysName2};
  my @stepOut = ();
  foreach my $row (@steps) {
    my @rowOut = map $row->{$_}, @stepOutFields;
    push @stepOut, \@rowOut;
  }

  my @ruleFields = qw{orgId pathway rule nHi nMed nLo score expandedPath path path2};
  my @rules = ReadTable("$gapsPre.rules", \@ruleFields);
  my @ruleOutFields = qw{orgId pathway rule nHi nMed nLo score expandedPath path path2};
  my @ruleOut = ();
  foreach my $row (@rules) {
    my @rowOut = map $row->{$_}, @ruleOutFields;
    push @ruleOut, \@rowOut;
  }

  my @warningFields = qw{orgId pathway rule requiredPath requiredRule requiredStep not comment};
  my @warnings = ReadTable($warningsFile, \@warningFields);
  my @warningOutFields = qw{orgId pathway rule requiredPath requiredRule requiredStep not comment};
  my @warningOut = ();
  foreach my $row (@warnings) {
    my @rowOut = map $row->{$_}, @warningOutFields;
    push @warningOut, \@rowOut;
  }

  my @markerSimFields = qw{orgId orgId2 identity nMarkers};
  my @markerSim = ReadTable($markerSimFile, \@markerSimFields);
  my @markerSimOutFields = qw{orgId orgId2 identity nMarkers};
  my @markerSimOut = ();
  foreach my $row (@markerSim) {
    my @rowOut = map $row->{$_}, @markerSimOutFields;
    push @markerSimOut, \@rowOut;
  }

  my $tmpDir = $ENV{TMP} || "/tmp";
  my $tmpDbFile = "$tmpDir/buildGapsDb.$$.db";
  print STDERR "Building temporary database $tmpDbFile\n";
  unlink($tmpDbFile);
  system("sqlite3 $tmpDbFile < $schemaFile") == 0
    || die "Error loading schema $schemaFile into $tmpDbFile -- $!";

  SqliteImport($tmpDbFile, "Candidate", \@candidateOut);
  SqliteImport($tmpDbFile, "StepScore", \@stepOut);
  SqliteImport($tmpDbFile, "RuleScore", \@ruleOut);
  SqliteImport($tmpDbFile, "RequirementNotMet", \@warningOut);
  SqliteImport($tmpDbFile, "MarkerSimilarity", \@markerSimOut);

  system("cp $tmpDbFile $outDb") == 0 || die "Copying $tmpDbFile to $outDb failed: $!";
  unlink($tmpDbFile);
  print STDERR "Built gaps database $outDb\n";
}

