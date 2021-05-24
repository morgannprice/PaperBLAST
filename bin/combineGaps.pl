#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};

my $usage = <<END
bin/combineGaps.pl -orgs orgSpec1 ... orgSpecN -out orgSpec -set aa
  Combine the GapMind results for org specifier 1 to N.
  Specifically, updates set.sum.{db,cand,steps,rules,warn,knownsim}

  Usually run after combineOrgs.pl.
END
;

my @orgs;
my ($set, $orgOut);

die $usage
  unless GetOptions('orgs=s{1,}' => \@orgs,
                    'out=s' => \$orgOut,
                    'set=s' => \$set)
  && @orgs > 0 && defined $set && defined $orgOut;

my $tmpDir = "$RealBin/../tmp";
my $useKnownSim = -e "$tmpDir/$orgs[0]/$set.sum.knownsim";
my $useKnownSimString = $useKnownSim ? "with" : "without";
print STDERR "Building $set for $orgOut\nfrom @orgs\n$useKnownSimString knownsim\n";

my @tableSuffix = ("cand","steps","rules","warn");
push @tableSuffix, "knownsim" if $useKnownSim;

my $stepsDb = "$tmpDir/path.$set/steps.db";
die "No such file: $stepsDb\n" unless -e $stepsDb;

foreach my $orgSpec (@orgs) {
  die "Input and output org specifiers should not match!\n"
    if $orgSpec eq $orgOut;
  my $pre = "$tmpDir/$orgSpec/orgs";
  foreach my $suffix ("org", "faa") {
    my $file = "$pre.$suffix";
    die "No such file: $file\n" unless -e $file;
  }
  my $sum = "$tmpDir/$orgSpec/$set.sum";
  foreach my $suffix ("db", @tableSuffix) {
    my $file = "$sum.$suffix";
    die "No file: $file\n" unless -e $file;
  }
}

my $outDir = "$tmpDir/$orgOut";
foreach my $suffix ("org", "faa") {
  my $file = "$outDir/orgs.$suffix";
  die "No such file: $file\n" unless -e $file;
}

my $combTables = "$RealBin/combineTables.pl";
die "No such executable: $combTables\n" unless -x $combTables;

foreach my $suffix (@tableSuffix) {
  my @inFiles = map "$tmpDir/$_/$set.sum.$suffix", @orgs;
  my $outFile = "$tmpDir/$orgOut/$set.sum.$suffix";
  system($combTables, "-in", @inFiles, "-out", $outFile) == 0
    || die "$combTables failed to build $outFile -- $!";
  print STDERR "Wrote $outFile\n";
}


my @buildCmd = ("$RealBin/buildGapsDb.pl",
                "-gaps", "$outDir/$set.sum",
                "-requirements", "$outDir/$set.sum.warn",
                "-steps", $stepsDb,
                "-out", "$outDir/$set.sum.db");
push @buildCmd, ("-markersim", "$outDir/$set.sum.knownsim")
  if $useKnownSim;
system(@buildCmd) == 0
  || die "@buildCmd\nfailed -- $!";
