#!/usr/bin/perl -w
# "Reverse search" -- compare candidates to curated proteins
# to see if they have better matches to other (unexpected) proteins.

use strict;
use Getopt::Long;
use List::Util qw{min max};
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{ReadTable ReadFastaEntry};
use Steps qw{ReadOrgTable ReadOrgProtein};

{
  my $nCPU = $ENV{MC_CORES} || 4;
  my @infields = qw{locusId type curatedId bits locusBegin locusEnd cBegin cEnd cLength identity};
  my @outfields = qw{locusId otherId bits locusBegin locusEnd otherBegin otherEnd otherIdentity};

  my $usage = <<END
Usage: gaprevsearch.pl -hits hitsfile -orgs orgsprefix -curated curated.udb -out outfile

Write to outfile, tab-delimited with 1 row per entry from the hits file.
The hits file must includes fields
  @infields
This script writes a new table with the fields
  @outfields
with all the hits for each hitId in the curated database

The curated file (in fasta or udb format) should be from
  running curatedFaa.pl with -curatedids

Optional arguments:
 -nCPU $nCPU -- number of CPUs to use (defaults to the MC_CORES
                environment variable, or 4)
END
;

  my ($hitsFile, $orgprefix, $curatedFile, $outFile);
  die $usage
    unless GetOptions('hitsFile=s' => \$hitsFile,
                      'orgs=s' => \$orgprefix,
                      'curated=s' => \$curatedFile,
                      'out=s' => \$outFile,
                      'nCPU=i' => \$nCPU)
      && @ARGV == 0
      && defined $hitsFile && defined $orgprefix && defined $outFile;
  foreach my $file ($hitsFile, $curatedFile) {
    die "No such file: $file\n" unless -e $file;
  }
  my $usearch = "$Bin/usearch";
  foreach my $x ($usearch) {
    die "No such executable: $x\n" unless -x $x;
  }
  die "Must specify at least 1 CPU\n" unless $nCPU >= 1;

  my @hits = ReadTable($hitsFile, \@infields);
  my %loci = ();
  foreach my $hit (@hits) {
    $loci{ $hit->{locusId} } = 1;
  }

  my @orgs = ReadOrgTable("$orgprefix.org");
  die "The organism table $orgprefix.org has no rows\n" unless @orgs > 0;

  my $aaIn = "$orgprefix.faa";
  die "No such file: $aaIn\n" unless -e $aaIn;
  open(my $fhIn, "<", $aaIn) || die "Cannot read $aaIn\n";
  my $faaCand = "/tmp/gaprevsearch.$$.db";
  open(my $fhCand, ">", $faaCand) || die "Cannot write to $faaCand\n";
  my $state = {};
  while (my $prot = ReadOrgProtein($fhIn,$state)) {
    my $locusId = $prot->{orgId} . ":" . $prot->{locusId};
    print $fhCand ">${locusId}\n$prot->{aaseq}\n"
      if exists $loci{$locusId};
  }
  close($fhIn) || die "Error reading $aaIn\n";
  close($fhCand) || die "Error writing to $faaCand\n";

  my $rhitsFile = "/tmp/gaprevsearch.$$.revhits";
  my $cmd = "$usearch -ublast $faaCand -db $curatedFile -id 0.3 -evalue 0.01 -blast6out $rhitsFile -threads $nCPU >& /dev/null";
  system($cmd) == 0 || die "Error running $cmd: $!\n";
  unlink($faaCand);
  open(my $fhRH, "<", $rhitsFile) || die "Cannot read from $rhitsFile\n";
  open(my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
  print $fhOut join("\t", @outfields)."\n";
  while (my $line = <$fhRH>) {
    chomp $line;
    my ($locusId, $otherId, $identity, $alen, $mm, $gap, $locusBegin, $locusEnd, $otherBegin, $otherEnd, $eval, $bits)
      = split /\t/, $line;
    print $fhOut join("\t", $locusId, $otherId, $bits, $locusBegin, $locusEnd, $otherBegin, $otherEnd, $identity)."\n";
  }
  close ($fhRH) || die "Error reading from $rhitsFile\n";
  close ($fhOut) || die "Error writing to $outFile\n";
  unlink($rhitsFile);
}
