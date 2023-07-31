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
  my $maxKeep = 8;
  # By default, use the (more sensitive) seg qmask to avoid slow performance on some queries, such as
  # S layer protein IFHOHEBM_00707. (In rare cases this causes some moderately-strong alignments
  # to other hits to be missed.)
  my $qmask = "seg";
  my @infields = qw{locusId type queryId bits locusBegin locusEnd qBegin qEnd qLength identity};
  my @outfields = qw{locusId otherId bits locusBegin locusEnd otherBegin otherEnd otherIdentity};

  my $usage = <<END
Usage: gaprevsearch.pl -hits hitsfile -orgs orgsprefix -curated curatedDb -out outfile

Write to outfile, tab-delimited with 1 row per entry from the hits file.
The hits file must includes fields
  @infields
This script writes a new table with the fields
  @outfields
with all the hits for each hitId in the curated database

The curated file should be in fasta or udb format if using usearch, or
should be a diamond database if using diamond. The sequences should be
from running curatedFaa.pl with -curatedids

Optional arguments:
-maxKeep $maxKeep -- maximum number of curated hits to retain
-qmask $qmask -- qmask option for usearch
  (older versions of gaprevsearch did not set this option when running
   usearch/ublast, which is equivalent to -qmask fastamino)
-nCPU $nCPU -- number of CPUs to use (defaults to the MC_CORES
                environment variable, or 4)
-diamond -- use diamond instead of usearch
END
;

  my ($hitsFile, $orgprefix, $curatedFile, $outFile, $useDiamond);
  die $usage
    unless GetOptions('hitsFile=s' => \$hitsFile,
                      'orgs=s' => \$orgprefix,
                      'curated=s' => \$curatedFile,
                      'out=s' => \$outFile,
                      'qmask=s' => \$qmask,
                      'nCPU=i' => \$nCPU,
                     'diamond' => \$useDiamond)
      && @ARGV == 0
      && defined $hitsFile && defined $orgprefix && defined $outFile;
  foreach my $file ($hitsFile, $curatedFile) {
    die "No such file: $file\n" unless -e $file;
  }
  my $usearch = "$Bin/usearch";
  my $diamond = "$Bin/diamond";
  my @exe;
  if (defined $useDiamond) {
    push @exe, $diamond;
  } else {
    push @exe, $usearch;
  }
  foreach my $x (@exe) {
    die "No such executable: $x\n" unless -x $x;
  }
  die "Must specify at least 1 CPU\n" unless $nCPU >= 1;

  my @hits = ReadTable($hitsFile, \@infields);
  my $rhitsFile = "/tmp/gaprevsearch.$$.revhits";
  # usearch fails if there are no candidates so:
  if (scalar(@hits) == 0) {
    system("touch", $rhitsFile);
  } else {
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

    my $cmd;
    if (defined $useDiamond) {
      $cmd = "$diamond blastp --query $faaCand --db $curatedFile --id 0.3 --evalue 0.01 --out $rhitsFile --very-sensitive --outfmt 6 --masking $qmask --threads $nCPU";
    } else {
      $cmd = "$usearch -ublast $faaCand -db $curatedFile -id 0.3 -evalue 0.01 -blast6out $rhitsFile -qmask $qmask -threads $nCPU";
    }
    $cmd .= " > /dev/null 2>&1";
    system($cmd) == 0 || die "Error running $cmd: $!\n";
    unlink($faaCand);
  }
  open(my $fhRH, "<", $rhitsFile) || die "Cannot read from $rhitsFile\n";
  open(my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
  print $fhOut join("\t", @outfields)."\n";
  my %nHits = (); # locusId => number of hits so far
  while (my $line = <$fhRH>) {
    chomp $line;
    my ($locusId, $otherId, $identity, $alen, $mm, $gap, $locusBegin, $locusEnd, $otherBegin, $otherEnd, $eval, $bits)
      = split /\t/, $line;
    $nHits{$locusId}++;
    print $fhOut join("\t", $locusId, $otherId, $bits, $locusBegin, $locusEnd, $otherBegin, $otherEnd, $identity)."\n"
      if $nHits{$locusId} <= $maxKeep;
  }
  close ($fhRH) || die "Error reading from $rhitsFile\n";
  close ($fhOut) || die "Error writing to $outFile\n";
  unlink($rhitsFile);
}
