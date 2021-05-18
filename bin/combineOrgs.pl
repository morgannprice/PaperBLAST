#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry};
use Steps qw{ReadOrgTable};

my $usage = <<END
combineOrgs.pl -in orgs1 orgs2 ... orgsN -out orgsCombined
  Combines one or more organism tables with .org and .faa suffixes
END
;

my @orgs;
my $outPre;
die $usage
  unless GetOptions('in=s{1,}' => \@orgs,
                    'out=s' => \$outPre)
  && defined $outPre && @orgs > 0 && @ARGV == 0;

my $comb = "$RealBin/combineTables.pl";
die "No such executable: $comb" unless -x $comb;
my $blastdir = "$RealBin/blast";
die "No such directory: $blastdir\n" unless -d $blastdir;
my $formatdb = "$blastdir/formatdb";
die "No such executable: $formatdb\n" unless -x $formatdb;

my @orgFiles = map { $_ . ".org" } @orgs;
my @faaFiles = map { $_ . ".faa" } @orgs;
foreach my $file (@orgFiles, @faaFiles) {
  die "No such file: $file" unless -e $file;
}

system($comb, "-out", "$outPre.org", "-in", @orgFiles) == 0
  || die "$comb failed: $!";

my %seen = (); # protein ids seen already
open(my $fhOut, ">", "$outPre.faa")
  || die "Cannot write to $outPre.faa";
foreach my $faaFile (@faaFiles) {
  open(my $fhIn, "<", $faaFile) || die "Cannot read $faaFile\n";
  my $state = {};
  my $nSeq = 0;
  while (my ($header,$sequence) = ReadFastaEntry($fhIn, $state)) {
    my $id = $header;
    $id =~ s/ .*//;
    die "Duplicate id $id from fasta input $faaFile"
      if exists $seen{$id};
    $seen{$id} = 1;
    print $fhOut ">$header\n$sequence\n";
    $nSeq++;
  }
  print STDERR "Read $nSeq proteins from $faaFile\n";
  close($fhIn) || die "Error reading $faaFile";
}
close($fhOut) || die "Error writing $outPre.faa";

system($formatdb, "-p", "T", "-i", "$outPre.faa", "-o", "T") == 0
  || die "$formatdb on $outPre.faa failed: $!\n";
