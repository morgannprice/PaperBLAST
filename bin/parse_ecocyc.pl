#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ParsePTools, ReadFasta

my $usage = <<END
Usage: parse_ecocyc.pl ecocyc_data_dir > ecocyc.curated_parsed

Reads proteins.dat, enzrxns.dat and protseq.fsa from the ecocyc data directory.

The output file is in curated_parsed format with fields ecocyc, ecocyc
identifier, bnumber, gene name, description, organism, comment
(blank), and pubmed ids (comma separated).
END
;

die $usage unless @ARGV == 1;
my ($dir) = @ARGV;
die "Not a directory: $dir\n$usage" unless -d $dir;

my $protfile = "$dir/proteins.dat";
my $enzrxnsfile = "$dir/enzrxns.dat";
my $faafile = "$dir/protseq.fsa";

# a hash, with keys of the form gnl|ECOLI|id
my $seqs = ReadFasta($faafile);
my $nNoSeq = 0;
my $nPrint = 0;

my %enzrxnName = (); # enzrxnId to common name
open(my $fhEnz, "<", $enzrxnsfile) || die "Cannot read $enzrxnsfile";
while (my $enzrxn = ParsePTools($fhEnz)) {
  my $id = $enzrxn->{"UNIQUE-ID"}[0]{value};
  my $name = $enzrxn->{"COMMON-NAME"}[0]{value};
  $enzrxnName{$id} = $name if defined $id && defined $name;
}
close($fhEnz) || die "Error reading $enzrxnsfile";

open(my $fh, "<", $protfile) || die "Cannot read $protfile";
while (my $prot = ParsePTools($fh)) {
  # does it have TYPES = Polypeptides (not necessarily as the first value)
  my @types = map { $_->{value} } @{ $prot->{TYPES} };
  my @poly = grep { $_ eq "Polypeptides" } @types;
  next unless @poly > 0;

  my $id = $prot->{"UNIQUE-ID"}[0]{value};
  die "No unique identifier" unless $id;
  my $key = "gnl|ECOLI|$id";
  if (!exists $seqs->{$key}) {
    # Some of the "proteins" are pseudogenes, or known only from biochemical studies,
    # and do not have sequences.
    $nNoSeq++;
    next;
  }
  my $seq = $seqs->{$key};

  my $bnum = "";
  my @dblinks = map { $_->{value} } @{ $prot->{"DBLINKS"} };
  my @blinks = grep m/ECOLIWIKI/, @dblinks;
  if (@blinks == 1) {
    $blinks[0] =~ m/"([bB]\d+)"/ || die "Cannot find b number in $blinks[0]";
    $bnum = $1;
  }

  my $name = $prot->{"ABBREV-NAME"}[0]{value} || $prot->{"SYNONYMS"}[0]{value} || "";
  my $desc = $prot->{"COMMON-NAME"}[0]{value} || "";
  if ($desc eq "" && exists $prot->{"CATALYZES"}) {
    my $enzrxnId = $prot->{"CATALYZES"}[0]{value};
    $desc = $enzrxnName{ $enzrxnId } if $enzrxnId && exists $enzrxnName{ $enzrxnId };
  }

  # Links to papers may be in the CITATIONS attribute or in the CITATIONS annotation of GO-TERMS
  # or in the COMMENT field as entries like |CITS: [8621577]| or |CITS: [2693216][1425658]| or |CITS: [9755155] [9388228]|
  my @pmLines = ();
  foreach my $value ( @{ $prot->{"CITATIONS"} }) {
    push @pmLines, $value->{"value"};
  }
  foreach my $value ( @{ $prot->{"GO-TERMS"} }) {
    push @pmLines, $value->{CITATIONS}
      if exists  $value->{CITATIONS};
  }
  my %pmIds = ();
  foreach my $line (@pmLines) {
    my $pmId = $1
      if $line =~ m/^(\d+):/ || $line =~ m/^(\d+)$/ || $line =~ m/^\[(\d+)\]$/;
    # Some citations are like Wang02, not sure why
    $pmIds{$pmId} = 1 if $pmId;
  }
  my $comment = $prot->{"COMMENT"}[0]{"value"} || "";
  while ($comment =~ m/[|]CITS: ([\[\]0-9 ]+)[|]/g) {
    # this should extract a field like [8801422] or [10501935][11872485] or [6706930] [7142155]
    my $ids = $1;
    $ids =~ s/[\[\]]/ /g;
    my @ids = split / +/, $ids;
    foreach my $id (@ids) {
      $pmIds{$id} = 1 if $id;
    }
  }
  my @pmIds = sort { $a <=> $b } keys %pmIds;

  print join("\t", "ecocyc", $id,
             $bnum, $name,
             $desc,
             "Escherichia coli K-12 substr. MG1655",
             $seq,
             "", # ocmment
             join(",",@pmIds))."\n";
  $nPrint++;
}
close($fh) || die "Error reading $protfile";
print STDERR "Found $nPrint polypeptides with sequences and $nNoSeq without\n";
