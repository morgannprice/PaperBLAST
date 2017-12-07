#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ParsePTools, ReadFasta

my $usage = <<END
parse_ecocyc.pl proteins.dat protseq.fsa > ecocyc.curated_parsed

The output file is in curated_parsed format with fields ecocyc, ecocyc
identifier, bnumber, gene name, description, organism, comment
(blank), and pubmed ids (comma separated).
END
;

die $usage unless @ARGV == 2;
my ($protfile, $faafile) = @ARGV;

# a hash, with keys of the form gnl|ECOLI|id
my $seqs = ReadFasta($faafile);
my $nNoSeq = 0;
my $nPrint = 0;

open(my $fh, "<", $protfile) || die "Cannot read $protfile";
while (my $prot = ParsePTools($fh)) {
  next unless $prot->{TYPES}[0]{value} eq "Polypeptides";
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
