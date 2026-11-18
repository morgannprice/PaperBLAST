#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry readPdbNames};

die "Usage: pdbDesc.pl pdbfile pdbNames.tsv > pdb.faa\n"
  unless @ARGV == 2;
my ($pdbFile, $namesFile) = @ARGV;
open(my $fhSeq, "<", $pdbFile) || die "Cannot read $pdbFile\n";
open(my $fhName, "<", $namesFile) || die "Cannot read $namesFile\n";

my %pdbDesc = %{ readPdbNames($namesFile) };

my $state = {};
while (my ($name, $seq) = ReadFastaEntry($fhSeq, $state)) {
  next unless $name =~ m/mol:protein/ && length($seq) >= 20;
  $name =~ m/^([0-9a-z]+)_([0-9A-Za-z]*) +mol:protein length:\d+ +(.*)$/
    || die "Cannot handle name $name";
  my ($id, $chain, $desc) = ($1,$2,$3);
  $desc = $pdbDesc{$id} if exists $pdbDesc{$id};
  print ">${id}_${chain} $desc\n$seq\n";
}
close($fhSeq)|| die "Error reading $pdbFile\n";
