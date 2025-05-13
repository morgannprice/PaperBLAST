#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry};

die "Usage: pdbDesc.pl pdbfile protnames.lst > pdb.faa\n"
  unless @ARGV == 2;
my ($pdbFile, $namesFile) = @ARGV;
open(my $fhSeq, "<", $pdbFile) || die "Cannot read $pdbFile\n";
open(my $fhName, "<", $namesFile) || die "Cannot read $namesFile\n";

my %pdbDesc;
while (my $line = <$fhName>) {
  chomp $line;
  # example input line:
  # 100d - 31-Mar-95 X-ray   1.900 Crystal structure of the highly distorted chimeric decamer r(c)d(cggcgccg)r(g)-spermine complex-spermine binding to phosphate only and minor groove tertiary base-pairing
  my @parts = split /\s+/, $line;
  if (@parts < 6) {
    print STDERR "Warning: cannot parse pdb names from $line\n"
      unless $line =~ m/^[a-z0-9]+ +[*-] +[-] +[-] *$/; # no name
    next;
  }
  my $id = $parts[0];
  $pdbDesc{$id} = join(" ", splice(@parts, 5));
}
close($fhName) || die "Error reading $namesFile\n";

my $state = {};
while (my ($name, $seq) = ReadFastaEntry($fhSeq, $state)) {
  next unless $name =~ m/mol:protein/ && length($seq) >= 20;
  $name =~ m/^([0-9a-z]+)_([0-9A-Za-z]+) +mol:protein length:\d+ +(.*)$/
    || die "Cannot handle name $name";
  my ($id, $chain, $desc) = ($1,$2,$3);
  $desc = $pdbDesc{$id} if exists $pdbDesc{$id};
  print ">${id}_${chain} $desc\n$seq\n";
}
close($fhSeq)|| die "Error reading $pdbFile\n";
