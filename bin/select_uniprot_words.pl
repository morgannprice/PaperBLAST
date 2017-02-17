#!/usr/bin/perl -w
# Given a list of words, select possible uniprot identifiers or accessions
# and output an input file for sprotToQuery.pl
use strict;
die "Run as a filter\n" unless @ARGV == 0;

while(<STDIN>) {
    chomp;
    my @F = split /\t/, $_;
    print $F[1]."\n"
      if @F >= 2
      && $F[1] =~ m/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$/;
}
