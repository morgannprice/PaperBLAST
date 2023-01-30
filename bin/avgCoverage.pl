#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use List::Util qw{sum};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry};

die "Run as a filter on fasta scaffolds to compute average coverage, weighted by scaffold lengths, (or NA)\n"
  unless @ARGV == 0;

my @cov = ();
my @len = ();
my $state = {};
while (my ($header, $sequence) = ReadFastaEntry(\*STDIN, $state)) {
  push @len, length($sequence);
  if ($header =~ m/_cov_([0-9.]+)$/) {
    push @cov, $1;
  } else {
    push @cov, "";
  }
}

my $avgCov = "NA";
if (scalar(@len) == 0) {
  # (not reachable given current implementation of ReadFastaEntry
  print STDERR "Warning: input is empty\n";
} elsif (scalar(grep $_ eq "", @cov) > 0) {
  print STDERR "Warning: cannot parse coverage from scaffold names\n";
} else {
  $avgCov = 0;
  for (my $i = 0; $i < scalar(@cov); $i++) {
    $avgCov +=$len[$i] * $cov[$i];
  }
  $avgCov /= sum(@len);
}
print $avgCov,"\n";
