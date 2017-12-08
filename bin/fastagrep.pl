#!/usr/bin/perl -w
use strict;

die "Usage: fastagrep.pl regexp < input.faa > subset.faa\n"
  unless @ARGV == 1;
my ($regexp) = @ARGV;

my $writing = 0;
while(my $line = <STDIN>) {
  if ($line =~ m/^>(.*)/) {
    my $desc = $1;
    $writing = $desc =~ $regexp;
  }
  print $line if $writing;
}
