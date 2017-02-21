#!/usr/bin/perl -w

use strict;

my $usage = <<END
Usage: removeDupQueries queryprot_files > unique_queries
    Removed queries from the arguments that match an earlier query.
    Issues warnings if queries have different sequences.
    If the description field is present, it maintains it.
END
;

die $usage if @ARGV == 0;
my @files = @ARGV;

my %known = (); # organism_:::_query => [id, sequence]
my $sep = "_:::_";

while(my $line = <>) {
  chomp $line;
  my ($organism, $query, $id, $seq, $desc) = split /\t/, $line;
  my $genus = $organism; $genus =~ s/ .*//;
  die "Invalid line:\n$line\nin $ARGV" unless defined $seq && $seq ne "";
  $desc = "" if !defined $desc;
  my $key = join($sep, $genus, $query);
  if (exists $known{$key}) {
    my ($oldid, $oldseq) = @{ $known{$key} };
    my $oldlen = length($oldseq);
    my $newlen = length($seq);
    print STDERR "Warning: sequence mismatch for $query in genus $genus -- $oldid vs. $id, $oldlen vs. $newlen a.a.\n"
      if $oldseq ne $seq;
  } else {
    print $line."\n";
    $known{$key} = [ $id, $seq ];
  }
}

