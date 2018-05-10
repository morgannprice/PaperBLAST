#!/usr/bin/perl -w
# Given the interesting sequences (from interesting_seqs.pl) and tab-delimited blast hits
# (with those as the input),
# report regions that have no hits
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ReadFastaEntry

my $minIdentity = 30;
my $minRegion = 50;
my $usage = <<END
interesting_nohits.pl [ -minIdentity $minIdentity ] [ -minRegion $minRegion ] -hits blastout.tab -int int.faa
END
;

my ($hitsfile, $intfile);
die $usage
  unless GetOptions('hits=s' => \$hitsfile,
                    'int=s' => \$intfile,
                    'minIdentity=f' => \$minIdentity,
                    'minRegion=i' => \$minRegion)
  && @ARGV == 0
  && defined $hitsfile && defined $intfile;
die "No such file: $hitsfile\n" unless -e $hitsfile;
die "No such file: $intfile\n" unless -e $intfile;

my %seq = (); # sequence id to sequence
my %seqlen = (); # sequence id to length
my %seqInt = (); # sequence id to #pubs

open(my $fhInt, "<", $intfile) || die "Cannot read $intfile";
my $state = {};
while(my ($header,$sequence) = ReadFastaEntry($fhInt,$state)) {
  my ($id, $n) = split / /, $header;
  die $header unless $n > 0;
  $seq{$id} = $sequence;
  $seqlen{$id} = length($sequence);
  $seqInt{$id} = $n;
}
close($fhInt) || die "Error reading $intfile";
print STDERR "Read " . scalar(keys %seqlen) . " sequences from $intfile\n";

my %ranges = (); # sequence id to list of [qBeg,qEnd]

open(my $fhHits, "<", $hitsfile) || die "Cannot read $hitsfile";
while (my $line = <$fhHits>) {
  chomp $line;
  my ($query, $subject, $identity, $alen, $mismatch, $gap, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $bits) = split /\t/, $line;
  die $line unless defined $bits;
  die "Unknown query $query" unless exists $seqlen{$query};
  my $qLen = $seqlen{$query};
  die "Invalid query end in $line" unless $qEnd <= $qLen;
  next unless $identity >= $minIdentity;
  push @{ $ranges{$query} }, [ $qBeg, $qEnd ];
}
close($fhHits) || die "Error reading $hitsfile";

# List no-hit regions, prioritized by #pubs
my %intSeqs = ();
while (my ($seq, $int) = each %seqInt) {
  push @{ $intSeqs{$int} }, $seq;
}

print join("\t", qw{id nPapers length begin end subseq})."\n";
my @intValues = sort { $b <=> $a } keys %intSeqs;
foreach my $int (@intValues) {
  foreach my $seq (sort @{ $intSeqs{$int} }) {
    my @out = (); # list of uncovered ranges
    if (!exists $ranges{$seq}) {
      push @out, [1, $seqlen{$seq} ]
        if $seqlen{$seq} >= $minRegion;
    } else {
      # figure out what regions, of at least minRegion in size, are not covered and print out those
      my @ranges = sort { $a->[0] <=> $b->[0] } @{ $ranges{$seq} };
      my $at = 1;
      foreach my $range (@ranges) {
        my ($beg,$end) = @$range;
        push @out, [ $at, $beg-1 ] if $beg >= $at + $minRegion;
        $at = $end+1 unless $at > $end+1;
      }
      push @out, [ $at, $seqlen{$seq} ]
        if $seqlen{$seq} >= $at + $minRegion;
    }
    foreach my $range (@out) {
      my ($beg,$end) = @$range;
      die if $end - $beg + 1 < $minRegion;
      print join("\t", $seq, $int, $seqlen{$seq},
                 $beg, $end,
                 substr($seq{$seq}, $beg-1, $end-$beg+1))."\n";
    }
  }
}
