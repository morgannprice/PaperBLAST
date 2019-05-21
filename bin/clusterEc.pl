#!/usr/bin/perl -w
use strict;
use lib "lib";
use pbutils qw{ReadFastaEntry};

# Split any clusters that have inconsistent EC numbers, and save
# the final clustering as a fasta file

die "Usage: clusterEC.pl curatedEc.faa curatedEc.uc > clusterEc.faa\n"
  unless @ARGV == 2;

my ($faaFile, $clustFile) = @ARGV;

my @clusters = (); # cluster number to list of members

open(my $fhClust, "<", $clustFile) || die "Cannot read $clustFile";
while (my $line = <$fhClust>) {
  chomp $line;
  my ($type, $cluster, undef, undef, undef, undef, undef, undef, $queryDesc, $hitDesc) = split /\t/, $line;
  die "Invalid cluster line, not enough fields:\n$line\n"
    unless defined $hitDesc;
  die "Invalid cluster $cluster" unless $cluster =~ m/^\d+$/;
  if ($type eq "S" || $type eq "H") {
    push @{ $clusters[$cluster] }, $queryDesc;
  } else {
    die "Invalid cluster line, unknown type\n$line\n" unless $type eq "C";
  }
}
close($fhClust) || die "Error reading $clustFile";

my @splitClusters = ();
foreach my $cluster (@clusters) {
  my %byEc = ();
  foreach my $member (@$cluster) {
    my %ec = ();
    # allow names like .n2 for the final EC number component
    while ($member =~ m/ EC ([0-9]+[.][0-9-]+[.][0-9-]+[.][0-9a-z-]+)/g) {
      $ec{$1} = 1;
    }
    my $ec = join("::", sort keys %ec);
    push @{ $byEc{$ec} }, $member;
  }
  foreach my $key (sort keys %byEc) {
    # add the reference to the cluster
    push @splitClusters, $byEc{$key};
  }
}
print STDERR "Split " . scalar(@clusters) . " clusters by EC to " . scalar(@splitClusters) . " clusters\n";

# and save the first representative of each cluster
my %keep = (); # header to 1 if we should print this sequence
foreach my $cluster (@splitClusters) {
  die "Empty cluster\n" if @$cluster==0;
  die "Duplicate cluster for " . $cluster->[0] if exists $keep{$cluster->[0]};
  $keep{$cluster->[0]} = 1;
}

open(my $fhFaa, "<", $faaFile) || die "Cannot read $faaFile";
my $state = {};
my $nWrite = 0;
while (my ($header, $seq) = ReadFastaEntry($fhFaa, $state)) {
  if (exists $keep{$header}) {
    print ">$header\n$seq\n";
    $nWrite++;
  }
}
close($fhFaa) || die "Error reading $faaFile";
die "Wrote $nWrite sequences instead of the expected " . scalar(keys %keep) . " sequences -- input files are inconsistent?\n"
  unless $nWrite == scalar(keys %keep);
