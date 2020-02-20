#!/usr/bin/perl
# Given a stream of SwissProt entries, and a curated parsed table,
# create a new table that includes the SUBUNIT entry as well as isHetero
use strict;
use lib "SWISS/lib";
use Getopt::Long;
use SWISS::Entry;

my $usage = <<END
Usage: zcat uniprot_sprot.dat.gz | sprotSubunit.pl -curated sprot.curated_parsed > withsubunit.tab

The output is tab-delimited with a header line, and with the fields
db, id, desc, organism, isHetero, subunit

The subunit field is the first sentence of the SUBUNIT comment (and is
often empty).  Because the subunit field is free text, the
identification of heteromeric proteins ("isHetero") is not perfectly
accurate.
END
;

my $curatedFile;
die $usage
  unless GetOptions('curated=s' => \$curatedFile)
  && @ARGV == 0
  && defined $curatedFile;

my @curatedRows = ();
open(my $fhC, "<", $curatedFile) || die "Cannot read $curatedFile\n";
while(my $line = <$fhC>) {
  my ($db, $id, undef, undef, $desc, $org) = split /\t/, $line;
  die "Invalid curated parsed line:\n$line"
    unless defined $org && $org ne "" && $desc ne "" && $id =~ m/^[A-Z0-9_]+$/;
  push @curatedRows, { 'db' => $db, 'id' => $id, 'desc' => $desc, 'org' => $org };
}
close($fhC) || die "Error reading $curatedFile\n";
my %curatedRows = map { $_->{id} => $_ } @curatedRows;

# Read an entire record at a time
local $/ = "\n//\n";

print join("\t", "db", "id", "desc", "organism", "isHetero", "subunit")."\n";

my $nWritten = 0;
my $nHetero = 0;
while(my $text = <STDIN>) {
  my $entry = SWISS::Entry->fromText($text);
  my $id = $entry->AC;
  next unless exists $curatedRows{$id};
  my $crow = $curatedRows{$id} || die;
  my @sub  = map { $_->comment} grep { $_->topic eq 'SUBUNIT'} $entry->CCs->elements;
  my $subunit = join("; ", @sub);
  $subunit =~ s/[.].*//; # keep just the first sentence

  # matching any of these indicates is a heteromer
  my $matchesHetero = $subunit =~ m/hetero[a-z-]*mer/i
    || $subunit =~ m/forms a hetero[a-z-]*mer/i;
  # matching any of these indicates it is a homomer
  my $matchesHomo = $subunit eq ""
    || $subunit =~ m/monomer/i
    || $subunit =~ m/homo[a-z]+mer/i
    || $subunit =~ m/^homo/i
    || $subunit =~ m/^interacts/i # this usually comes after any heteromore information
    || $subunit =~ m/^dimer$/i
    || $subunit =~ m/^dimer of (identical|homo|dimer|trimer|tetramer)/i
    || $subunit =~ m/^binds/i;
  my $isHetero = $matchesHetero && ! $matchesHomo;
  $nHetero++ if $isHetero;
  print join("\t", $crow->{db}, $crow->{id}, $crow->{desc}, $crow->{org},
             $isHetero ? 1 : 0,
             $subunit)."\n";
  $nWritten++;
}
print STDERR "Wrote " . scalar($nWritten) . " entries. Found $nHetero heteromers. Curated parsed had " . scalar(@curatedRows) . " entries\n";
