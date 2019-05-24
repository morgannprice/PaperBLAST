#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use lib "$Bin/../lib/SWISS";
use SWISS::Entry;

die "Usage: zcat uniprot_sprot.dat.gz | sprotCuratedEc.pl > curatedEc.faa\n"
  unless @ARGV == 0;

# Read an entire record at a time
local $/ = "\n//\n";

while(my $text = <STDIN>) {
  next if $text =~ m/Full=Uncharacterized/i || $text =~ m/Full=Putative uncharacterized/i;
  next unless $text =~ m/Bacteria/ || $text =~ m/Archaea/;
  my $entry = SWISS::Entry->fromText($text);
  my @OCs = $entry->OCs->elements;
  next unless $OCs[0] eq "Bacteria" || $OCs[0] eq "Archaea";
  my $desc = join("; ", map { $_->text } $entry->DEs->elements);
  next unless $desc =~ m/ EC /;
  my $FTs = $entry->FTs->list;
  my $isFragment = 0;
  foreach my $ft (@$FTs) {
    if ($ft->[0] eq "NON_TER") {
      print STDERR $entry->AC . " is a fragment, skipped\n";
      $isFragment = 1;
    }
  }
  next if $isFragment;
  my $hasCaution = 0;
  foreach my $CC ($entry->CCs->elements) {
    if ($CC->topic eq "CAUTION") {
      print STDERR $entry->AC . " has a caution, skipped\n";
      $hasCaution = 1;
    }
  }
  next if $hasCaution;
  print ">" . $entry->AC . " $desc\n" . $entry->SQ . "\n";
}
