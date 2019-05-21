#!/usr/bin/perl -w
use strict;
use lib "/usr2/people/mprice/workbook/PaperBLAST/SWISS/lib";
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
  print ">" . $entry->AC . " $desc\n" . $entry->SQ . "\n";
}
