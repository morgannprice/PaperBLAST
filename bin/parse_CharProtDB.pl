#!/usr/bin/perl -w
# Parse the CharProtDB database (in SwissProt format) and convert into
# the PaperBLAST "curated_parsed" format.
# Ignores trusted_uniprot entries (as redundant with SwissProt)
# and trusted_aspgd entries (there were only a few of these and they
# seemed less functionally specific).
#
# The output format is:
# CharProtDB, protein identifier, secondary identifier, gene name, description, organism,
# sequence, comment (?), pmids

use strict;
use lib "SWISS/lib";
use SWISS::Entry;

die "Run as a filter: parse_CharProtDB.pl < charprotdb_trusted.sp > CharProtDB.curated_parsed\n"
  if @ARGV > 0;

# Read an entire record at a time
local $/ = "\n//\n";

while(<STDIN>) {
  # Ignore trusted_uniprot or trusted_aspgd entries
  next if m/ trusted_uniprot;/ || m/ trusted_aspgd;/;
  my $text = $_;
  my $entry = SWISS::Entry->fromText($text);
  my $id = $entry->ID;
  my $acc = $entry->AC;
  my @acc = split / /, $acc;
  shift @acc if @acc > 0 && $acc[0] eq $id;
  shift @acc if @acc > 0 && $acc[0] eq "|";

  my @DE = @{ $entry->DEs->list };
  my @desc = map { $_->text } grep { $_->category eq "RecName" } @DE;
  my $desc = join("; ", @desc);
  my $comment = "";
  my @pmids = ();
  foreach my $ref (@{ $entry->Refs->list }) {
    push @pmids, @{ $ref->RX->{PMID} }
      if exists $ref->RX->{PMID};
  }
  my @os = @{ $entry->OSs->list() };
  my @cc = map $_->comment(), $entry->CCs->elements;
  my $cc = join(";_", @cc); # in practice there is just one and there is no clear topic
  print join("\t", $id, $acc[0] || "",
             $entry->GNs->getFirst() || "",
             $desc,
             $os[0]->text() || "",
             $entry->SQ,
             $cc,
             join(",", @pmids)
            )."\n";
}
