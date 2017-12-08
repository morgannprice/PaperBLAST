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
    if (exists $ref->RX->{PMID}) {
      foreach my $part (@{ $ref->RX->{PMID} }) {
        $part =~ s/ //g; # occasionally a space is in there
        push @pmids, split /,/, $part;
      }
    }
  }
  my @filtered= ();
  foreach my $pmid (@pmids) {
    # remove common prefixes or suffixes
    $pmid =~ s/^PMID://;
    $pmid =~ s/[|]PMID$//;
    if ($pmid =~ m/^\d+$/) {
      push @filtered, $pmid;
    } else {
      print STDERR "Skipping invalid pubmed id $pmid for $id\n";
    }
  }

  my %pmids = map { $_ => 1 } @filtered;
  @pmids = sort { $a <=> $b } keys %pmids;

  my @os = @{ $entry->OSs->list() };
  my @cc = map $_->comment(), $entry->CCs->elements;
  my $cc = join("; ", @cc); # in practice there is just one and there is no clear topic
  # Fixing up multiple comments
  $cc =~ s/\nCC +/ /g;
  $cc =~ s/[\r\n\t]/ /g;

  print join("\t",
             "CharProtDB",
             $id, $acc[0] || "",
             $entry->GNs->getFirst() || "",
             $desc,
             $os[0]->text() || "",
             $entry->SQ,
             $cc,
             join(",", @pmids)
            )."\n";
}
