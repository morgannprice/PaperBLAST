#!/usr/bin/perl
# Given a large stream of SwissProt lines, produce a tab-delimited file
# that includes characterized items with links to rhea.
# It also writes out subunit information, to help resolve heteromeric complexes
# (but, does not actually do so)

use strict;
use lib "SWISS/lib";
use SWISS::Entry;

die "Run as a filter on swiss-prot entries\n" unless @ARGV==0;

# Read an entire record at a time
local $/ = "\n//\n";

while(my $text = <STDIN>) {
    my $entry = SWISS::Entry->fromText($text);
    my $DEs = $entry->DEs;
    my @desc = map { $_->text } grep { $_->type ne "EC" } $DEs->elements;
    # Add EC numbers
    my @ec = map { $_->text } grep { $_->type eq "EC" } $DEs->elements;
    foreach my $sublist ($DEs->Includes->elements) {
      push @ec, map { $_->text } grep { $_->type eq "EC" } $sublist->elements;
    }
    foreach my $sublist ($DEs->Contains->elements) {
      push @ec, map { $_->text } grep { $_->type eq "EC" } $sublist->elements;
    }
    # remove duplicates from @ec
    my %ecSoFar = ();
    @ec = grep { my $e = !exists $ecSoFar{$_}; $ecSoFar{$_} = 1; $e } @ec;
    # the EC elements already include an "EC" prefix. Just add to the ; separated list
    push @desc, @ec;
    my $desc = join("; ", @desc);

    my $organisms = join ("; ", map { $_->text } $entry->OSs->elements);
    # Each comment is a CC object with fields topic and comment
    # References to pubmed normally show up as something like
    # {ECO:0000269|PubMed:16731932, ECO:0000269|PubMed:19861680}
    # However sometimes the evidence code is in the special blocks section?
    my @cc = map $_->toString(), $entry->CCs->elements;
    my @cc2 = map { s/^CC +//; s/\nCC +/ /g; s/\n+//g; s/^-!- *//; $_; } @cc;
    @cc2 = grep !m/Copyrighted by/, @cc;
    my $cc = join("_:::_", @cc2);
    # which produces multiple fields like FUNCTION: ....
    # joined by the arbitrary _:::_ join.
    # Note that this includes all of the fields, and not just the most functionally relevant ones
    # (such as SUBUNIT, FUNCTION,COFACTOR, CATALYTIC, ENZYME, DISRUPTION).
    # Also it still includes the references.

    # Parse the pubmed identifiers from the comments lines
    my @pmIds = $cc =~ m!ECO:0000269[|]PubMed:(\d+)!g;
    my %pmIds = map { $_ => 1 } @pmIds;
    @pmIds = sort { $a <=> $b } keys %pmIds;

    # Find the rhea entries with experimental evidence
    my %rhea = ();
    foreach my $comment (@cc) {
      if ($comment =~ m/^CATALYTIC ACTIVITY/ && $comment =~ m/ECO:0000269/) {
        while($comment =~ /Xref=Rhea:(RHEA:\d+)/g) {
          $rhea{$1} = 1;
        }
      }
    }
    my $subunit = join("; ", grep m/^SUBUNIT/, @cc);
    $subunit =~ s/-!- SUBUNIT: +//g;
    my @rhea = sort keys %rhea;
    # AC = accession
    # If there are multiple accessions (due to merging), then AC will report just the first (primary) one,
    # which desirable (I think)
    # The first element of the IDs list is the primary ID
    # SQ = sequence
    foreach my $rhea (@rhea) {
      print join("\t", $rhea, $entry->AC, @{ $entry->IDs->list }[0] || "", $desc, $organisms,
                 $entry->SQ, join(",",@pmIds),
                 $subunit)."\n";
    }
}
