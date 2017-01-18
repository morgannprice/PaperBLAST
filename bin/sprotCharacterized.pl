#!/usr/bin/perl
# Given a large stream of SwissProt lines, produce a tab-delimited file
# that includes only the characterized ones.
# The fields are:
# accession, description, organism, sequence,
# and comments (may be multiple fields, each starting with TOPIC::)


use strict;
use lib "SWISS/lib";
use SWISS::Entry;

die "Run as a filter\n" unless @ARGV==0;

# Read an entire record at a time
local $/ = "\n//\n";

while(<>) {
    my $text = $_;
    my $entry = SWISS::Entry->fromText($text);
    my $desc = join ("; ", map { $_->text } $entry->DEs->elements);
    my $organisms = join ("; ", map { $_->text } $entry->OSs->elements);
    # Each comment is a CC object with fields topic and comment
    # References to pubmed normally show up as
    # {ECO:0000269|PubMed:18346083}
    # However sometimes the evidence code is in the special blocks section?
    # So instead do
    my @cc = map $_->toString(), $entry->CCs->elements;
    my @cc2 = map { s/^CC +//; s/\nCC +/ /g; s/\n+//g; s/^-!- *//; $_; } @cc;
    my $cc = join("_:::_", @cc2);
    # which produces multiple fields like FUNCTION: ....
    # joined by the arbitrary _:::_ join.

    # AC = accession
    # If there are multiple accessions (due to merging), then AC will report just the first (primary) one,
    # which desirable (I think)
    # SQ = sequence
    print join("\t", $entry->AC, $desc, $organisms, $entry->SQ, $cc)."\n"
        if $cc =~ m/ECO:0000269/;
}
