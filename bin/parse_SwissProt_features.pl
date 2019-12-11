#!/usr/bin/perl -w
# Extract evidence-based features from the SwissProt file
use strict;
use lib "SWISS/lib";
use SWISS::Entry;

# If an entry has experimental evidence for only these ignored types,
#   ignore it.
# Allow NP_BIND because usually describes just a few residues
# (unlike DNA_BIND or ZN_FING)
my @ignoreTypes = qw{REGION VARIANT UNSURE CONFLICT VAR_SEQ
                     INIT_MET SIGNAL TRANSIT
                     TOPO_DOM TRANSMEM INTRAMEM COILED
                     PROPEP PEPTIDE CHAIN
                     ZN_FING DNA_BIND DOMAIN
                     NON_TER NON_CONS COMPBIAS};
# (Not sure why COMPBIAS would have experimental evidence, but it doesn't)
my %ignoreTypes = map { $_=> 1 } @ignoreTypes;

my $usage = <<END
Run as a filter on swiss-prot entries.
Output is tab-delimited with fields
  accession, id, description,
  feature type, feature begin, feature end, feature comment

Only features with experimental evidence (ECO:0000269) are included.

SwissProt entries whose only evidence-based features are of the
following types are ignored:
  @ignoreTypes
END
;

die $usage unless @ARGV == 0;

# Read an entire record at a time
local $/ = "\n//\n";

while(<>) {
  my $text = $_;
  my $entry = SWISS::Entry->fromText($text);

  # each entry is an array of: feature type or "key", from_position, to_position,
  # description, qualifier, FTid, evidence tag

  my @ftEvidence = grep { $_->[6] =~ m/ECO:0000269/ } $entry->FTs->elements;
  next unless @ftEvidence > 0;
  my @ftEvidence2 = grep { !exists $ignoreTypes{$_->[0]} } @ftEvidence;
  next unless @ftEvidence2 > 0;

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
  push @desc, @ec;
  my $desc = join("; ", @desc);
  foreach my $ft (@ftEvidence) {
    my ($ftKey, $ftFrom, $ftTo, $ftDesc) = @$ft;
    print join("\t",
               $entry->AC,
               @{ $entry->IDs->list }[0] || "",
               $desc,
               $ftKey, $ftFrom, $ftTo, $ftDesc)."\n";
  }
}
