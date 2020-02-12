#!/usr/bin/perl -w
# Extract evidence-based features from the SwissProt file
use strict;
use lib "SWISS/lib";
use SWISS::Entry;
use SWISS::Ref;
use SWISS::FTs;
use SWISS::TextFunc;

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
  accession, id, description, organism
  feature type, feature begin, feature end, feature comment,
  pubmedIds (comma delimited)

Only features with experimental evidence (ECO:0000269) are included.

SwissProt entries whose only evidence-based features are of the
following types are ignored:
  @ignoreTypes
END
;
my $debug;
if (@ARGV > 0 && $ARGV[0] eq "-debug") {
  $debug = 1;
  shift @ARGV;
}
die $usage unless @ARGV == 0;

# Read an entire record at a time
local $/ = "\n//\n";

while(<>) {
  my $text = $_;
  my $entry = SWISS::Entry->fromText($text);

  # each entry is an array of: feature type or "key", from_position, to_position,
  # description, qualifier, FTid, evidence tag

  my $fts = $entry->FTs->list;
  # Older UniProt files had evidence information in [6] vs. [7]
  my @ftEvidence = grep { $_->[6] =~ m/ECO:0000269/ || $_->[7] =~ m/ECO:0000269/ } @$fts;
  print join("\t", $entry->AC, scalar(@$fts), scalar(@ftEvidence))."\n" if $debug;
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
    # Swiss-Prot files from 2019, parsed with Swisskinfe 1.73, had the
    # evidence code in [6] ($ev)
    # As of February 2020, need to use Swisskinfe 1.79 to parse the
    # evidence code, and it is in [7] ($ev2)
    my ($ftKey, $ftFrom, $ftTo, $ftDesc, undef, undef, $ev, $ev2) = @$ft;
    my $evidence = $ev2 || $ev;
    my @pmIds = $evidence =~ m!ECO:0000269[|]PubMed:(\d+)!g;
    my %pmId = map { $_ => 1 } @pmIds;
    @pmIds = sort {$a <=> $b} keys %pmId;
    print join("\t",
               $entry->AC,
               @{ $entry->IDs->list }[0] || "",
               $desc,
               join ("; ", map { $_->text } $entry->OSs->elements),
               $ftKey, $ftFrom, $ftTo, $ftDesc,
              join(",",@pmIds))."\n";
  }
}
