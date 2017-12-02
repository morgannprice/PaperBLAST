#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ParsePTools, ReadFastaEntry

my $usage = <<END
Usage: parse_metacyc.pl [ -debug ] metacyc_data_directory uniprot_fasta trembl_fasta > metacyc.curated_parsed

The data directory must contain proteins.dat and enzrxn.dat
in attribute-value format.

The uniprot and trembl fasta filenames can be .gz files instead.

The output file is in curated_parsed format, with fields metacyc,
metacyc identifier, uniprot identifier, short name, description, organism,
sequence, blank (comment), pubmed ids (comma separated).

END
;
my $debug = 0;
if ($ARGV[0] eq "-debug") {
  shift @ARGV;
  $debug = 1;
}
die $usage unless @ARGV == 3;
my ($metadir, $faa1, $faa2) = @ARGV;
die "No such directory: $metadir\n" unless -d $metadir;
die "No such file: $faa1\n" unless -e $faa1;
die "No such file: $faa2\n" unless -e $faa2;

# A list of hashes, one per protein. Each entry has fields
# UNIQUE-ID, ABBREV-NAME, COMMON-NAME, CATALYZES, UNIPROT, CITATIONS (as a list of pmIds),
# and ultimately SEQUENCE
my @prot = ();

# Find each protein that has a link to sequence (a uniprot id), links to paper(s),
# and either a useful description ("COMMON-NAME") or a link to an enzymatic reaction
# (where we can get a description).

my $pfile = "$metadir/proteins.dat";
my $nProt = 0;
open(my $fhProt, "<", $pfile) || die "Cannot read $pfile\n";
while (my $prot = ParsePTools($fhProt)) {
  $nProt++;

  # Note that these objects contain SPECIES attributes of the form
  # TAX-nnnn, where the number is an NCBI taxonomy id, but that is
  # ignored in favor of the OS field in the uniprot fasta file.

  # Try to link to sequence
  my @dblinks = map { $_->{"value"} } @{ $prot->{"DBLINKS"} };
  my @uniprot = grep m/UNIPROT/, @dblinks;
  next unless @uniprot > 0; # skip if no sequence
  my $uniprotLine = $uniprot[0];
  my $uniprotId;
  if ($uniprotLine =~ m/[(]UNIPROT "([0-9A-Z_]+)"/) {
    $uniprotId = $1;
  } else {
    print STDERR "Cannot parse uniprot id from DBLINKS entry: $uniprotLine\n";
    next;
  }

  # Links to papers may be in the CITATIONS attribute or in the CITATIONS annotation of GO-TERMS
  # or in the COMMENT field as entries like |CITS: [8621577]| or |CITS: [2693216][1425658]| or |CITS: [9755155] [9388228]|
  my @pmLines = ();
  foreach my $value ( @{ $prot->{"CITATIONS"} }) {
    push @pmLines, $value->{"value"};
  }
  foreach my $value ( @{ $prot->{"GO_TERMS"} }) {
    push @pmLines, $value->{CITATIONS}
      if exists  $value->{CITATIONS};
  }
  my %pmIds = ();
  foreach my $line (@pmLines) {
    my $pmId = $1
      if $line =~ m/^(\d+):/ || $line =~ m/^(\d+)$/ || $line =~ m/^\[(\d+)\]$/;
    # Some citations are like Wang02, not sure why
    $pmIds{$pmId} = 1 if $pmId;
  }
  my $comment = $prot->{"COMMENT"}[0]{"value"} || "";
  while ($comment =~ m/[|]CITS: ([\[\]0-9 ]+)[|]/g) {
    # this should extract a field like [8801422] or [10501935][11872485] or [6706930] [7142155]
    my $ids = $1;
    $ids =~ s/[\[\]]/ /g;
    my @ids = split / +/, $ids;
    foreach my $id (@ids) {
      $pmIds{$id} = 1 if $id;
    }
  }
  my @pmIds = sort { $a <=> $b } keys %pmIds;
  next unless @pmIds > 0;

  my $save = { 'CITATIONS' => \@pmIds, 'UNIPROT' => $uniprotId };
  foreach my $key (qw/UNIQUE-ID ABBREV-NAME COMMON-NAME CATALYZES/) {
    $save->{$key} = $prot->{$key}[0]{"value"} || "";
  }
  # Only a handful of proteins lack a COMMON-NAME *and* are not associated with an enzymatic reaction.
  # Ignore those.
  push @prot, $save if $save->{"COMMON-NAME"} || $save->{"CATALYZES"};
  last if @prot >= 100 && $debug;
}
close($fhProt) || die "Error reading $pfile\n";
print STDERR "Read $nProt protein entries and saved " . scalar(@prot) . "\n";

my %enzrxn = (); # unique id => common name

# Parse the enzymatic reactions file to get descriptions for the catalyzed reactions
my $enzfile = "$metadir/enzrxns.dat";
open(my $fhEnz, "<", $enzfile) || die "Cannot read $enzfile";
while (my $enzrxn = ParsePTools($fhEnz)) {
  my $id = $enzrxn->{"UNIQUE-ID"}[0]{"value"};
  my $desc = $enzrxn->{"COMMON-NAME"}[0]{"value"};
  $enzrxn{$id} = $desc if $id && $desc;
}
close($fhEnz) || die "Error reading $fhEnz";
print STDERR "Parsed descriptions for " . scalar(keys %enzrxn) . " enzymatic reactions\n";

foreach my $prot (@prot) {
  if ($prot->{"COMMON-NAME"} eq "" && $prot->{"CATALYZES"}) {
    my $enzrxnId = $prot->{"CATALYZES"};
    $prot->{"COMMON-NAME"} = $enzrxn{$enzrxnId} if exists$enzrxn{$enzrxnId};
  }
}
@prot = grep { $_->{"COMMON-NAME"} } @prot;
print STDERR "Reduced to " . scalar(@prot) . " proteins\n";

# Which sequences we are looking for; empty string for as-yet unknown sequence
my %uniprot = map { $_->{UNIPROT} => "" } @prot;
my %uniprotOrg = ();

# And now fetch sequences from the fasta files
foreach my $faafile ($faa1, $faa2) {
  my $fh;
  if ($faafile =~ m/[.]gz$/) {
    open($fh, "zcat $faafile |") || die "Cannot run zcat on $faafile";
  } else {
    open($fh, "<", $faafile);
  }
  my $state = {};
  while(my ($header,$seq) = ReadFastaEntry($fh,$state)) {
    my ($ids) = split / /, $header;
    my @ids = split /[|]/, $ids;
    shift @ids; # ignore sp or tr at front
    die "Invalid header $header" if @ids == 0;
    foreach my $id (@ids) {
      if (exists $uniprot{$id}) {
        $uniprot{$id} = $seq;
        if ($header =~ m/ OS=([^=]+) /) {
          $uniprotOrg{$id} = $1;
        }
      }
    }
  }
}

my $nWritten = 0;
foreach my $prot (@prot) {
  my $uniprotId = $prot->{"UNIPROT"};
  next unless exists $uniprot{$uniprotId};
  my $seq = $uniprot{$uniprotId};
  my $org = $uniprotOrg{$uniprotId} || "";
  if ($seq eq "") {
    print STDERR "Warning: no sequence found for $uniprotId\n";
  } else {
    print join("\t",
               "metacyc", $prot->{"UNIQUE-ID"}, $uniprotId,
               $prot->{"ABBREV-NAME"}, $prot->{"COMMON-NAME"}, $org,
               $seq, "", join(",", @{ $prot->{"CITATIONS"} })
              )."\n";
  }
}
