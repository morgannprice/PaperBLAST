#!/usr/bin/perl -w
# Given the hits from querying EuropePMC, look for any expected hits from the
# EuropePMC-provided links to UniProt and to PubMed. For these, create a new
# hits file. Calls pubmedFields.pl to get the metadata about the papers

use strict;
use FindBin qw{$Bin};
use Getopt::Long;

my $usage = <<END
addPMCLinks.pl -query comb.query -papers hits.papers -in ind -out extra_links

The query input file should be tab delimited with fields organism,
query term , protein identifier, protein sequence, and description (as
from oaquery.pl).

The papers file should be tab-delimited with fields queryId,
query_term, pmcId, pmId, doi, title, authors, journal, year, isOpen
(as from parseEuropePMCHits.pl; only the queryId and pmId fields are
used).

The input directory should contain comma-delimited tables
RefSeq_PMC.csv and UniProt_PMC.csv, with fields RefSeq or UniProt,
PMCID, and PMID.

Writes three files (out.queries, out.papers, and out.faa), in the same
format as parseEuropePMCHits.pl. These will include only the hits that
were not already in the hits.papers file. Any links for genes that
are not in the query input file are ignored.
END
;

my ($queryfile, $papersfile, $indir, $out, $test);
die $usage
  unless GetOptions('query=s' => \$queryfile,
                    'papers=s' => \$papersfile,
                    'in=s' => \$indir,
                    'out=s' => \$out,
                    'test' => \$test )
  && @ARGV == 0
  && defined $queryfile && defined $papersfile && defined $indir && defined $out;
die "No such file: $queryfile\n" unless -e $queryfile;
die "No such file: $papersfile\n" unless -e $papersfile;
die "No such directory: $indir\n" unless -d $indir;
my $inRefSeq = "$indir/RefSeq_PMC.csv";
die "No such file: $inRefSeq\n" unless -e $inRefSeq;
my $inUniProt = "$indir/UniProt_PMC.csv";
die "No such file: $inUniProt\n" unless -e $inUniProt;

my $pmFields = "$Bin/pubmedFields.pl";
die "No such executable: $pmFields\n" unless -x $pmFields;

print STDERR "Test mode -- not actually running pubmedFields.pl\n"
  if defined $test;

my %link = (); # queryId => pmId => 1 if in the papers file, 0 otherwise
my $nLink = 0;
my $nRefSeq = 0;
open(REFSEQ, "<", $inRefSeq) || die "Cannot read $inRefSeq";
while(my $line = <REFSEQ>) {
  chomp $line;
  my ($refseqId, $pmcId, $pmId) = split /,/, $line;
  die "Invalid input in $inRefSeq: $line" unless defined $pmId;
  next unless $pmId =~ m/^\d+$/; # occassionally missing; also skip header
  $refseqId =~ s/[.]\d+//;
  $link{$refseqId}{$pmId} = 0;
  $nRefSeq++;
  $nLink++;
}
close(REFSEQ) || die "Error reading $inRefSeq";
print STDERR "Read $nRefSeq useful entries from $inRefSeq\n";

my $nUniProt = 0;
open(UNIPROT, "<", $inUniProt) || die "Cannot read $inUniProt";
while(my $line = <UNIPROT>) {
  chomp $line;
  my ($uniprotId, $pmcId, $pmId) = split /,/, $line;
  next unless $pmId =~ m/^\d+$/;
  # UniProt ids sometimes show up with a version number, which we ignore
  # Also note that some of the UniProt ids are secondary accession numbers, which
  # may fail to match (unless sprotToQuery.pl is updated)
  $uniprotId =~ s/[.]\d+$//;
  $link{$uniprotId}{$pmId} = 0;
  $nLink++;
  $nUniProt++;
}
close(UNIPROT) || die "Error reading $inUniProt";
print STDERR "Read $nUniProt useful entries from $inUniProt\n";

# proteinId => list of queryId, organism, protein sequence, description
# For RefSeq, the query id will have the version number, but the query
# term will not. The mention in PMC often lacks the version number.
# So, we index by either queryId or queryTerm.
my %query = ();
open(QUERY, "<", "$queryfile") || die "Cannot read $queryfile";
while(my $line = <QUERY>) {
  chomp $line;
  my ($org, $queryTerm, $queryId, $seq, $desc) = split /\t/, $line;
  die "Not enough fields in $queryfile: $line" unless defined $desc;
  my $row = [ $queryId, $org, $seq, $desc ];
  # For RefSeq, the identifier 
  $query{$queryId} =  $row if exists $link{$queryId};
  $query{$queryTerm} = $row if exists $link{$queryTerm};
}
close(QUERY) || die "Error reading $queryfile";
print STDERR "Found " . scalar(keys %query) . " relevant queries in $queryfile, of " . scalar(keys %link) . "\n";

# Save the unknown queries to out.unknown
my $unknownFile = "$out.unknown";
open(UNKNOWN, ">", "$out.unknown") || die "Cannot write to $out.unknown";
foreach my $queryId (keys %link) {
  print UNKNOWN "$queryId\n" if !exists $query{$queryId};
}
close(UNKNOWN) || die "Error writing to $out.unknown";
print STDERR "Wrote $out.unknown\n";

# Identify the links that were already known
open(PAPERS, "<", $papersfile) || die "Cannot read $papersfile";
my $nLinkKnown;
while(my $line = <PAPERS>) {
  chomp $line;
  my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
  die "Invalid input in $papersfile: $line" unless defined $isOpen && $isOpen =~ m/^\d+/;
  next if $pmId eq "";
  if (exists $link{$queryId}{$pmId} || exists $link{$queryTerm}{$pmId}) {
    $link{$queryId}{$pmId} = 1  if exists $link{$queryId}{$pmId};
    $link{$queryTerm}{$pmId} = 1  if exists $link{$queryTerm}{$pmId};
    $nLinkKnown++;
  }
}
close(PAPERS) || die "Error reading $papersfile";
print STDERR "$nLinkKnown of $nLink PMC links are already known\n";

# Make a list of pmIds to fetch metadata for
my %pmId = ();
while (my ($queryTerm, $hash) = each %link) {
  next unless exists $query{$queryTerm};
  my ($queryId, $org, $seq, $desc) = @{ $query{$queryTerm} };
  while (my ($pmId, $known) = each %$hash) {
    $pmId{$pmId} = 1 unless $known;
  }
}
print STDERR "Identified " . scalar(keys %pmId) . " pubmed ids to fetch metadata for\n";
open(PMID, ">", "$out.pmid") || die "Cannot write to $out.pmid";
foreach my $pmId (sort keys %pmId) {
  print PMID $pmId."\n";
}
close(PMID) || die "Error writing to $out.pmid";
print STDERR "Wrote $out.pmid\n";

exit(0) if defined $test;

system("$pmFields < $out.pmid > $out.pmid.tab") == 0
  || die "Error running $pmFields: $!";

my %pm = (); # pmId to list of pmcId, DOI, title, author_list, journal_name, year
open(PM, "<", "$out.pmid.tab") || die "Cannot read $out.pmid.tab";
while(my $line = <PM>) {
  chomp $line;
  my @F = split /\t/, $line, -1; # -1 means do not drop empty last field (no year)
  my $pmId = shift @F;
  die "Unexpected pubmed id $pmId" unless exists $pmId{$pmId};
  $pm{$pmId} = \@F;
}
close(PM) || die "Error reading $out.pmid.tab";
print STDERR "Read metadata for " . scalar(keys %pm) . " pubmed ids\n";

# And now, make the output files:
# $out.queries has queryId, organism, protein_length, description
# $out.faa is in FASTA format with the queryId as the item name
# $out.papers has queryId, query_term, pmcId, pmId, doi, title, authors, journal, year, isOpen
my %queryWritten = (); # queryId => 1 if output already to query and FAA files
open(QUERIES, ">", "$out.queries") || die "Cannot write to $out.queries\n";
open(PAPERS, ">", "$out.papers") || die "Cannot write to $out.papers\n";
open(FAA, ">", "$out.faa") || die "Cannot write to $out.faa\n";

my $nLinksWritten = 0;
while (my ($queryTerm, $hash) = each %link) {
  next unless exists $query{$queryTerm};
  my ($queryId, $org, $seq, $desc) = @{ $query{$queryTerm} };
  while (my ($pmId, $known) = each %$hash) {
    if (!$known && exists $pm{$pmId}) {
      # print out the information
      if (!exists $queryWritten{$queryId}) {
        print FAA ">$queryId\n$seq\n";
        print QUERIES join("\t", $queryId, $org, length($seq), $desc)."\n";
        $queryWritten{$queryId} = 1;
      }
      my ($pmcId, $doi, $title, $authors, $journal, $year) = @{ $pm{$pmId} };
      die "Invalid metadata for pmId $pmId: no year" if !defined $year;
      my $isOpen = 1; # assume all PMC links are open
      print PAPERS join("\t",
                        $queryId, $queryTerm,
                        $pmcId, $pmId, $doi, $title, $authors, $journal, $year,
                        $isOpen)."\n";
      $nLinksWritten++;
    }
  }
}
close(QUERIES) || die "Error writing to $out.queries";
close(PAPERS) || die "Error writing to $out.papers";
close(FAA) || die "Error writing to $out.faa";

print STDERR "Wrote $nLinksWritten links and " . scalar(keys %queryWritten) . " previously uncovered proteins\n";
print STDERR "Success: wrote $out.{queries.papers,faa}\n";
