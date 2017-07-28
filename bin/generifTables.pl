#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<END
generifTables.pl -rif generifs_basic -prot generifs_prot
                 -query generif.queries -known hits.queries
                 -papers generifs_pmid.tab -out generif_tab

Given the query file for the protein ids in GeneRIF, the gene to
protein mapping, the query file with information about those proteins,
a set of known queries (to avoid duplicates), and metadata about the
relevant papers (from pubmedFields.pl), writes out the new queries,
the gene-paper mappings, and the snippets. These are written to:
out.queries out.papers out.faa out.rif

The first three files are in the same format as the output of
ParseEuropePMC.pl:

out.queries is tab-delimited with protein id, organism name, number of
amino acids, and description. Only proteins not in the known file are
included.

out.papers is tab-delimited with protein id, query term (short protein
id), pmcId, pmId, doi, title, authors, journal, year (and blank
instead of open access or not).

out.faa has the fasta sequences of the new proteins.

out.rif has the same format as a snippets file, with pmcId, pmid,
queryTerm (short protein id), protein id, and GeneRIF
END
  ;

my ($riffile, $protfile, $queryfile, $knownfile, $papersfile, $out);

die $usage
  unless GetOptions('rif=s' => \$riffile,
                    'prot=s' => \$protfile,
                    'query=s' => \$queryfile,
                    'known=s' => \$knownfile,
                    'papers=s' => \$papersfile,
                    'out=s' => \$out)
  && @ARGV == 0
  && defined $riffile && defined $protfile
  && defined $queryfile && defined $knownfile
  && defined $papersfile && defined $out;

foreach my $file ($riffile, $protfile, $queryfile, $knownfile, $papersfile) {
  die "No such file: $file\n" unless -e $file;
}

my %known = (); # queryId => length
open(KNOWN, "<", $knownfile) || die "Cannot read $knownfile";
while (my $line = <KNOWN>) {
  chomp $line;
  my ($org, $queryTerm, $queryId, $seq, $desc) = split /\t/, $line;
  die "Invalid query in $knownfile" unless defined $seq;
  $known{$queryId} = length($seq);
}
close(KNOWN) || die "Error reading $knownfile";

open(QUERIES, "<", $queryfile) || die "Cannot read $queryfile";
open(FAA, ">", "$out.faa") || die "Cannot write to $out.faa";
open(OUT, ">", "$out.queries") || die "Cannot write to $out.queries";
my $nNew = 0;
while(my $line = <QUERIES>) {
  chomp $line;
  my ($org, $queryTerm, $queryId, $seq, $desc) = split /\t/, $line;
  die "Invalid query in $queryfile" unless defined $desc;
  my $len = length($seq);
  if (exists $known{$queryId}) {
    print STDERR "Warning: inconsistent lengths for $queryId, $len vs. $known{$queryId}\n"
      if $known{$queryId} != $len;
    next;
  }
  $nNew++;
  $known{$queryId} = $len;
  print FAA ">$queryId\n$seq\n";
  print OUT join("\t", $queryId, $org, $len, $desc)."\n";
}
close(QUERIES) || die "Error reading $queryfile";
close(FAA) || die "Error writing to $out.faa";
close(OUT) || die "Error writing to $out.queries";
print STDERR "Wrote $out.queries ($nNew entries)\n";
print STDERR "Wrote $out.faa\n";

# map short protein ids (without the .version suffix) to protein ids
my %shortToProt = ();
foreach my $protId (keys %known) {
  my $short = $protId; $short =~ s/[.]\d+$//;
  $shortToProt{$short} = $protId;
}

my %gene = (); # gene id to protein id
open(PROT, "<", $protfile) || die "Error reading $protfile";
while(my $line = <PROT>) {
  chomp $line;
  my ($geneId, $protId) = split /\t/, $line;
  die "Invalid input in $protfile" unless defined $protId;
  die "Invalid gene id $geneId" unless $geneId =~ m/^\d+$/;
  $protId = $shortToProt{$protId} if exists $shortToProt{$protId};
  if (exists $gene{$geneId}) {
    my $protId2 = $gene{$geneId};
    if (exists $known{$protId}) {
      print STDERR "Warning: inconsistent lengths for $geneId -- $protId $known{$protId} vs. $protId2 $known{$protId2}\n"
        unless $known{$protId} == $known{$protId2};
    }
  } elsif (exists $known{$protId}) {
    $gene{$geneId} = $protId;
  } else {
    print STDERR "Warning: unknown protein $protId for gene $geneId\n";
  }
}
close(PROT) || die "Error reading $protfile";
print STDERR "Mapped " . scalar(keys %gene) . " gene ids to known protein ids\n";

my %papers = (); # pmId to list of pmcId, doi, title, authors, journal, year
open(PAPERS, "<", $papersfile) || die "Cannot read $papersfile";
while(my $line = <PAPERS>) {
  chomp $line;
  my ($pmId, $pmcId, $doi, $title, $authors, $journal, $year) = split /\t/, $line;
  if (!defined $year) {
    print STDERR "Warning: invalid paper from pubmed with ids $pmId $pmcId $doi\n";
    next;
  }
  $papers{$pmId} = [ $pmcId, $doi, $title, $authors, $journal, $year ];
}
close(PAPERS) || die "Error reading $papersfile";
print STDERR "Read metadata for " . scalar(keys %papers) . " pubmed ids\n";

my %known_prot_pmId = (); # if proteinId::pmId is set, suppress duplicate papers line
open(RIF, "<", $riffile) || die "Error reading $riffile";
open(PAPERS, ">", "$out.papers") || die "Cannot write to $out.papers";
open(OUT ,">", "$out.rif") || die "Cannot write to $out.rif";
while(my $line = <RIF>) {
  chomp $line;
  my ($taxId, $geneId, $pmIds, $timestamp, $rif) = split /\t/, $line;
  next if $taxId =~ m/^#/; # comment
  die "Invalid rif line" unless defined $rif;
  next if $rif eq "";
  next unless exists $gene{$geneId};
  my $protId = $gene{$geneId};
  my $shortId = $protId; $shortId =~ s/[.]\d+$//;
  $pmIds =~ m/^[0-9,]+$/ || die "Invalid pubmed specifier $pmIds";
  my @pmIds = split /,/, $pmIds;
  my $pmId = $pmIds[0]; # arbitrarily pick the first one (almost always is just one)
  if (!exists $papers{$pmId}) {
    print STDERR "Warning: unknown pubmed id $pmId, skipping it\n";
    next;
  }
  my ($pmcId, $doi, $title, $authors, $journal, $year) = @{ $papers{$pmId} };
  my $key = join("::", $protId, $pmId);
  if (!exists $known_prot_pmId{$key}) {
    $known_prot_pmId{$key} = 1;
    print PAPERS join("\t", $protId, $shortId, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, "")."\n";
  }
  print OUT join("\t", $pmcId, $pmId, $shortId, $protId, $rif)."\n";
}

