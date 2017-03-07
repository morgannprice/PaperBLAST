#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use Miner;

my $usage = <<END
Usage: run as a filter, with the pubmedid as the first column.

Output files are:
pubmedId, pmcId, DOI, title, author_list, journal_name, year
END
;

die $usage unless @ARGV == 0;

my @pmid = ();
while(<STDIN>) {
  chomp;
  my @F = split /\t/, $_;
  my $pmid = shift @F;
  next unless defined $pmid && $pmid =~ m/^\d+$/;
  push @pmid, $pmid;
  if (@pmid >= 100) {
    &Handle(@pmid);
    @pmid = ();
  }
}
&Handle(@pmid);

sub Handle {
  my @pmid = @_;
  my @dom = Miner::FetchPubMed(@pmid);
  foreach my $dom (@dom) {
    my $pmid = $dom->find("MedlineCitation/PMID");
    if (! $pmid) {
      print STDERR "Warning, no pmid found in results for $pmid[0] to $pmid[-1]\n";
      next;
    }
    my @ids = $dom->findnodes("PubmedData/ArticleIdList/ArticleId");
    my %ids = ();
    foreach my $node (@ids) {
      my $idType = $node->attributes()->getNamedItem("IdType")->textContent();
      $ids{$idType} = $node->textContent();
    }
    my $pmcid = $ids{"pmc"} || "";
    my $doi = $ids{"doi"} || "";
    my $journal = $dom->find("MedlineCitation/Article/Journal/Title") || "";
    my $title = $dom->find("MedlineCitation/Article/ArticleTitle") || "";
    my $year = $dom->find("MedlineCitation/Article/Journal/JournalIssue/PubDate/Year") || "";
    my @authors = $dom->findnodes("MedlineCitation/Article/AuthorList/Author");
    my @authStrings = ();
    foreach my $author (@authors) {
      push @authStrings, $author->find("LastName") . " " . $author->find("Initials");
    }
    my $line = join("\t", $pmid, $pmcid, $doi, $title, join(", ", @authStrings), $journal, $year)."\n";
    utf8::encode($line);
    print $line;
  }
}
