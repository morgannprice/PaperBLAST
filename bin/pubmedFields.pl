#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use Getopt::Long;
use Miner;

my $usage = <<END
pubmedFields.pl < pubmed_ids > pubmed_ids.fields
or
pubmedFields.pl -xml input.xml < pubmed_ids > pubmed_ids.fields
or
pubmedFields.pl -gz input.xml.gz < pubmed_ids > pubmed_ids.fields

If no xml or gz file is specified, sends queries to NCBI's entrez
(which is a bit slow).

Standard input should contain one pubmed id per line. Any additional
fields (tab-delimited) are ignored.

Output fields are:
pubmedId, pmcId, DOI, title, author_list, journal_name, year
END
;

sub DomToPubmedId($);
sub PrintDom($);

my ($xmlfile, $gzfile) = @_;
die $usage unless
  GetOptions('xml=s' => \$xmlfile,
             'gz=s' => \$gzfile)
  && @ARGV == 0;
die "Cannot specify both -xml and -gz\n" if defined $xmlfile && defined $gzfile;
die "No such file: $xmlfile\n" if defined $xmlfile && ! -e $xmlfile;
die "No such file: $gzfile\n" if defined $gzfile && ! -e $gzfile;

# read in the list of pubmed ids
my @pmid = ();
while(<STDIN>) {
  chomp;
  my @F = split /\t/, $_;
  my $pmid = shift @F;
  next unless defined $pmid && $pmid =~ m/^\d+$/;
  push @pmid, $pmid;
}
print STDERR "Looking for " . scalar(@pmid) . " pubmed ids\n";

if (defined $xmlfile || defined $gzfile) {
  if (defined $xmlfile) {
    open(IN, "<", $xmlfile) || die "Cannot read $xmlfile\n";
  } elsif (defined $gzfile) {
    open(IN, "zcat $gzfile |") || die "Cannot zcat $gzfile: $!\n";
  }
  my %pmid = map { $_ => 1 } @pmid;
  my $nArticles = 0;
  my $nRelevant = 0;
  my @lines = ();
  while (my $line = <IN>) {
    chomp $line;
    next if $line =~ m/^<[?]xml.*>$/
      || $line =~ m/^<!DOCTYPE .*>$/
        ||  $line =~ m/^<PubmedArticleSet>$/;
    push @lines, $line;
    if ($line =~ m!</PubmedArticle>$!) {
      $lines[0] =~ s/^ +//; # parser does not like leading whitespace
      my $article = XML::LibXML->load_xml(string => join("\n",@lines), recover => 1);
      my @children = $article->childNodes;
      $article = $children[0] if @children == 1 && $children[0]->nodeName eq "PubmedArticle";
      @lines = ();
      $nArticles++;
      my $pmId = DomToPubmedId($article);
      die "No pubmed id" unless $pmId;
      if (exists $pmid{$pmId}) {
        $nRelevant++;
        PrintDom($article);
      }
    }
  }
  if (scalar(@lines) > 0) {
    if ($lines[0] =~ m/<DeleteCitation>/ || $lines[-1] eq "</PubmedArticleSet>") {
      ; # ok, ignoring trailer
    } else {
      print STDERR "Warning: ignored " . scalar(@lines) . " lines at end\n";
    }
  }
  close(IN) || die "Error reading " . ($xmlfile || $gzfile);
  print STDERR "Processed $nRelevant articles from $nArticles articles\n";
} else {
  while(@pmid > 0) {
    my @work = ();
    while (@work < 100 && @pmid > 0) {
      push @work, shift @pmid;
    }
    my @dom = Miner::FetchPubMed(@work);
    foreach my $dom (@dom) {
      my $pmid = DomToPubmedId($dom);
      if (! $pmid) {
        print STDERR "Warning, no pmid found in results for $work[0] to $work[-1]\n";
        next;
      }
      &PrintDom($dom);
    }
  }
}

sub PrintDom($) {
  my ($dom) = @_;
  my $pmid = DomToPubmedId($dom);
  die unless $pmid;
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
  # remove leading and trailing white space from title (occasionally get new lines)
  $title =~ s/^\s+//;
  $title =~ s/\s+$//;
  my $line = join("\t", $pmid, $pmcid, $doi, $title, join(", ", @authStrings), $journal, $year);
  $line =~ s/\n//g; # just in case some sneak through
  utf8::encode($line);
  print $line."\n";
}

sub DomToPubmedId($) {
  my ($dom) = @_;
  return $dom->find("MedlineCitation/PMID");
}
