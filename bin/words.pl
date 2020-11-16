#!/usr/bin/perl -w
# Given XML articles file(s) from PMC Europe,
# find all the locus tag like words
use strict;
use Getopt::Long;
use XML::LibXML;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner;

sub ProcessArticle($);

{
  my @files = ();
  GetOptions('in=s{1,}' => \@files)
    || die "Run as a filter, or with -in file1 ... flieN\n";
  my $n = 0;
  if (@files > 0) {
    print STDERR "Processing " . scalar(@files) . " files\n";
    foreach my $file (@files) {
      open(my $fh, "<", $file) || die "Error reading $file";
      my $state = {};
      while (my $dom = ReadXMLArticle($fh, $state)) {
        $n++;
        print STDERR "Parsed $n articles\n" if ($n % 1000) == 0;
        ProcessArticle($dom);
      }
    }
  } else {
    my $fh = *STDIN;
    my $state = {};
    while (my $dom = ReadXMLArticle($fh, $state)) {
      $n++;
      ProcessArticle($dom);
      print STDERR "Parsed $n articles\n" if ($n % 1000) == 0;
    }
    print STDERR "Processed $n articles\n";
  }
}

sub ProcessArticle($) {
    my ($dom) = @_;
    my $pmcid = &DomToPMCId($dom);
    $pmcid = "" if !defined $pmcid;

    # For starters, just print the text
    my $text = &RemoveWideCharacters( &NodeText($dom) );
    my @words = split /\s+/, $text;
    my %seen = ();
    foreach my $word (@words) {
      # Looking for potential locus tags with a letter, ending with at least 3 numbers in a row,
      # and removing any trailing punctuation
      # This pattern also works for UniProt accessions and for RefSeq protein ids (without the version), but
      # not for UniProt entry names, which generally start with digits or end with letters
      $word =~ s/[.,;-]+$//;
      $word =~ s/^\[//;
      $word =~ s/^\(//;
      $word =~ s/\]$//;
      $word =~ s/\)$//;
      # remove version numbers from refseq ids
      $word =~ s/[.]\d+$// if $word =~ m/^[A-Z][P]_\d+[.]\d+$/;
      next unless $word =~ m/^[a-zA-Z][a-zA-Z90-9_]+\d\d\d$/;
      print "$pmcid\t$word\n" unless exists $seen{$word};
      $seen{$word} = 1;
    }
}
