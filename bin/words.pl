#!/usr/bin/perl -w
# Given XML articles file(s) from PMC Europe,
# find all the locus tag like words
use strict;
use XML::LibXML;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner;

sub ProcessArticle($);

{
    die "Run as a filter\n" unless @ARGV==0;
    my $firstline = <STDIN>;
    chomp $firstline;

    my $n = 0;
    if ($firstline =~ m/^<!DOCTYPE/) {
        # Handle a single article starting with a <!DOCTYPE ...> tag followed by
        # an <article> tag, as from the author's manuscripts part of EuropePMC
        # Also handles a series of such files concatenated together
        my @lines = ($firstline);
        while(my $line = <STDIN>) {
            chomp $line;
            push @lines, $line;
            if ($line eq "</article>") {
                my $dom = XML::LibXML->load_xml(string => join("\n",@lines), recover => 1);
                $n++;
                print STDERR "Parsed $n articles\n" if ($n % 1000) == 0;
                ProcessArticle($dom);
                @lines = ();
            }
        }
    } elsif ($firstline eq "<articles>") {
        # Handle many articles, each with an <article> tag, with an <articles> tag above them,
        # as from the open access publications part of EuropePMC
        while(my $line = <STDIN>) {
            chomp $line;
            next if $line eq "";
            if ($line eq "</articles>") {
                next;
                last;
            }
            $line =~ m/^<article[> ]/ || die "Does not begin with an article: " . substr($line, 0, 100);
            my @lines = ( $line );
            while ($lines[-1] !~ m!</article>!) {
                $line = <STDIN>;
                if (defined $line) {
                    chomp $line;
                    push @lines, $line;
                } else {
                    last;
                }
            }
            if ($lines[-1] !~ m!</article>!) {
                print STDERR "Last article is truncated? Starts with " . substr($lines[0], 0, 100);
            }
            # without recover, get errors like
            # :1: parser error : Premature end of data in tag p line 1
            # l 2004 news article &#x0201c;Reaching across the Border with the SBRP&#x0201d; [

            my $dom = XML::LibXML->load_xml(string => join("\n",@lines), recover => 1);
            $n++;
            print STDERR "Parsed $n articles\n" if ($n % 1000) == 0;
            ProcessArticle($dom);
        }
    } else {
        die "First line should be either be\n<articles>\nor should be a <!DOCTYPE indicator\n";
    }
    print STDERR "Processed $n articles\n";
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
