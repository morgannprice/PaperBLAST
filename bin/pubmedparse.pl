#!/usr/bin/perl -w
use strict;
use XML::LibXML;

die "Run as a filter\n" if @ARGV > 0;
binmode(STDOUT,":utf8");

my @lines = ();
my $nArticles = 0;
while (my $line = <STDIN>) {
    chomp $line;
    next if $line =~ m/^<[?]xml.*>$/
        || $line =~ m/^<!DOCTYPE .*>$/
        ||  $line =~ m/^<PubmedArticleSet>$/;
    push @lines, $line;
    if ($line =~ m!</PubmedArticle>$!) {
        $lines[0] =~ s/^ +//; # parser does not like leading whitespace
        my $article = XML::LibXML->load_xml(string => join("\n",@lines), recover => 1);
        @lines = ();
        $nArticles++;

        my @pmidNodes = $article->findnodes("//MedlineCitation/PMID");
        my $pmId = "";
        if (@pmidNodes == 1) {
            $pmId = $pmidNodes[0]->textContent;
        }
        # There can be multiple sections in the abstract
        my @abstractNodes = $article->findnodes("//MedlineCitation/Article/Abstract/AbstractText");
        my $abstract = join(" ", map { $_->textContent } @abstractNodes);
        # occasionaly newline characters are within the text nodes
        $abstract =~ s/\s/ /g;
        print join("\t", $pmId, $abstract)."\n" if $abstract ne "";
    }
}
if (scalar(@lines) > 0) {
    if ($lines[0] =~ m/<DeleteCitation>/ || $lines[-1] eq "</PubmedArticleSet>") {
        ; # ok, ignoring trailer
    } else {
        print STDERR "Warning: ignored " . scalar(@lines) . " lines at end\n";
    }
}
print STDERR "Processed $nArticles articles\n";
