#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use LWP::Simple;
use XML::LibXML;

my $snippetBefore = 20;
my $snippetAfter = 35;
my $maxChar = 150;

my $usage = <<END
crossrefSnippets.pl -tokenfile token -papers epmc_pop.papers -dir cache
	The token file should contain a single line.
	Caches the fetched papers in cache/pubmed_*.xml
END
    ;

{
    my ($tokenfile, $papersfile, $cachedir);
    die $usage
        unless GetOptions('tokenfile=s' => \$token,
                          'papers=s' => \$papersfile,
                          'cache=s' => \$cachedir)
        && defined $tokenfile && defined $papersfile && defined $cachedir;

    die "No such file: $tokenfile\n" unless -e $tokenfile;
    die "No such file: $papersfile\n" unless -e $papersfile;
    die "No such directory: $cachedir\n" unless -d $cachedir;

    open(TOKEN, "<", $tokenfile) || die "Cannot read $tokenfile";
    my $token = <TOKEN>;
    chomp $token;
    die "Invalid token: $token" unless $token =~ m/^[a-zA-Z0-9-]+$/;
    close(TOKEN) || die "Error reading $tokenfile";

    my %fail = (); # pmId => 1 if failed to fetch it
    open(PAPERS, "<", $papersfile) || die "Cannot read $papersfile";
    while(my $line = <PAPERS>) {
        chomp $line;
        my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
        die "Invalid input in $papersfile" unless defined $isOpen;
        die "Invalid input in $papersfile: isOpen = $isOpen" unless $isOpen eq "1" || $isOpen eq "0";
        next unless $isOpen eq "0" && $pmId ne "" && $doi ne "";
        next if exists $fail{$pmId};
        my ($type, $URL) = &CrossrefToURL($doi);
        $fail{$pmId}  = 1 if ! $type;
        print join("\t", $journal, $doi, $type, $URL) . "\n" if $type;
    }
}
