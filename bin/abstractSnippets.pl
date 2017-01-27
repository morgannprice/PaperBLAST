#!/usr/bin/perl -w

use strict;
use XML::LibXML;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner; # for TextSnippets()

my $snippetBefore = 15;
my $snippetAfter = 30;
my $maxChar = 1000;

my $usage = <<END
abstractSnippets.pl -in abstracts -out pubmed.snippets oa.papers ... refseq.papers

	The input file should be tab-delimited with the 1st column
	having the pubmedId and the second column having the abstract
	(as from pubmedparse.pl)

	The papers files are from parseEuropePMCHits.pl

	Output will be tab-delimited with the fields
	pmcId, pmId, queryTerm, queryId, snippet

	Optional arguments control the size of each snippet in words or characters:
	-before $snippetBefore
        -after $snippetAfter
        -maxChar $maxChar

	There is a limit of 2 snippets per gene per paper
END
    ;

sub ProcessArticle($$$);

{
    my ($infile, $listfile, $outfile);
    my @paperfiles;
    die $usage
        unless GetOptions( 'in=s' => \$infile,
                           'list=s' => \$listfile,
                           'out=s' => \$outfile,
                           'before=i' => \$snippetBefore,
                           'after=i' => \$snippetAfter,
                           'maxChar=i' => \$maxChar )
        && (defined $infile xor defined $listfile)
        && defined $outfile
        && scalar(@ARGV) > 0;
    @paperfiles = @ARGV;
    die "-before cannot be negative" if $snippetBefore < 0;
    die "-after must be positive" if $snippetAfter < 1;
    die "-maxChar must be at least 50" if $maxChar < 50;

    my %papers = (); # pmId => queryTerm => list of queryId
    my %pm2pmc = (); # pmId => pmcId
    foreach my $paperfile (@paperfiles) {
        die "No such file: $paperfile\n" unless -e $paperfile;
    }
    foreach my $paperfile (@paperfiles) {
        open(PAPERS, "<", $paperfile) || die "Error reading $paperfile";
        while(my $line = <PAPERS>) {
            chomp $line;
            my ($queryId, $queryTerm, $pmcId, $pmId) = split /\t/, $line;
            die "Not enough fields in line\n$line\nin $paperfile" unless defined $pmId;
            next if $pmId eq "";
            push @{ $papers{$pmId}{$queryTerm} }, $queryId;
            $pm2pmc{$pmId} = $pmcId;
        }
        close(PAPERS) || die "Error reading $paperfile";
    }
    print STDERR "Read " . scalar(keys %papers) . " pubmedIds to search\n";

    my $nRelevant = 0;
    open(IN, "<", $infile) || die "Cannot read $infile";
    open(OUT, ">", $outfile) || die "Error writing to $outfile";
    while(my $line = <IN>) {
        chomp $line;
        my ($pmId, $abstract) = split /\t/, $line;
        next unless $abstract;
        next unless exists $papers{$pmId};
        $nRelevant++;
        my $hash = $papers{$pmId};
        my $pmcId = $pm2pmc{$pmId};
        while (my ($queryTerm, $queryIds) = each %$hash) {
            next unless $abstract =~ m/$queryTerm/;
            my @snippets = TextSnippets($abstract, $queryTerm, $snippetBefore, $snippetAfter, $maxChar);
            my %queryIds = map { $_ => 1 } @$queryIds;
            my @queryIds = keys(%queryIds);
            foreach my $queryId (@queryIds) {
                foreach my $snippet (@snippets) {
                    print OUT join("\t", $pmcId, $pmId, $queryTerm, $queryId, $snippet)."\n";
                }
            }
        }
    }
    print STDERR "Processed $nRelevant relevant abstracts\n";
}

