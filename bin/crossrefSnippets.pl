#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner;

my $snippetBefore = 20;
my $snippetAfter = 35;
my $maxChar = 150;

my $usage = <<END
crossrefSnippets.pl -tokenfile token -papers epmc_pop.papers -dir cache
	The token file should contain a single line.
	Caches the fetched papers in cache/pubmed_*.xml

	Output will be tab-delimited with the fields
	pmcId, pmId, queryTerm, queryId, snippet

Optional arguments:
	-cache-only -- ignore any non-cached items (for testing)
END
    ;

{
    my ($tokenfile, $papersfile, $cachedir);
    my $cache_only;
    die $usage
        unless GetOptions('tokenfile=s' => \$tokenfile,
                          'papers=s' => \$papersfile,
                          'dir=s' => \$cachedir,
                          'cache-only' => \$cache_only )
        && defined $tokenfile && defined $papersfile && defined $cachedir;

    die "No such file: $tokenfile\n" unless -e $tokenfile;
    die "No such file: $papersfile\n" unless -e $papersfile;
    die "No such directory: $cachedir\n" unless -d $cachedir;
    print STDERR "Only processing cached files\n" if defined $cache_only;

    open(TOKEN, "<", $tokenfile) || die "Cannot read $tokenfile";
    my $token = <TOKEN>;
    chomp $token;
    die "Invalid token: $token" unless $token =~ m/^[a-zA-Z0-9_-]+$/;
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
        my ($type, $URL);
        if (defined $cache_only) {
            $type = "other"; # and check for other cases
        } else {
            ($type, $URL) = &CrossrefToURL($doi, $journal);
            if (! $type) {
                $fail{$pmId}  = 1;
                next;
            }
        }
        my $prefix = "$cachedir/pubmed_$pmId";
        my $file = "$prefix.$type";
        my $file2 = "$prefix.xml";
        my $file3 = "$prefix.pdf";
        if (-e $file) {
            print STDERR "Using cached $file for $doi $journal\n";
        } elsif (-e $file2) {
            print STDERR "Using cached $file2 for $doi $journal\n";
            $file = $file2;
            $type = "xml";
        } elsif (-e $file3) {
            print STDERR "Using cached $file3 for $doi $journal\n";
            $file = $file3;
            $type = "pdf";
        } elsif (defined $cache_only) {
            print STDERR "Skipping non-cached $doi $journal\n";
            next;
        } else {
            my $newtype = &CrossrefFetch($URL, $type, $file, $journal);
            if ($newtype) {
                # Figure out what type it is, and rename if necessary
                if ($newtype eq "application/xml") {
                    $newtype = "xml";
                } elsif ($newtype eq "application/pdf") {
                    $newtype = "pdf";
                } elsif ($newtype eq "text/plain") {
                    $newtype = "other";
                } else {
                    print STDERR "Warning: unrecognized content type $newtype\n";
                    $newtype = "other";
                }
                my $newfile = "$prefix.$newtype";
                if ($newfile ne $file) {
                    rename($file,$newfile) || die "Error renaming $file to $newfile";
                    $file = $newfile;
                }
                if ($type eq "pdf") {
                    my $filereport = `file $newfile`;
                    if (! $filereport =~ m/PDF document/) {
                        print STDERR "Error: fetched $newfile for $doi but it is not in PDF format\n";
                        unlink($file);
                        $fail{$pmId} = 1;
                        next;
                    }
                }
                print STDERR "Fetched $doi into $file from $journal\n";
            } else {
                print STDERR "Failed to fetch $doi from $journal\n";
                $fail{$pmId} = 1;
                next;
            }
        }
        my $text;
        if ($type eq "other") {
            $text = TextFileToText($file);
        } elsif ($type eq "xml") {
            $text = XMLFileToText($file);
        } elsif ($type eq "pdf") {
            $text = PDFFileToText($file);
        }
        if (! $text) {
            print STDERR "Error: failed to fetch data from type $type file $file\n";
            next;
        }
        my @snippets = TextSnippets($text, $queryTerm, $snippetBefore, $snippetAfter, $maxChar);
        foreach my $snippet (@snippets) {
            print join("\t", $pmcId, $pmId, $queryTerm, $queryId, $snippet)."\n";
        }
    }
}
