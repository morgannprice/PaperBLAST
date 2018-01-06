#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner;

my $snippetBefore = 20;
my $snippetAfter = 35;
my $maxChar = 160;

my $usage = <<END
crossrefSnippets.pl -tokenfile token -papers epmc_pop.papers -dir cache -out snippets
	The token file should contain a single line.
	Caches the fetched papers in cache/pubmed_*.xml

	Output will be tab-delimited with the fields
	pmcId, pmId, queryTerm, queryId, snippet

	Also writes to output.access, which is tab-delimited with the fields
	pmcId, pmdId, "full"
	if it was able to read the paper.

Optional arguments:
	-cache-only -- ignore any non-cached items (for testing)
END
    ;

{
    my ($tokenfile, $papersfile, $cachedir, $outfile);
    my $cache_only;
    die $usage
        unless GetOptions('tokenfile=s' => \$tokenfile,
                          'papers=s' => \$papersfile,
                          'dir=s' => \$cachedir,
                          'out=s' => \$outfile,
                          'cache-only' => \$cache_only )
        && @ARGV == 0;
    die "Required argument is missing:\n$usage"
        unless defined $tokenfile && defined $papersfile && defined $cachedir && defined $outfile;

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
    open(OUT, ">", $outfile) || die "Error writing to $outfile";
    my $accfile = "$outfile.access";
    open(ACC, ">", $accfile) || die "Error writing to $accfile";
    # if the pmId is the same as the last line, save time by reusing the resulting text
    my $last_pmId = undef;
    my $last_text = undef;
    while(my $line = <PAPERS>) {
        chomp $line;
        my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
        die "Invalid input in $papersfile" unless defined $isOpen;
        die "Invalid input in $papersfile: isOpen = $isOpen" unless $isOpen eq "1" || $isOpen eq "0";
        next unless $isOpen eq "0" && $pmId ne "" && $doi ne "";
        next if exists $fail{$pmId};
        $queryTerm =~ s/^"//;
        $queryTerm =~ s/"$//;
        my $text;
        if (defined $last_pmId && $pmId eq $last_pmId) {
            $text = $last_text;
        } else {
            $last_pmId = undef;
            $last_text = undef;
            my $prefix = "$cachedir/pubmed_$pmId";
            my ($type, $file, $URL);
            if (-e "$prefix.pdf") {
              $type = "pdf";
              $file = "$prefix.pdf";
            } elsif (-e "$prefix.xml") {
              $type = "xml";
              $file = "$prefix.xml";
            } elsif (-e "$prefix.other") {
              $type = "other"; # text
              $file = "$prefix.other";
            }
            ($type, $URL) = &CrossrefToURL($doi, $journal) unless $file;
            if (! $type) {
              $fail{$pmId}  = 1;
              next;
            }
            if (defined $cache_only && ! $file) {
              print STDERR "Skipping non-cached $doi $journal\n";
              next;
            } elsif ($URL) {
              my $tmpfile = "$prefix.tmp";
              my $newtype = &CrossrefFetch($URL, $type, $tmpfile, $journal);
              if ($newtype) {
                # Figure out what type it is, and rename if necessary
                if ($newtype eq "application/xml") {
                  $newtype = "xml";
                } elsif ($newtype eq "application/pdf") {
                  $newtype = "pdf";
                } elsif ($newtype eq "application/octet-stream") {
                  # Download in another way -- not sure why, but these never worked.
                  unlink($tmpfile);
                  if (system("wget","-nv","-O",$tmpfile,$URL) == 0) {
                    $newtype = "pdf";
                  } else {
                    unlink($tmpfile);
                    print STDERR "Failed to handle octet stream\n";
                    next;
                  }
                } elsif ($newtype eq "text/plain") {
                  $newtype = "other";
                } else {
                  print STDERR "Warning: unrecognized content type $newtype\n";
                  $newtype = "other";
                }
                $file = "$prefix.$newtype";
                rename($tmpfile,$file) || die "Error renaming $tmpfile to $file";

                if ($type eq "pdf") {
                  my $filereport = `file $file`;
                  if (! $filereport =~ m/PDF document/) {
                    print STDERR "Error: fetched $file for $doi but it is not in PDF format\n";
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
            print ACC join("\t", $pmcId, $pmId, "full")."\n";
            $last_pmId = $pmId;
            $last_text = $text;
        }
        die unless defined $text;
        my @snippets = TextSnippets($text, $queryTerm, $snippetBefore, $snippetAfter, $maxChar);
        foreach my $snippet (@snippets) {
            print OUT join("\t", $pmcId, $pmId, $queryTerm, $queryId, $snippet)."\n";
        }
    }
    close(OUT) || die "Error writing to $outfile";
    close(ACC) || die "Error writing to $accfile";
}
