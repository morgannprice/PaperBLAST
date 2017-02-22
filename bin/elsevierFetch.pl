#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner;

my $usage = <<END
elsevierSnippets.pl -keyfile key -papers epmc_pop.papers -journals elsevier_journals2 -dir cache
	The keyfile should contain a single line with the APIkey.
	Places the fetched papers in cache/pubmed_*.xml
END
    ;

{
    my ($keyfile, $papersfile, $journalsfile, $cachedir);
    die $usage
        unless GetOptions('keyfile=s' => \$keyfile,
                          'papers=s' => \$papersfile,
                          'journals=s' => \$journalsfile,
                          'dir=s' => \$cachedir)
        && defined $keyfile && defined $papersfile && defined $journalsfile && defined $cachedir;

    die "No such file: $keyfile\n" unless -e $keyfile;
    die "No such file: $papersfile\n" unless -e $papersfile;
    die "No such file: $journalsfile\n" unless -e $journalsfile;
    die "No such directory: $cachedir\n" unless -d $cachedir;

    open(KEY, "<", $keyfile) || die "Cannot read $keyfile";
    my $key = <KEY>;
    chomp $key;
    die "Invalid key: $key" unless $key =~ m/^[a-zA-Z0-9]+$/;
    close(KEY) || die "Error reading $keyfile";

    my %journals = ();
    open(JOURNALS, "<", $journalsfile) || die "Cannot read $journalsfile";
    while (my $line = <JOURNALS>) {
        chomp $line;
        my @F = split /\t/, $line;
        $journals{$F[0]} = 1 if @F > 0;
    }
    close(JOURNALS) || die "Error reading $journalsfile";
    print STDERR "Read " . scalar(keys %journals) . " journals from $journalsfile\n";

    my $nTry = 0;
    my $nFetch = 0;
    my $nCache = 0;
    my $nFail = 0;
    my %fail = (); # pmId => 1 if failed to fetch it
    open(PAPERS, "<", $papersfile) || die "Cannot read $papersfile";
    while(my $line = <PAPERS>) {
      chomp $line;
      my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
      die "Invalid input in $papersfile" unless defined $isOpen;
      die "Invalid input in $papersfile: isOpen = $isOpen" unless $isOpen eq "1" || $isOpen eq "0";
      next unless $isOpen eq "0" && $pmId ne "" && exists $journals{$journal};
      my $file = "$cachedir/pubmed_$pmId.xml";
      my $file2 = "$cachedir/pubmed_$pmId.pdf";
      $nTry++;
      if (exists $fail{$pmId}) {
        $nFail++;
      } elsif (-e $file || -e $file2) {
        $nCache++;
      } elsif (! &ElsevierFetch($pmId, $key, $file, $journal)) {
        $fail{$pmId} = 1;
        $nFail++;
      } else {
        my $text = &XMLFileToText($file);
        if (!defined $text) {
          print STDERR "Warning: error reading XML for pubmed id $pmId\n";
          unlink($file);
          $nFail++;
        } else {
          $nFetch++;
        }
      }
    }
    print STDERR "Tried $nTry fetched $nFetch cached $nCache failed $nFail\n";
}
