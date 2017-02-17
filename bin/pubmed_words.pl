#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<END
Usage: pubmed_words.pl -list list_of_files > words
    The input file should have a list of files, each
    of which is the output of pubmedparse.pl
END
    ;

my $listfile;
die $usage
  unless GetOptions('list=s' => \$listfile)
  && defined $listfile
  && @ARGV == 0;

my %seen = ();
open(LIST, "<", $listfile) || die "Cannot read $listfile";
while (my $file = <LIST>) {
  chomp $file;
  open(IN, "<", $file) || die "Cannot read $file";
  while (my $line = <IN>) {
    chomp $line;
    my ($pmId, $abstract) = split /\t/, $line;
    my @words = split /\s+/, $abstract;
    foreach my $word (@words) {
        $word =~ s/[.,;-]+$//;
        if ($word =~ m/^[a-zA-Z][a-zA-Z90-9_]+\d\d\d$/) {
          print "$pmId\t$word\n" unless exists $seen{$word};
          $seen{$word} = 1;
        }
      }
  }
  close(IN) || die "Error reading $file";
}
close(LIST) || die "Error reading in files from $listfile";
