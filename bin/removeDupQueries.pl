#!/usr/bin/perl -w
# This script is obsolete -- see uniqueQueries.pl instead
use strict;

my $usage = <<END
Usage: removeDupQueries queryprot_files < input_queries > cleaned_queries
    Removed queries from STDIN that match to queries in the queryprot_files.
    Issues warnings if queries have different sequences.
    If the description field is present, it maintains it.
END
;

die $usage if @ARGV == 0;

my %known = (); # organism_:::_query => sequence
my $sep = "_:::_";
foreach my $file (@ARGV) {
    open(FILE, "<", $file) || die "Cannot read $file\n";
    while (my $line = <FILE>) {
        chomp $line;
        my ($organism, $query, undef, $seq) = split /\t/, $line;
        die "Invalid line in $file: $line\n" unless defined $seq && $seq ne "";
        $known{$organism . $sep . $query} = $seq;
    }
    close(FILE) || die "Error reading $file";
}
print STDERR "Read " . scalar(keys %known) . " known queries to filter on\n";

my $nLines = 0;
my $nKept = 0;
my $nStartDiff = 0;
my $nIncon = 0;

while(my $line = <STDIN>) {
    chomp $line;
    my ($organism, $query, undef, $seq) = split /\t/, $line;
    die "Invalid line in STDIN: $line\n" unless defined $seq && $seq ne "";
    $nLines++;
    my $key  = $organism . $sep . $query;
    if (exists $known{$key}) {
        if ($seq ne $known{$key}) {
            # is this a real inconsistency, or just a difference in start codons?
            my $seq2 = $known{$key};
            my $len1 = length($seq);
            my $len2 = length($seq2);
            my $minlen = $len1 < $len2 ? $len1 : $len2;
            # after trimming to matching lengths, kip the first character because
            # of non-ATG start codons
            my $sub1 = substr($seq, length($seq)-$minlen+1);
            my $sub2 = substr($seq2, length($seq2)-$minlen+1);
            if ($sub1 eq $sub2) {
                $nStartDiff++;
            } else {
                $nIncon++;
                print STDERR "Warning: inconsistent sequences for $query in $organism\n";
            }
        }
    } else {
        print $line."\n";
        $nKept++;
    }
}
print STDERR "Kept $nKept of $nLines lines\n";
print STDERR "Among removed sequences, found $nStartDiff start codon disagreements and $nIncon other inconsistencies\n";
