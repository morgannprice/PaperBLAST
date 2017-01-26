#!/usr/bin/perl -w
use strict;

my $usage = "Usage: combineSnippets.pl snippetsFile ... snippetsFile > uniqueSnippets\n";
die $usage if @ARGV == 0;

my @files = @ARGV;
foreach my $file (@files) {
    die "No such file: $file\n\n$usage" unless -e $file;
}

# pmcId::pmId::queryId => ifile for the first file that has snippets for this
my %known = ();
foreach my $ifile (0..(scalar(@files)-1)) {
    my $file = $files[$ifile];
    open(FILE, "<", $file) || die "Error reading $file";
    my $nSkip = 0;
    my $nPrint = 0;
    while(my $line = <FILE>) {
        chomp $line;
        my ($pmcId, $pmId, $queryTerm, $queryId, $snippet) = split /\t/, $line;
        die "Invalid input in $file" unless defined $snippet;
        my $key = join("::", $pmcId, $pmId, $queryId);
        if (exists $known{$key} && $known{$key} < $ifile) {
            $nSkip++;
        } else {
            $known{$key} = $ifile;
            print $line."\n";
            $nPrint++;
        }
    }
    close(FILE) || die "Error reading $file";
    print STDERR join("\t","Processed",$file,"Print",$nPrint,"Skip",$nSkip)."\n";
}
print STDERR "Processed " . scalar(@files) . " snippets files\n";
