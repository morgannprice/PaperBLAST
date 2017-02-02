#!/usr/bin/perl -w
use strict;
# If snippets were identified in an earlier file, ignore snippets in later files
# (I.e., prefer snippets from full text to snippets from just the abstract)
# Also ignores repeating lines, which can happen in rare caess
# (i.e. redundant abstract entries in pubmed)

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
    my $lastline = undef;
    while(my $line = <FILE>) {
        chomp $line;
        my ($pmcId, $pmId, $queryTerm, $queryId, $snippet) = split /\t/, $line;
        die "Invalid input in $file" unless defined $snippet;
        my $key = join("::", $pmcId, $pmId, $queryId);
        if (exists $known{$key} && $known{$key} < $ifile) {
            $nSkip++;
        } else {
            $known{$key} = $ifile;
            print $line."\n" unless defined $lastline && $line eq $lastline;
            $lastline = $line;
            $nPrint++;
        }
    }
    close(FILE) || die "Error reading $file";
    print STDERR join("\t","Processed",$file,"Print",$nPrint,"Skip",$nSkip)."\n";
}
print STDERR "Processed " . scalar(@files) . " snippets files\n";
