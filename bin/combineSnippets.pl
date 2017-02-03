#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<END
Usage: combineSnippets.pl -out uniqueSnippets snippetsFile ... snippetsFile
    Expects the snippets to be tab delimited with
    pmcId, pmId, queryTerm, queryId, snippet

    If snippets were identified in an earlier file, ignore snippets in later files
    (I.e., prefer snippets from full text to snippets from just the abstract, assuming
	the input files are ordered correctly.)
    Also ignores repeating lines, which can happen in rare caess
    (i.e. redundant abstract entries in pubmed)

    For each snippet file, there should also be a snippet.access file, with fields
    pmcId, pmId, and either "full" or "abstract". This reports the access used, and warns
    if the access used is not the best access.
END
    ;

my $out;
die $usage unless
    GetOptions('out=s' => \$out)
    && defined $out;
die "No input files\n$usage" if @ARGV == 0;

my @files = @ARGV;
foreach my $file (@files) {
    die "No such file: $file\n\n$usage" unless -e $file;
    die "No access file: $file.access\n\n$usage" unless -e "$file.access";
}

# pmcId::pmId::queryId => ifile for the first file that has snippets for this
open(OUT, ">", $out) || die "Cannot write to $out";
open(ACCOUT, ">", "$out.access") || die "Cannot write to $out.access";
my %known = (); # pmcId::pmId::queryTerm => 1 if seen in a previous file
my %knownPair = (); # pmcId::pmId => 1 if written out already
foreach my $ifile (0..(scalar(@files)-1)) {
    my %keyToAcc = (); # pmcId::pmId => access
    my $file = $files[$ifile];
    my $accfile = "$file.access";
    open(ACC, "<", $accfile) || die "Error reading $accfile";
    while(my $line = <ACC>) {
        chomp $line;
        my ($pmcId, $pmId, $acc) = split /\t/, $line;
        die "Invalid line in $accfile: $line" unless defined $acc;
        die "No pmcId or pmId in $line in $accfile" if $pmcId eq "" && $pmId eq "";
        my $key = join("::", $pmcId, $pmId);
        if (exists $keyToAcc{$key}) {
            die "Inconsistent access for pmc $pmcId pm $pmId in $accfile\n"
                unless $keyToAcc{$key} eq $acc;
        } else {
            $keyToAcc{$key} = $acc;
            print ACCOUT join("\t", $pmcId, $pmId, $keyToAcc{$key})."\n"
                unless exists $knownPair{$key};
            $knownPair{$key} = 1;
        }
    }

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
            print OUT $line."\n" unless defined $lastline && $line eq $lastline;
            $lastline = $line;
            $nPrint++;
        }
    }
    close(FILE) || die "Error reading $file";
    print STDERR join("\t","Processed",$file,"Print",$nPrint,"Skip",$nSkip)."\n";
}
close(OUT) || die "Error writing to $out";
close(ACCOUT) || die "Error writing to $out.access";
print STDERR "Processed " . scalar(@files) . " snippets files\n";
