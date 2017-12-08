#!/usr/bin/perl -w
# This script is obsolete -- use parse_ecocyc.pl instead
use strict;
# Output format is tab-delimited with:
# Ecocyc protein id (i.e., G6760-MONOMER)
# b number (i.e., b1449)
# short protein name (i.e., YncB)
# protein description (may contain HTML entities or <sup> tags)
# comma-delimited list of PubMed ids (or empty)

die "Run as a filter from the EcoCyc proteins.dat file\n"  unless @ARGV == 0;

sub ProcessRecord;

my @lines = ();
while(<STDIN>) {
    next if m/^#/;
    chomp;
    push @lines, $_;
    if ($_ eq "//") {
        ProcessRecord(@lines);
        @lines = ();
    }
}
ProcessRecord(@lines);

sub ProcessRecord {
    my @lines = @_;
    return if @lines == 0;
    my @polypeptides = grep m/^TYPES - Polypeptides/, @lines;
    return unless @polypeptides > 0;

    my @uniqlines = grep m/^UNIQUE-ID/, @lines;
    die "No UNIQUE-ID entry" unless @uniqlines > 0;
    die "More than one unique-id $uniqlines[0] $uniqlines[1]" unless @uniqlines == 1;
    $uniqlines[0] =~ m/^UNIQUE-ID - (.*)$/ || die "Cannot parse $uniqlines[0]";
    my $uniqId = $1;

    my @blinks = grep m/^DBLINKS/ && m/ECOLIWIKI/, @lines;
    my $bnum = "";
    if (@blinks == 1) {
        $blinks[0] =~ m/"([bB]\d+)"/ || die "Cannot find b number in $blinks[0]";
        $bnum = $1;
    }
    
    my @syn = grep m/^ABBREV-NAME|^SYNONYMS/, @lines;
    my $name = "";
    if (@syn > 0 && $syn[0] =~ m/^[A-Z-]+ - ([a-zA-Z0-9_.-]+)$/) {
        $name = $1;
    }

    # extract "COMMON NAME" or description
    my $common_name = "";
    my @commlines = grep m/^COMMON-NAME/, @lines;
    if (@commlines > 0 && $commlines[0] =~ m/^COMMON-NAME - (.*)$/) {
        $common_name = $1;
    }

    # and extract pubmed ids
    my @pubmedIds = ();
    foreach my $line (@lines) {
        while ($line =~ m/[|]CITS: \[(\d+)\][|]/g) {
            push @pubmedIds, $1;
        }
        if ($line =~ m/^\^CITATIONS - (\d+):/) {
            push @pubmedIds, $1;
        }
    }
    my %pubmedIds = map { $_ => 1 } @pubmedIds;
    print join("\t", $uniqId, $bnum, $name, $common_name,
               join(",", sort { $a <=> $b } keys %pubmedIds))."\n";

}
