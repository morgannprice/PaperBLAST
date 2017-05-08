#!/usr/bin/perl -w
use strict;
use lib "$ENV{HOME}/Genomics/Perl/modules";
use GenomicsUtils;

die "Usage: moIds.pl < words > words.filtered\n"
    unless @ARGV==0;

{
    GenomicsUtils::connect('-host' => 'pub.microbesonline.org',
                           '-user' => 'guest',
                           '-pass' => 'guest',
                           '-dbname' => 'genomics')
        || die "Cannot connect to pub.microbesonline.org";
    my $rows;
    # type 1 for standard locus tag
    # type 3 for NCBI accession
    # type 4 for alternate locus tag
    # type 6 for deprecated locus tag
    $rows = query("SELECT locusId,Synonym.type,name FROM Locus JOIN Synonym USING (locusId,version)
                          WHERE priority=1 AND Synonym.type IN (1,3,4,6);");
    my %nameToLocus = ();
    foreach my $row (@$rows) {
        my ($locusId,$type,$name) = @$row;
        # remove version number from NCBI accessions
        $name =~ s/[.]\d+// if $type == 3;
        push @{ $nameToLocus{$name} }, $locusId;
    }

    print STDERR "Loaded " . scalar(keys %nameToLocus) . " ids from MicrobesOnline\n";

    while(<STDIN>) {
        chomp;
        my ($pmcid,$term) = split /\t/, $_, -1;
        if (exists $nameToLocus{$term}) {
            # often have multiple entries for the same one so uniqify
            my %locusIds = map {$_ => 1} @{$nameToLocus{$term}};
            print join("\t", $pmcid, $term, join(",", sort keys %locusIds))."\n";
        }
    }
}
