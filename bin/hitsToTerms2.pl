#!/usr/bin/perl -w
# This script was used to help compute the list of popular genomes. It is obsolete.
use strict;
use JSON;
use lib "$ENV{HOME}/Genomics/Perl/modules";
use GenomicsUtils;
use Gene;

die "Run as a filter: hitsToTerms2.pl < id_and_json > terms\n"
    unless @ARGV==0;

GenomicsUtils::connect("localhost") || die "Cannot connect to MicrobesOnline";

my $nEmpty = 0;
my $nError = 0;
my $nResults = 0;
my $nWithHits = 0;
while(<STDIN>) {
    chomp;
    my ($moId, $jsonstring) = split /\t/, $_;
    die "Not enough fields in input" unless defined $jsonstring;
    if ($jsonstring eq "") {
        $nEmpty++;
        next;
    }
    my $json = from_json($jsonstring);
    if (!defined $json || !exists $json->{hitCount}) {
        $nError++;
        next;
    }
    $nResults++;
    if ($json->{hitCount} > 0) {
        $nWithHits++;
        my $query = $json->{request}{query};
        if (!defined $query) {
            $nError++;
            next;
        }
        my @parts = split / AND /, $query;
        die "Invalid query: $query" unless @parts == 2;
        my ($locustag,$genus) = @parts;
        my $gene = Gene::new('locusId' => $moId);
        if (!defined $gene) {
            print STDERR "No gene found in MicrobesOnline for $moId\n";
            next;
        }
        next unless $gene->type() eq "1";
        my $aaseq = $gene->protein();
        my $organism = $gene->taxonomyName();
        print join("\t", $organism, $locustag, "VIMSS".$moId, $aaseq)."\n";
    }
}
print STDERR "Warning: skipped $nEmpty empty lines\n" if $nEmpty > 0;
print STDERR "Warning: skipped $nError lines with JSON parse errors or search failures\n" if $nError > 0;
print STDERR "Parsed $nResults ($nWithHits with hits)\n";
