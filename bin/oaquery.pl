#!/usr/bin/perl -w
use strict;
use lib "$ENV{HOME}/Genomics/Perl/modules";
use GenomicsUtils;
use Gene;

die "Run as a filter\n" unless @ARGV==0;

my %seen = ();

GenomicsUtils::connect('-host' => 'pub.microbesonline.org',
                       '-user' => 'guest',
                       '-pass' => 'guest',
                       '-dbname' => 'genomics')
  || die "Cannot connect to pub.microbesonline.org";

while(<STDIN>) {
    chomp;
    my ($pmcId, $locustag, $moIds) = split /\t/, $_;
    die "Wrong number of rows" unless defined $moIds && $moIds =~ m/\d/;
    my @moIds = split /,/, $moIds;
    foreach my $moId (@moIds) {
        next if exists $seen{$moId};
        $seen{$moId} = 1;
        my $gene = Gene::new('locusId' => $moId);
        if (!defined $gene) {
            print STDERR "No gene found in MicrobesOnline for $moId\n";
            next;
        }
        next unless $gene->type() eq "1";
        my $aaseq = $gene->protein();
        my $organism = $gene->taxonomyName();
        print join("\t", $organism, $locustag, "VIMSS".$moId, $aaseq, $gene->description())."\n";
    }
}
