#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{getMicrobesOnlineDbh};

die "Run as a filter\n" unless @ARGV==0;

my $dbh = getMicrobesOnlineDbh() || die "No local copy of MicrobesOnline\n";

my %seen = ();

while(<STDIN>) {
    chomp;
    my ($pmcId, $locustag, $moIds) = split /\t/, $_;
    die "Wrong number of rows" unless defined $moIds && $moIds =~ m/\d/;
    my @moIds = split /,/, $moIds;
    foreach my $moId (@moIds) {
        next if exists $seen{$moId};
        $seen{$moId} = 1;
        my $gene = $dbh->selectrow_hashref("SELECT * from Locus JOIN Scaffold USING (scaffoldId)
                                             JOIN Taxonomy USING (taxonomyId)
                                             WHERE locusId = ? AND priority = 1 AND isActive = 1",
                                            { Slice => {} }, $moId);
        if (!defined $gene) {
          print STDERR "No gene found in MicrobesOnline for $moId\n";
          next;
        }
        next unless $gene->{type} eq "1";
        my ($aaseq) = $dbh->selectrow_array("SELECT sequence from AASeq
                                             WHERE locusId = ? AND version = ?",
                                            {}, $moId, $gene->{version});
        next unless $aaseq;
        my $organism = $gene->{name};
        my ($description) = $dbh->selectrow_array("SELECT description FROM Description
                                                   WHERE locusId = ? AND version = ?",
                                                  {}, $moId, $gene->{version});
        print join("\t", $organism, $locustag, "VIMSS".$moId, $aaseq,
                   $description || "")."\n";
    }
}
