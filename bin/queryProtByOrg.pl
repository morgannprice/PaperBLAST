#!/usr/bin/perl -w
use strict;
use lib "$ENV{HOME}/Genomics/Perl/modules";
use GenomicsUtils;
use Genome;
use Gene;

die "Run as a filter on a tab-delimited file with organism names\n"
  unless @ARGV == 0;

GenomicsUtils::connect('-host' => 'pub.microbesonline.org',
                       '-user' => 'guest',
                       '-pass' => 'guest',
                       '-dbname' => 'genomics')
  || die "Cannot connect to pub.microbesonline.org";

my %nameToTaxId = ();
my $q = GenomicsUtils::query("SELECT name, taxonomyId FROM Taxonomy");
foreach my $row (@$q) {
  my ($taxName,$taxId) = @$row;
  $nameToTaxId{$taxName} = $taxId;
}
print STDERR "Loaded " . scalar(keys %nameToTaxId) . " taxonomy names\n";

my $nShow = 0;
while (<STDIN>) {
  chomp;
  my ($taxname) = split /\t/, $_;
  die "Invalid input" unless defined $taxname && $taxname ne "";
  my $taxId = $nameToTaxId{$taxname};
  if (!defined $taxId) {
    print STDERR "Skipping unknown taxon $taxname\n";
    next;
  }
  my $genome = Genome::new('taxonomyId' => $taxId);
  my @genes = $genome->genes;
  push @genes, $genome->nongenomicGenes;
  foreach my $gene (@genes) {
    next unless $gene->type eq "1";
    my @locustags = ();
    push @locustags, $gene->synonym(1);
    push @locustags, $gene->synonym(4);
    push @locustags, $gene->synonym(6);
    my %locustags = map { $_ => 1 } @locustags;
    @locustags = keys %locustags;
    foreach my $locustag (@locustags) {
      my $prot = $gene->protein;
      print join("\t", $taxname, $locustag, "VIMSS".$gene->locusId(), $gene->protein, $gene->description())."\n"
        if $prot;
    }
    $nShow++;
  }
}
print STDERR "Printed $nShow genes\n";
