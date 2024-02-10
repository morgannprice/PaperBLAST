#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{getMicrobesOnlineDbh};

die "Run as a filter on a tab-delimited file with organism names\n"
  unless @ARGV == 0;

my $dbh = getMicrobesOnlineDbh() || die "No local copy of MicrobesOnline\n";

my %nameToTaxId = ();
my $q = $dbh->selectall_arrayref("SELECT name, taxonomyId FROM Taxonomy");
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
  my $genes = $dbh->selectall_arrayref("SELECT * from Scaffold JOIN Locus USING (scaffoldId)
                                        WHERE taxonomyId = ? AND isActive=1 AND priority=1 AND type=1",
                                       { Slice => {} }, $taxId);
  foreach my $gene (@$genes) {
    my $syns = $dbh->selectall_arrayref("SELECT type,name from Synonym WHERE locusId = ? AND version = ?",
                                       {}, $gene->{locusId}, $gene->{version});
    my @locustags = ();
    foreach my $row (@$syns) {
      my ($type, $name) = @$row;
      push @locustags, $name if $type == 1 || $type == 4 || $type == 6;
    }
    my %locustags = map { $_ => 1 } @locustags;
    @locustags = keys %locustags;
    my ($aaseq) = $dbh->selectrow_array("SELECT sequence FROM AASeq WHERE locusId = ? AND version = ?",
                                           {}, $gene->{locusId}, $gene->{version});
    next unless $aaseq;
    my ($desc) = $dbh->selectrow_array("SELECT description FROM Description WHERE locusId = ? AND version = ?",
                                       {}, $gene->{locusId}, $gene->{version});
    foreach my $locustag (@locustags) {
      print join("\t", $taxname, $locustag, "VIMSS".$gene->{locusId}, $aaseq, $desc || "")."\n";
    }
    $nShow++;
  }
}
print STDERR "Printed $nShow genes\n";
