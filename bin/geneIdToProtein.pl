#!/usr/bin/perl -w
use strict;
use LWP::Simple;
use lib ".";
use Bio::ASN1::EntrezGene;

my $usage = <<END
Run as a filter with the first column being the gene id.
END
;
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $batch_size = 50;
sub HandleBatch;

my %known = ();
my @batch = ();
while(my $line = <STDIN>) {
  chomp $line;
  my ($geneId) = split /\t/, $line;
  die "No geneId in $line" unless defined $geneId;
  next if $geneId =~ m/^gene/i; # skip any header line
  die "Invalid gene $geneId" unless $geneId =~ m/^\d+$/;

  next if exists $known{$geneId};
  $known{$geneId} = 1;
  push @batch, $geneId;
  if (@batch >= $batch_size) {
    &HandleBatch(@batch);
    @batch = ();
  }
}
&HandleBatch(@batch);

sub HandleBatch {
  my @geneIds = @_;
  return if @geneIds == 0;
  my $db = 'Gene';
  my $url = $base . "efetch.cgi?db=$db&id=" . join(",",@geneIds);
  my $results = get($url);
  my @parts = split /Entrezgene ::= /, $results;
  shift @parts; # the first one should be empty
  print STDERR "Warning: " . scalar(@geneIds) . " queries, $geneIds[0] to $geneIds[-1], but got " . scalar(@parts) . " results\n"
    unless scalar(@parts) == scalar(@geneIds);
  my $nErr = 0;
  foreach my $asn1 (@parts) {
    my $parser = Bio::ASN1::EntrezGene->new();
    my $parsed;
    eval { $parsed = $parser->parse($asn1); };
    if ($@) {
      print STDERR "Parsing failed for " . join(",", @geneIds)."\n";
      $nErr++;
      my $saveFile = "/tmp/$$.$geneIds[0].$nErr";
      print STDERR "Saving the unparseable item to $saveFile\n";
      open(ERR, ">", $saveFile) || die "Cannot write to $saveFile";
      print ERR $asn1;
      close(ERR) || die "Error writing to $saveFile";
      next;
    }
    my $gene = $parsed->[0];
    die unless keys(%$gene) > 0;
    my $info = $gene->{"track-info"};
    die "Parsing failed: no track-info" unless defined $info;
    my $geneId = $info->[0]{geneid};
    die "Parsing failed: no gene id" unless defined $geneId && $geneId =~ m/^\d+$/;
    if ($gene->{type} eq "protein-coding") {
      if (!exists $gene->{locus}) {
        print STDERR "Warning: $geneId is protein-coding but has no locus entries\n";
      } elsif (ref($gene->{locus}) ne 'ARRAY') {
        print STDERR "Warning: $geneId is protein-coding but locus is not an array\n";
      } else {
        my @loci = @{ $gene->{locus} };
        @loci = grep { exists $_->{products} } @loci;
        if (@loci == 0) {
          print STDERR "Warning: $geneId is protein-coding but has no locus/products entries\n";
        } else {
          my $products = $loci[0]{products};
          if (@$products == 0) {
            print STDERR "Warning: $geneId is protein-coding with empty productn";
          } else {
            my $product = $products->[0];
            if (! $product->{accession}) {
              print STDERR "Warning: $geneId is protein-coding with no accession\n";
            } else {
              print join("\t", $geneId, $product->{accession})."\n";
            }
          }
        }
      }
    } else {
      print STDERR "Not a protein: $geneId\n";
    }
  }
}

