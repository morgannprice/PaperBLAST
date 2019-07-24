#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{ReadTable};

my $nCPU = $ENV{MC_CORES} || 4;
my $minNMarkers = 10;
my $minMedian = 75.0;
my $minIdentity = 50;
my $usage = <<END
Usage: orgsVsMarkers.pl -orgs orgs -vs markers.faa -out out.tsv

The output table includes orgId (from orgs), orgId2 (from the marker
file), identity (up to 100), and nMarkers (the number of protein
alignments used to compute the median identity)

Optional arguments:
 -nCPU $nCPU -- number of CPUs to use (defaults to the MC_CORES
	environment variable, or 4)
 -minNMarkers $minNMarkers -- number of marker proteins that must have good
	alignments in order to estimate percent identity
 -minMedian $minMedian -- lowest median identity to report
 -minIdentity $minIdentity -- minimum pairwise identity of marker
	proteins to consider

Limitations -- the use of usearch with maxaccepts may cause all hits
to be missed if the reference database includes many strains that are
at about the same distance from the query genome.
END
;

my ($orgsPre, $markerFaa, $outFile);
die $usage unless
  GetOptions('orgs=s' => \$orgsPre,
             'vs=s' => \$markerFaa,
             'out=s' => \$outFile,
             'nCPU=i' => \$nCPU,
             'minNMarkers=i' => \$minNMarkers,
             'minMedian=f' => \$minMedian,
             'minIdentity=f' => \$minIdentity)
  && @ARGV == 0
  && defined $orgsPre && defined $markerFaa && defined $outFile;

my $orgsFaa = "$orgsPre.faa";
foreach my $file ($orgsFaa, $markerFaa) {
  die "No such file: $file\n" unless -e $file;
}
my @orgs = ReadTable("$orgsPre.org", qw{orgId gdb gid genomeName});
die "No organisms in $orgsPre.org\n" unless @orgs > 0;
my %orgs = map { $_->{orgId} => $_ } @orgs;

my $usearch = "$Bin/usearch";
die "No such executable: $usearch\n" unless -x $usearch;

my $tmpFile = "$outFile.$$.tmp";
my @cmd = ($usearch, "-usearch_global", $orgsFaa, "-db", $markerFaa,
           "-id", $minIdentity/100, "-maxaccepts", 20, "-maxrejects", 20,
           "-query_cov", 0.7, "-target_cov", 0.7, "-blast6out", $tmpFile,
           "-threads", $nCPU, "--quiet");
system(@cmd) == 0 || die "usearch failed: $!";

my %locusHits = (); # orgId => locusId => orgId2 => list of [marker,identity]
my %orgHits = (); # orgId => orgId2 => marker => list of [locusId,identity]

open(my $fhIn, "<", $tmpFile) || die "Cannot read $tmpFile";
while (my $line = <$fhIn>) {
  chomp $fhIn;
  my ($query, $subject, $identity) = split /\t/, $line;
  $query =~ s/ .*//; # remove descrption (if any)
  my ($orgId, $locusId) = split /:/, $query;
  die "Invalid query $query" unless defined $locusId && $locusId ne "";
  die "Invalid query $query -- unknown orgId $orgId" unless exists $orgs{$orgId};
  my ($orgId2, $marker) = split /:/, $subject;
  die "Invalid subject $subject" unless $marker;
  push @{ $locusHits{$orgId}{$locusId}{$orgId2} }, [$marker, $identity];
  push @{ $orgHits{$orgId}{$orgId2}{$marker} }, [$locusId, $identity];
}
close($fhIn) || die "Error reading $tmpFile";
unlink($tmpFile);

open(my $fhOut, ">", $tmpFile) || die "Cannot write to $tmpFile";
print$fhOut join("\t", "orgId", "orgId2", "identity", "nMarkers")."\n";
foreach my $org (@orgs) {
  my $orgId = $org->{orgId};
  my @hits = (); # each row has orgId2, identity (median%), and n (#proteins that match)
  while (my ($orgId2, $hash) = each %{ $orgHits{$orgId} }) {
    next unless scalar(keys %$hash) >= $minNMarkers;
    # compute actual list of unique markers, and then the median %identity
    my @identity = ();
    while (my ($marker, $list) = each %$hash) {
      next if @$list > 1; # >1 alignment for this marker
      my ($locusId, $identity) = @{ $list->[0] };
      push @identity, $identity
        if scalar(@{ $locusHits{$orgId}{$locusId}{$orgId2} }) == 1; # locus matches exactly 1 marker
    }
    my $n = scalar(@identity);
    next unless $n >= $minNMarkers && $n >= 1;
    @identity = sort { $a <=> $b } @identity; # numeric sort
    # i.e. n=2 => avg($id[0], $id[1]); n=3 => $id[1]
    my $med = $n % 2 == 0 ? ($identity[$n/2 - 1] + $identity[$n/2])/2 : $identity[int($n/2)];
    push @hits, [$orgId2, $med, $n]
      if $med >= $minMedian;
  }
  @hits = sort { $b->[1] <=> $a->[1] } @hits; # highest %identity first
  foreach my $hit (@hits) {
    my ($orgId2, $median, $n) = @$hit;
    print $fhOut join("\t", $orgId, $orgId2, $median, $n)."\n";
  }
}
close($fhOut) || die "Error writing to $tmpFile\n";
rename($tmpFile, $outFile) || die "Cannot rewrite to $outFile\n";
