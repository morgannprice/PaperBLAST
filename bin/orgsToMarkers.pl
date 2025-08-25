#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use File::Basename;
use lib "$RealBin/../lib";
use pbutils qw{ReadTable ReadFastaEntry NewerThan};

my $markersFile = "$RealBin/../gaps/markers.tab";
my $minMarkerFrac = 0.7;
my $usage = <<END
orgsToMarkers.pl -org orgs -out ref.fasta

The input table would be orgs.org and must include gdb and gid. The
fasta input would be orgs.faa and must have headers of the form
gdb__gid:locusid (with optional whitespace after that which is
ignored). Genomes are automatically classified as bacteria or archaea
based on which set of markers have a higher proportion of hits.

Note that hmm analyses are cached in the parent directory of the orgs,
under the file names hmm.hits

Optional arguments:
-markers $markersFile -- must include Domain (Bacteria or Archaea), Name, and acc
-verbose -- additional logging output
-minMarkerFrac $minMarkerFrac -- warn if a genome has hits to less than this
  fraction of markers for both bacteria and archaea.
END
;

my ($orgsIn, $faaOut, $verbose);
die $usage
  unless GetOptions('org=s' =>\$orgsIn,
                    'out=s' => \$faaOut,
                    'markers=s' => \$markersFile,
                    'verbose' => \$verbose,
                    'minMarkerFrac=f' => \$minMarkerFrac)
  && @ARGV == 0
  && defined $orgsIn && defined $faaOut;
die "Invalid minMarkerFrac $minMarkerFrac\n"
  unless $minMarkerFrac > 0 && $minMarkerFrac <= 1;
my @orgs = ReadTable("$orgsIn.org", qw{gdb gid});
my $cacheDir = dirname($orgsIn);

my %orgs = map { $_->{gdb} . "__" . $_->{gid} => $_ } @orgs;
print STDERR "Working on " . scalar(keys %orgs) . " organisms\n";

my @markers = ReadTable($markersFile, qw{Domain Name acc});
foreach my $m (@markers) {
  die "Unknown domain $m->{Domain} in $markersFile"
    unless $m->{Domain} eq "Bacteria" || $m->{Domain} eq "Archaea";
  die "Invalid name $m->{Name}" if $m->{Name} eq "" || $m->{Name} =~ m/\s/;
  die "Invalid acc $m->{acc}" if $m->{acc} eq "" || $m->{acc} =~ m/\s/;
}

my $hmmsearch = "$RealBin/hmmsearch";
die "No such executable: $hmmsearch" unless -x $hmmsearch;

my %seq = (); # organism to locusId to sequence
open(my $fhIn, "<", "$orgsIn.faa") || die "Cannot read $orgsIn.faa\n";
my $state = {};
while (my ($header, $seq) = ReadFastaEntry($fhIn, $state)) {
  $header =~ s/ .*//;
  $header =~ m/^([^:_]+__[^:]+):(.*)$/ || die "Cannot parse gdb__gid:locusId from $header";
  my ($orgId,$locusId) = ($1,$2);
  next unless exists $orgs{$orgId};
  die "Duplicate sequence for $orgId:$locusId" if exists $seq{$orgId}{$locusId};
  $seq{$orgId}{$locusId} = $seq;
}
close($fhIn) || die "Error reading $orgsIn.faa\n";
print STDERR "Found sequences for " . scalar(keys %seq) . " organisms\n";
foreach my $orgId (sort keys %orgs) {
  die "No sequences for $orgId in $orgsIn.faa\n"
    unless exists $seq{$orgId};
}

my @bacAcc = map $_->{acc}, grep $_->{Domain} eq "Bacteria", @markers;
die "No bacterial markers\n" unless @bacAcc > 0;
my @arcAcc = map $_->{acc}, grep $_->{Domain} eq "Archaea", @markers;
die "No archaeal markers\n" unless @arcAcc > 0;

my %acc = (); # accession to list of rows
# (some markers are shared between Bacteria and Archaea)
foreach my $m (@markers) {
  push @{ $acc{$m->{acc}} }, $m;
}

my $tmpDir = $ENV{TMPDIR} || "/tmp";
my $tmpPre = "$tmpDir/orgsToMarkers.$$.";

my %marker = (); # orgId => marker => list of locusIds
foreach my $acc (sort keys %acc) {
  my $hmmFile = "$RealBin/../tmp/path.aa/$acc.hmm";
  die "No hmm for $acc -- $hmmFile does not exist" unless -e $hmmFile;
  my $hitsFile = "$cacheDir/$acc.hits";
  if (-e $hitsFile && NewerThan($hitsFile, $hmmFile) && NewerThan($hitsFile, "$orgsIn.faa")) {
    print STDERR "Using cached hits $hitsFile\n";
  } else {
    my @cmd = ($hmmsearch,
               "--cpu", $ENV{MC_CORES} || 1,
               "--cut_tc",
               "-o", "/dev/null",
               "--tblout", $hitsFile,
               $hmmFile, "$orgsIn.faa");
    print STDERR "Running: @cmd\n" if defined $verbose;
    system(@cmd) == 0 || die "hmmsearch failed: $!";
  }
  open(my $fhHits, "<", $hitsFile) || die "Cannot read $hitsFile";
  my @loci = ();
  while (my $line = <$fhHits>) {
    chomp $line;
    next if $line =~ m/^#[ -]/;
    next unless $line =~ m/\s(TIGR|PF)/;
    $line =~ s/ .*//;
    push @loci, $line;
  }
  close($fhHits) || die "Error reading $hitsFile";
  foreach my $locus (@loci) {
    $locus =~ m/^([^_]+)__([^:]+):(.*)/ || die "Cannot parse locus $locus in $hitsFile";
    my ($gdb, $gid, $locusId) = ($1,$2,$3);
    my $orgId = $gdb . "__" . $gid;
    push @{ $marker{$orgId}{$acc} }, $locusId;
    print STDERR "$orgId $acc $locusId\n" if defined $verbose;
  }
}

foreach my $orgId (sort keys %orgs) {
  print STDERR "Warning, no marker genes found for $orgId\n"
    unless exists $marker{$orgId};
}

open(my $fhOut, ">", $faaOut) || die "Cannot write to $faaOut\n";
foreach my $orgId (sort keys %marker) {
  my $m = $marker{$orgId};
  my @bacHas = grep exists $m->{$_}, @bacAcc;
  my @arcHas = grep exists $m->{$_}, @arcAcc;
  my $fBac = scalar(@bacHas) / scalar(@bacAcc);
  my $fArc = scalar(@arcHas) / scalar(@arcAcc);
  if ($fBac < $minMarkerFrac && $fArc < $minMarkerFrac) {
    print STDERR "Warning: low fraction of markers found for $orgId -- bacterial $fBac archaeal $fArc\n";
  }
  my @accUsed;
  if ($fBac > $fArc) {
    @accUsed = @bacHas;
  } else {
    @accUsed = @arcHas;
  }
  foreach my $acc (@accUsed) {
    my $loci = $marker{$orgId}{$acc};
    die unless defined $loci && @$loci > 0;
    if (@$loci == 1) {
      my $locusId = $loci->[0];
      my $seq = $seq{$orgId}{$locusId};;
      die unless defined $seq;
      my $mrow = $acc{$acc}[0];
      my $mname = $mrow->{Name};
      die $acc unless defined $mname && $mname ne "";
      print $fhOut ">$orgId:$mname\n$seq\n";
    }
  }
}
close($fhOut) || die "Error writing to $faaOut\n";
print STDERR "Wrote $faaOut\n";
