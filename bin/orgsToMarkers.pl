#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{ReadTable ReadFastaEntry};

my $markersFile = "$Bin/../gaps/markers.tab";
my $usage = <<END
orgsToMarkers.pl -table in -faa orgsDef.faa -out ref.fasta

The input table must include gdb, gid, and gtax (which is used to
determine if a genome is bacteria or archaea).

The fasta input must have headers of the form gdb__gid:locusid (with
optional whitespace after that which is ignored)

Optional arugments:
-markers $markersFile -- must include Domain (Bacteria or Archaea), Name, and acc
-verbose -- additional logging output
END
;

my ($table, $faaIn, $faaOut, $verbose);
die $usage
  unless GetOptions('table=s' =>\$table,
                    'faa=s' => \$faaIn,
                    'out=s' => \$faaOut,
                   'markers=s' => \$markersFile,
                   'verbose' => \$verbose)
  && @ARGV == 0
  && defined $table && defined $faaIn && defined $faaOut;

my @orgs = ReadTable($table, qw{gdb gid gtax});
my %orgs = map { $_->{gdb} . "__" . $_->{gid} => $_ } @orgs;
print STDERR "Working on " . scalar(keys %orgs) . " organisms\n";

my @markers = ReadTable($markersFile, qw{Domain Name acc});
foreach my $m (@markers) {
  die "Unknown domain $m->{Domain} in $markersFile"
    unless $m->{Domain} eq "Bacteria" || $m->{Domain} eq "Archaea";
  die "Invalid name $m->{Name}" if $m->{Name} eq "" || $m->{Name} =~ m/\s/;
}

my $hmmsearch = "hmmer/hmmsearch";
die "No such executable: $hmmsearch" unless -x $hmmsearch;

my %seq = (); # organism to locusId to sequence
open(my $fhIn, "<", $faaIn) || die "Cannot read $faaIn\n";
my $state = {};
while (my ($header, $seq) = ReadFastaEntry($fhIn, $state)) {
  $header =~ s/ .*//;
  $header =~ m/^([^:_]+__[^:]+):(.*)$/ || die "Cannot parse gdb__gid:locusId from $header";
  my ($orgId,$locusId) = ($1,$2);
  next unless exists $orgs{$orgId};
  die "Duplicate sequence for $orgId:$locusId" if exists $seq{$orgId}{$locusId};
  $seq{$orgId}{$locusId} = $seq;
}
close($fhIn) || die "Error reading $faaIn\n";
print STDERR "Found sequences for " . scalar(keys %seq) . " organisms\n";

my %marker = (); # orgId => marker => list of locusIds

while (my ($orgId, $seqs) = each %seq) {
  my $domain;
  my $gtax = $orgs{$orgId}{gtax};
  if ($gtax =~ m/^d__Bacteria/) {
    $domain = "Bacteria";
  } elsif ($gtax =~ m/^d__Archaea/) {
    $domain = "Archaea";
  } else {
    die "Cannot determine Bacteria or Archaea from taxonomy string: $gtax";
  }
  print STDERR "Analyzing $orgId from domain $domain with " . scalar(keys %$seqs) . " proteins\n";

  # Build the temporary input file
  my $faaTmp = "/tmp/orgsToMarkers.$$.tmp";
  open(my $fhTmp, ">", $faaTmp) || die "Cannot write to $faaTmp\n";
  while (my ($locusId,$seq) = each %$seqs) {
    print $fhTmp ">$locusId\n$seq\n";
  }
  close($fhTmp) || die "Error writing to $faaTmp";

  foreach my $m (@markers) {
    next unless $m->{Domain} eq $domain;
    my $hmmFile = "path.aa/$m->{acc}.hmm";
    die "No hmm for $m->{acc} -- $hmmFile does not exist" unless -e $hmmFile;
    my @cmd = ($hmmsearch,
               "--cpu", $ENV{MC_CORES} || 1,
               "--cut_tc",
               "-o", "/dev/null",
               "--tblout", "$faaTmp.hits",
               $hmmFile, $faaTmp);
    print STDERR "Running: @cmd\n" if defined $verbose;
    system(@cmd) == 0 || die "hmmsearch failed: $!";
    open(my $fhHits, "<", "$faaTmp.hits") || die "Cannot read $faaTmp.hits";
    my %loci = ();
    while (my $line = <$fhHits>) {
      chomp $line;
      next if $line =~ m/^#[ -]/;
      next unless $line =~ m/\s(TIGR|PF)/;
      $line =~ s/ .*//;
      $loci{$line} = 1;
    }
    my @loci = keys %loci;
    print STDERR "$orgId $m->{Name} @loci\n" if defined $verbose;
    $marker{$orgId}{$m->{Name}} = \@loci;
  }
  unlink($faaTmp);
  unlink("$faaTmp.hits");
  print STDERR "Warning, no marker genes found for $orgId\n"
    unless exists $marker{$orgId};
}

# This is necessary because some marker names are in both the bacterial and archaeal sections
my %written = (); # orgId => marker => 1
open(my $fhOut, ">", $faaOut) || die "Cannot write to $faaOut";
foreach my $orgId (sort keys %marker) {
  foreach my $m (@markers) {
    my $mname = $m->{Name};
    if (exists $marker{$orgId}{$mname} && !exists $written{$orgId}{$mname}) {
      my $list = $marker{$orgId}{$mname};
      if (@$list == 1) {
        my $locusId = $list->[0];
        my $seq = $seq{$orgId}{$locusId};;
        die unless defined $seq;
        print $fhOut ">$orgId:$mname\n$seq\n";
        $written{$orgId}{$mname} = 1;
      }
    }
  }
}
close($fhOut) || die "Error writing to $faaOut\n";
print STDERR "Wrote $faaOut\n";
