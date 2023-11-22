#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};

my $identity = 0.7;
my $coverage = 0.8;

my $usage = <<END
Usage: clusterRegPrecise.pl -in regprecise.curated_parsed > regprecise_clustered.curated_parsed

Given the output from fetchRegPrecise.pl, removes ambiguous items (two VIMSS ids mapped
to the same TF in the same organism, or the same VIMSS id mapped to two TFs), and
then cluster (using usearch) and output a smaller table.

Optional arguments:
-identity $identity
-coverage $coverage
END
;

my $inFile;
die $usage
  unless GetOptions('in=s' => \$inFile,
                    'identity=f' => \$identity)
  && defined $inFile
  && @ARGV == 0;
die "Invalid identity $identity: must be between 0 and 1\n"
  unless $identity > 0 && $identity < 1;

my @rows = ();
open(my $fhIn, "<", $inFile) || die "Cannot read $inFile\n";
while (my $line = <$fhIn>) {
  chomp $line;
  # id is the regulon id; id2 is the VIMMS id; name is the TF name
  my ($db, $id, $id2, $name, $desc, $org, $seq, $comment, $pmIds) = split /\t/, $line;
  die "Not enough fields\n" unless defined $pmIds;
  die "Invalid sequence $seq\n" unless $seq =~ m/^[A-Z]+$/;
  push @rows, { 'db' => $db, 'id' => $id, 'id2' => $id2, 'name' => $name,
                'desc' => $desc, 'org' => $org, 'seq' => $seq, 'comment' => $comment, 'pmIds' => $pmIds };
}
close($fhIn) || die "Error reading $inFile";
print STDERR "Read " . scalar(@rows) . " entries\n";

my $usearch = "$RealBin/usearch";
die "No such executable: $usearch\n" unless -x $usearch;

# Remove non-unique VIMSS ids or non-unique TF/organism combinations
my %vimssRows = ();
my %orgNameRows = ();
foreach my $row (@rows) {
  my $vimss = $row->{id2};
  my $org = $row->{org};
  my $name = $row->{name};
  push @{ $vimssRows{$vimss} }, $row;
  push @{ $orgNameRows{$org}{$name} }, $row;
}

my @filtered = ();
foreach my $row (@rows) {
  my $vimss = $row->{id2};
  my $org = $row->{org};
  my $name = $row->{name};
  push @filtered, $row
    if scalar(@{ $vimssRows{$vimss} }) == 1
      && scalar(@{ $orgNameRows{$org}{$name} }) == 1;
}
print STDERR "Filtered to " . scalar(@filtered). " entries\n";

my $tmpDir = $ENV{TMPDIR} || "/tmp";
die "Not a  directory: $tmpDir\n" unless -d $tmpDir;
my $tmpPre = "$tmpDir/clusterRegPrecise.$$";
my $tmpFaa = "$tmpPre.faa";
my $tmpClust = "$tmpPre.clust";

open (my $fhFaa, ">", $tmpFaa) || die "Cannot write to $tmpFaa\n";
foreach my $row (@filtered) {
  print $fhFaa ">" . $row->{id2} . "\n" . $row->{seq} . "\n";
}
close($fhFaa) || die "Error writing to $tmpFaa\n";

my %vimssRow = map { $_->{id2} => $_ } @filtered;

my $cmd = "$usearch -cluster_fast $tmpFaa -quiet -uc $tmpClust"
  . " -id $identity -query_cov $coverage -target_cov $coverage ";
system($cmd) == 0 || die "usearch failed:\n$cmd\n$!";

my %clusters = (); # clusterId to list
open(my $fhClust, "<", $tmpClust) || die "Cannot read $tmpClust";
while(my $line = <$fhClust>) {
  chomp $line;
  my @F = split /\t/, $line;
  my $type = $F[0];
  my $clusterId = $F[1];
  my $vimss = $F[8];
  die "Unknown vimss $vimss" unless exists $vimssRow{$vimss};
  push @{ $clusters{$clusterId} }, $vimss;
}
close($fhClust) || die "Error reading $tmpClust";

unlink($tmpFaa);
unlink($tmpClust);

print STDERR "Read " . scalar(keys %clusters) . " clusters\n";
my $nOut = 0;
foreach my $clusterId (sort {$a<=>$b} keys %clusters) {
  # Output entries only if they are a new name
  my %namesSeen = ();
  foreach my $vimss (@{ $clusters{$clusterId} }) {
    my $row = $vimssRow{$vimss};
    die $vimss unless defined $row;
    my $name = $row->{name};
    next if exists $namesSeen{$name};
    $namesSeen{$name} = 1;
    print join("\t", map $row->{$_}, qw{db id id2 name desc org seq comment pmIds})."\n";
    $nOut++;
  }
}
print STDERR "Kept $nOut entries\n";
