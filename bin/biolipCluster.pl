#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};

my $identity = 0.9;
my $minCoverage = 0.8;
my $nThreads = 8;
my $cdhit = "$Bin/cd-hit";
my $usage = <<END
Usage: biolipCluster.pl -in BioLip.txt -curated biolip.curated_parsed > BioLip_clustered.txt

Given the biolip file, uses cd-hit to cluster the sequences, and then
selects representatives so that every sequence in curated is included,
and so that every ligand for that cluster is included. That is, if a
protein binds A and B, there might be one structure with A, and
another structure with B, with a very similar sequence; both will be
retained. Outputs the subset of the input lines for the selected
representatives.

Optional arguments
-cdhit -- defaults to $cdhit
-identity -- defaults to $identity
-minCoverage -- minimum coverage of shorter sequence by cluster representative
  defaults to $minCoverage
-threads -- defaults to $nThreads
END
;

my ($inFile, $curatedFile);
die $usage
  unless GetOptions('in=s' => \$inFile,
                    'curated=s' => \$curatedFile,
                    'cd-hit=s' => \$cdhit,
                    'identity=f' => \$identity,
                    'minCoverage=f' => \$minCoverage,
                    'threads=i' => \$nThreads)
  && @ARGV == 0
  && defined $inFile && defined $curatedFile;
die "No such executable: $cdhit\n" unless -x $cdhit;
die "Identity must be between 0.5 and 1\n" unless $identity >= 0.5 && $identity <= 1.0;
die "Coverage must be between 0 and 1\n" unless $minCoverage >= 0 && $minCoverage <= 1;
$nThreads = 1 if $nThreads < 0;

my $tmpPre = ($ENV{TMPDIR} || "/tmp") . "/biolipCluster.$$";
my $tmpFaa = "$tmpPre.faa";
open(my $fhFaa, ">", $tmpFaa) || die "Cannot write to $tmpFaa";

open(my $fhIn, "<", $inFile) || die "Cannot read $inFile\n";
my %idLigands = (); # pdbIdchain => ligand => 1
while (my $line = <$fhIn>) {
  chomp $line;
  my ($pdbId, $chain, $resolution, $siteId, $ligandId, $ligandChain, $ligandNo,
      $bindingResiduesPDB, $bindingResidues,
      $siteResiduesPDB, $siteResidues,
      $EC, $GO,
      $affinityPM, $affinityMOAD, $affinityPDBBind, $affinityBindingDB,
      $uniprotId, $pubmedIds, $seq) = split /\t/, $line;
  die "Invalid line $line" unless $seq =~ m/^[A-Za-z]+$/;
  next if $seq =~ m/[a-z]/; # ignore any lower-case sequence (rare)
  my $id = $pdbId.$chain;
  if (!exists $idLigands{$id}) {
    print $fhFaa ">" . $id . "\n" . $seq . "\n";
  }
  $idLigands{$id}{$ligandId} = 1;
}
close($fhIn) || die "Error reading $inFile\n";
close($fhFaa) || die "Error writing to $tmpFaa";

my %curated = (); # id => 1 if in the curated-parsed file
open(my $fhCurated, "<", $curatedFile) || die "Cannot read $curatedFile";
my $nUnexpected = 0;
while (my $line = <$fhCurated>) {
  chomp $line;
  my @F = split /\t/, $line;
  my $id = $F[1];
  die "Invalid id $id in $curatedFile" unless $id =~ m/^[0-9a-zA-Z_]+$/;
  $nUnexpected++ if !exists $idLigands{$id};
  $curated{$id} = 1;
}
close($fhCurated) || die "Error reading $curatedFile";
print STDERR "Warning: $nUnexpected ids in $curatedFile are not in $inFile\n"
  if $nUnexpected > 0;

my $logFile = "$tmpPre.log";
my $clusterFile = "$tmpPre.cluster";
# -M 0 -- no memory limit
# -d 5 -- word length 5 (appropriate for higher %identities)
# -d 0 -- stop at 1st space in defline
my $cmd = "$cdhit -i $tmpPre.faa -o $clusterFile -c $identity -M 0 -T $nThreads -n 5 -d 0 -aS $minCoverage >& $logFile";
system($cmd) == 0 || die "cdhit failed\n$cmd\n$!";

open (my $fhClust, "<", "$clusterFile.clstr")
  || die "Cannot read $clusterFile.clstr\nfrom $cmd";
my $iClust;
my @clusters = (); # cluster index => list of members
my @rep = (); # cluster index => representative
while(my $line = <$fhClust>) {
  chomp $line;
  if ($line =~ m/^>Cluster (\d+)$/) {
    $iClust = $1;
    die "Repeated cluster $iClust" if exists $clusters[$iClust];
    $clusters[$iClust] = [];
  } elsif ($line =~ m/^\d+\s+\d+aa,\s+>([0-9A-Za-z._-]+)[.][.][.]\s+(.*)/) {
    my ($id, $atString) = ($1,$2);
    die "Unexpected member line in $clusterFile.clstr" unless defined $iClust;
    die "Unexpected id $id in $clusterFile.clstr" unless exists $idLigands{$id};
    push @{ $clusters[$iClust] }, $id;
    $rep[$iClust] = $id if $atString eq "*";
  } else {
    die "Unexpected line $line";
  }
}
close($fhClust) || die "Error reading $clusterFile.clstr";
print STDERR scalar(keys %idLigands) . " sequences to " . scalar(@clusters) . " clusters\n";

my %keep = ();
foreach my $iClust (0..(scalar(@clusters)-1)) {
  my $clust = $clusters[$iClust];
  my $rep = $rep[$iClust];
  die "No representative for $iClust" unless defined $rep;
  my @keepClust = grep exists $curated{$_}, @$clust;
  push @keepClust, $rep if scalar(@keepClust) == 0;
  my %keepClust = map { $_ => 1 } @keepClust;
  my %ligandsSoFar = ();
  foreach my $id (@keepClust) {
    my @ligands = keys %{ $idLigands{$id} };
    foreach my $ligandId (@ligands) {
      $ligandsSoFar{$ligandId} = 1;
    }
  }
  foreach my $id (@$clust) {
    next if exists $keep{$id};
    my @ligands = keys %{ $idLigands{$id} };
    my @newLigands = grep !exists $ligandsSoFar{$_}, @ligands;
    if (@newLigands > 0) {
      $keepClust{$id} = 1;
      foreach my $ligandId (@ligands) {
        $ligandsSoFar{$ligandId} = 1;
      }
    }
  }
  foreach my $id (keys %keepClust) {
    $keep{$id} = 1;
  }
}

print STDERR "Selected " . scalar(keys %keep) . " representatives to cover all cluster x ligand combinations\n";

unlink($clusterFile);
unlink("$clusterFile.clstr");
unlink($tmpFaa);

# In rare cases, chain identifiers differ only by case, i.e. 1rxs:A and 1rxs:a
# These will cause problems downstreai with BLAST tools.
# So, suppress any ids already seen as lower case
my %lcidSeen = ();
my %keep2 = ();
foreach my $id (sort keys %keep) {
  next if exists $lcidSeen{lc($id)};
  $lcidSeen{lc($id)} = 1;
  $keep2{$id} = 1;
}
print STDERR "Reduced to " . scalar(keys %keep2) . " representatives to ensure that chain identifiers differ by more than just case\n";

open($fhIn, "<", $inFile) || die "Error reading $inFile";
while (my $line = <$fhIn>) {
  chomp $line;
  my ($pdbId, $chain) = split /\t/, $line;
  my $id = $pdbId.$chain;
  print $line."\n" if exists $keep2{$id};
}
