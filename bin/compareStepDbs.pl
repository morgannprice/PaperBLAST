#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
use FindBin qw{$RealBin};

my $minCoverage = 0.7;
my $minIdentity = 0.35;
my $usage = <<END
compareStepDbs.pl -in oldDir newDir

Compares two directories that contain curated.db, curated.faa.udb, and
steps.db. This is most useful to see the impact of changing to a new
curated database (i.e., incorporating recently curated items from
SwissProt, BRENDA, or MetaCyc). Most of the report consider changes to
which genes are associated with a step (the StepQuery table), but
changes to curated2 (SwissProt EC numbers without experimental
evidence) are not considered. Also, any changes to the rules are not
considered.

Optional arguments:
-debug -- more status info, and do not remove intermediate files
-minIdentity $minIdentity
-minCoverage $minCoverage
  sequences that are at least this similar to a known
  sequence, and cover at least this fraction of it, 
  are not considered novel
END
;

my $usearch = "$RealBin/usearch";
die "No such executable: $usearch\n" unless -x $usearch;
my $debug; # debug mode if defined

# Given a pair of sequence lists,
# identify sequences that are unique to each list and not very similar
# to a sequence in the other list.
# Returns two references to lists of sequences, unique in 1 and unique in 2
sub CompareSequenceLists($$);

sub CuratedIdsToPFamString($$);

my @in = ();
die $usage
  unless GetOptions('in=s{2,2}' => \@in,
                    'minIdentity=f' => \$minIdentity,
                    'minCoverage=f' => \$minCoverage,
                    'debug' => \$debug)
  && @in == 2;
die "Invalid minIdentity\n" unless $minIdentity >= 0 && $minIdentity < 1;
die "Invalid minCoverage\n" unless $minCoverage >= 0 && $minCoverage < 1;

foreach my $dir (@in) {
  die "No such directory: $dir\n" unless -d $dir;
  foreach my $file ("steps.db", "curated.db", "curated.faa.udb") {
    die "No such file: $dir/$file\n" unless -e "$dir/$file";
  }
}

my ($dir1,$dir2) = @in;
# C1 is a handle for the old curated.db
my $dbhC1 = DBI->connect("dbi:SQLite:dbname=$dir1/curated.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $dbhC2 = DBI->connect("dbi:SQLite:dbname=$dir2/curated.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
# S1 is a handle for the old steps.db
my $dbhS1 = DBI->connect("dbi:SQLite:dbname=$dir1/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $dbhS2 = DBI->connect("dbi:SQLite:dbname=$dir2/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $pathways1 = $dbhS1->selectcol_arrayref("SELECT pathwayId FROM Pathway;");
my $pathways2 = $dbhS2->selectcol_arrayref("SELECT pathwayId FROM Pathway;");
my %pathways1 = map { $_ => 1 } @$pathways1;
my %pathways2 = map { $_ => 1 } @$pathways2;
foreach my $pathwayId (@$pathways1) {
  print STDERR "Pathway in first db only: $pathwayId\n"
    unless exists $pathways2{$pathwayId};
}
foreach my $pathwayId (@$pathways2) {
  print STDERR "Pathway in second db only: $pathwayId\n"
    unless exists $pathways1{$pathwayId};
}

print join("\t", qw{pathwayId stepId only id pfams desc seq})."\n";
foreach my $pathwayId (@$pathways1) {
  next unless exists $pathways2{$pathwayId};
  my $stepObj1 = $dbhS1->selectall_arrayref("SELECT * from Step WHERE pathwayId = ?",
                                            { Slice => {} }, $pathwayId);
  my $stepObj2 = $dbhS2->selectall_arrayref("SELECT * from Step WHERE pathwayId = ?",
                                            { Slice => {} }, $pathwayId);
  my %steps1 = map { $_->{stepId} => 1 } @$stepObj1;
  my %steps2 = map { $_->{stepId} => 1 } @$stepObj2;
  foreach my $stepId (sort keys %steps1) {
    print STDERR "Step in first db only: $pathwayId $stepId\n"
      unless exists $steps2{$stepId};
  }
  foreach my $stepId (sort keys %steps2) {
    print STDERR "Step in second db only: $pathwayId $stepId\n"
      unless exists $steps1{$stepId};
  }
  foreach my $stepId (sort keys %steps1) {
    next unless exists $steps2{$stepId};
    my $queries1 = $dbhS1->selectall_arrayref("SELECT * from StepQuery WHERE pathwayId = ? AND stepId = ?",
                                            { Slice => {} }, $pathwayId, $stepId);
    my $queries2 = $dbhS2->selectall_arrayref("SELECT * from StepQuery WHERE pathwayId = ? AND stepId = ?",
                                            { Slice => {} }, $pathwayId, $stepId);
    # First, compare the curated/predicted/uniprot items
    # Index by sequence in case the id(s) associaed have changed
    # For simplicity, ignore the difference between predicted or high-confidence,
    # and ignore HMM entries
    my %seq1; # sequence to list of rows
    foreach my $queryObj (@$queries1) {
      if ($queryObj->{queryType} eq "curated"
          || $queryObj->{queryType} eq "uniprot"
          || $queryObj->{queryType} eq "predicted") {
        my $seq = $queryObj->{seq};
        die unless $seq;
        push @{ $seq1{$seq} }, $queryObj;
      }
    }
    my %seq2;
    foreach my $queryObj (@$queries2) {
      if ($queryObj->{queryType} eq "curated"
          || $queryObj->{queryType} eq "uniprot"
          || $queryObj->{queryType} eq "predicted") {
        my $seq = $queryObj->{seq};
        die unless $seq;
        push @{ $seq2{$seq} }, $queryObj;
      }
    }
    my @seq1 = sort keys %seq1;
    my @seq2 = sort keys %seq2;
    my ($uniq1, $uniq2) = CompareSequenceLists(\@seq1, \@seq2);
    foreach my $seq (@$uniq1) {
      my $queryObj = $seq1{$seq}[0];
      print join("\t", $pathwayId, $stepId, "1",
                 $queryObj->{curatedIds} || $queryObj->{uniprotId},
                 CuratedIdsToPFamString($dbhC1, $queryObj->{curatedIds} || $queryObj->{uniprotId}),
                 $queryObj->{desc},
                 $seq)."\n";
    }
    foreach my $seq (@$uniq2) {
      my $queryObj = $seq2{$seq}[0];
      print join("\t", $pathwayId, $stepId, "2",
                 $queryObj->{curatedIds} || $queryObj->{uniprotId},
                 CuratedIdsToPFamString($dbhC2, $queryObj->{curatedIds} || $queryObj->{uniprotId}),
                 $queryObj->{desc},
                $seq)."\n";
    }
  }
}

sub WriteFastaList($$) {
  my ($list, $file) = @_;
  die unless @$list > 0;
  open(my $fh, ">", $file) || die "Cannot write to $file\n";
  for (my $i = 0; $i < @$list; $i++) {
    print $fh ">$i\n" . $list->[$i] . "\n";
  }
  close($fh) || die "Error writing $file\n";
}

sub CompareSequenceLists($$) {
  my ($list1, $list2) = @_;
  return [], $list2  if @$list1 == 0;
  return $list1, [] if @$list2 == 0;
  my $tmpPre = "/tmp/compareStepDbs.$$";
  my $faa1 = "$tmpPre.faa1";
  my $faa2 = "$tmpPre.faa2";
  WriteFastaList($list1, $faa1);
  WriteFastaList($list2, $faa2);
  my $hitsFile = "$tmpPre.hits";
  my $usearchCmd = "$usearch -ublast $faa1 -db $faa2 -evalue 1e-3 -blast6out $tmpPre.hits -quiet";
  print "Running: $usearchCmd\n" if defined $debug;
  system($usearchCmd) == 0
    || die "usearch failed: $!\n$usearchCmd\n";
  my @hits;
  open(my $fh, "<", $hitsFile) || die "Cannot read $hitsFile\n";
  while(my $line = <$fh>) {
    chomp $line;
    my @F = split /\t/, $line;
    die "Not enough columns from usearch"
      unless @F >= 12;
    push @hits, \@F;
  }
  close($fh) || die "Error reading $hitsFile\n";

  unless (defined $debug) {
    foreach my $file ($faa1, $faa2, $hitsFile) {
      unlink($file);
    }
  }

  # For each index, set if there is a homolog above the thresholds in the other set
  my %hasSim1;
  my %hasSim2;
  foreach my $row (@hits) {
    my ($i1, $i2, $identity, $alen, $mismatch, $nGaps,
        $qBegin, $qEnd, $sBegin, $sEnd, $evalue, $bits) = @$row;
    # usearch reports identity as a percentage, not a fraction
    next unless $identity / 100.0 >= $minIdentity;
    die "Invalid $i1" unless defined $list1->[$i1];
    die "Invalid $i2" unless defined $list2->[$i2];
    my $len1 = length($list1->[$i1]);
    my $len2 = length($list2->[$i2]);
    my $covOf1 = ($qEnd-$qBegin+1) / $len1;
    my $covOf2 = ($sEnd-$sBegin+1) / $len2;
    $hasSim1{$i1} = 1 if $covOf2 >= $minCoverage;
    $hasSim2{$i2} = 1 if $covOf1 >= $minCoverage;
  }

  # Occasionally usearch will fail to find the self hit for a short sequence,
  # so ensure that identical sequences are not returned as unique
  my %seq1 = map { $list1->[$_] => $_ } (0..(scalar(@$list1)-1));
  my %seq2 = map { $list2->[$_] => $_ } (0..(scalar(@$list2)-1));

  my (@noSim1, @noSim2);
  for (my $i = 0; $i < @$list1; $i++) {
    push @noSim1, $list1->[$i]
      if !exists $hasSim1{$i} && !exists $seq2{ $list1->[$i]};
  }
  for (my $i = 0; $i < @$list2; $i++) {
    push @noSim2, $list2->[$i]
      if !exists $hasSim2{$i} && !exists $seq1{ $list2->[$i] };
  }
  return (\@noSim1, \@noSim2);
}

# Given a curated database handle and a curatedIds, return a string summarizing
# which pfam domains(s) are present. (If the id is not from the curated database,
# returns an empty string.)
sub CuratedIdsToPFamString($$) {
  my ($dbh, $curatedIds) = @_;
  my $pfamHits = $dbh->selectall_arrayref("SELECT * FROM CuratedPFam WHERE curatedIds == ? ORDER by seqFrom",
                                          { Slice => {} }, $curatedIds);
  my @pfamNames = map $_->{hmmName}, @$pfamHits;
  return join(", ", @pfamNames);
}
