#!/usr/bin/perl -w
# Find the relevant queries for a set of gaps
# Operates on the PaperBLAST database

use strict;
use Getopt::Long;
use DBI;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Steps qw{ReadSteps FetchUniProtSequence};
use pbutils qw{ReadTable ReadFastaEntry CuratedWordMatch};

# The arguments are a handle to the curated database, and the query term
# They perform word-based matching.
# They return a reference to a list of rows (hashes), either from CuratedInfo or Curated2
sub FindCuratedMatching($$);
sub FindCurated2Matching($$);

# database handle and curatedId to curatedIds, or undef if unknown
sub CuratedIdToIds($$);

# database handle and curatedIds to a row from CuratedInfo, or undef if unkown
sub CuratedIdsToInfo($$);

my $debug;

{
  my @hmmIn = qw{tigrinfo pfam.tab TIGRFAMs.hmm Pfam-A.hmm};
  my $usage = <<END
Usage: gapquery.pl -hmmdir hmmdir -steps stepsfile -outdir dir

stepsfile includes the following lines:
comment lines (that start with #)

step lines, which are tab-delimited with fields step name, step description,
  and attributes of the form type:value. Examples of attributes:
    EC:1.1.1.1	                   term:search string
    hmm:PF00001	hmm:TIGR0001       ignore_hmm:TIGR00519
    uniprot:Q8TPT3_METAC	   curated:reanno::Miya:8499492
    ignore:CharProtDB::CH_123581   ignore_other:EC 2.2.1.6

  The step name and step description are arbitrary.
  EC numbers are expanded to matching TIGRFams as well.

  To obtain the list of queries for finding candidates by BLAST, all
  curated items matching the attributes are combined, and then those
  in the ignore list are removed. (The ignore list also affects
  scoring of reverse-hits later on.)

rules lines starting with name: followed by a (space-delimited) list
of step names or rule names.

import lines of the form
import otherpathway.steps:rule-or-step1,rule-or-step2,...rule-or-stepN
  (Recursive imports are not supported -- if an imported rule requires steps
   from yet another steps file, import those steps first.)

blank (whitespace only lines)

This script considers only the step lines.

Writes to dir/steps.query, a tab-delimited file with the fields
step, type, query, desc, file (optional), sequence (optional).
  type is one of: curated (which usually means characterized), curated2
    (just curated), hmm, uniprot, or ignore.
  query is the sequence identifier(s), HMM identifier, or the uniprot
    identifier (If multiple curated items have the same sequence, all
    identifiers are joined by ",") sequence is not set for hmm items;
    those have file set instead (i.e., PF02965.17.hmm)
Also creates the HMM files in dir.

The hmm directory must contain these files:
@hmmIn

Optional arguments:
-curated curated.faa -- the file of characterized sequences. By
  default, it looks for this in dir/curated.faa
-curatedDb curated.db -- the file of characterized
  sequences. By default, it looks for this in dir/curated.db
-uniprot dir/uniprot.tsv -- tab-delimited file with cache of uniprot
  sequences.
-debug
END
;

  my ($stepsFile, $outDir, $hmmDir, $curatedFaa, $curatedDb, $uniprotFile);
  die $usage
    unless GetOptions('steps=s' => \$stepsFile,
                      'outdir=s' => \$outDir,
                      'hmmdir=s' => \$hmmDir,
                      'curated=s' => \$curatedFaa,
                      'curatedDb=s' => \$curatedDb,
                      'uniprot=s' => \$uniprotFile,
                      'debug' => \$debug)
      && defined $stepsFile && defined $outDir && defined $hmmDir;
  foreach my $dir ($outDir,$hmmDir) {
    die "No such directory: $dir\n" unless -d $dir;
  }
  $curatedFaa = "$outDir/curated.faa" unless defined $curatedFaa;
  $curatedDb = "$outDir/curated.db" unless defined $curatedDb;
  $uniprotFile = "$outDir/uniprot.tsv" unless defined $uniprotFile;

  @hmmIn = map "$hmmDir/$_", @hmmIn;
  foreach my $file ($stepsFile, @hmmIn, $curatedFaa, $curatedDb) {
    die "No such file: $file\n" unless -e $file;
  }
  my $hmmfetch = "$Bin/hmmfetch";
  foreach my $x ($hmmfetch) {
    die "No such executable: $x\n" unless -x $x;
  }

  print STDERR "Reading $stepsFile\n" if $debug;
  my $st = ReadSteps($stepsFile);
  my $steps = $st->{steps};
  my $rules = $st->{rules};

  print STDERR "Reading HMM info\n" if $debug;
  my @tigrinfo = ReadTable("$hmmDir/tigrinfo", ["tigrId","ec","definition","type"]);
  my %tigrinfo = map { $_->{tigrId} => $_ } @tigrinfo;
  my %ecTIGR = ();
  foreach my $row (@tigrinfo) {
    my $ec = $row->{ec};
    push @{ $ecTIGR{$ec} }, $row if $ec ne "";
  }

  my @pfam = ReadTable("$hmmDir/pfam.tab", ["acc","name"]);
  my %pfam = map { $_->{acc} => $_ } @pfam;
  my %pfToAcc = ();
  foreach my $acc (keys %pfam) {
    $acc =~ m/^(PF\d+)[.]\d+$/ || die "Invalid pfam accession $acc in $hmmDir/pfam.tab";
    my $pf = $1;
    $pfToAcc{$pf} = $acc;
  }

  my $dbhC = DBI->connect("dbi:SQLite:dbname=$curatedDb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my @uniprot;
  if (-e $uniprotFile) {
    print STDERR "Reading uniprot cache from $uniprotFile (remove it if broken)\n"
      if $debug;
    @uniprot = ReadTable($uniprotFile, ["uniprotId", "desc", "seq"])
  }
  my %uniprot = map { $_->{uniprotId} => $_ } @uniprot;

  # Check that each step or rule is used as a dependency, except for the rule named all
  my %dep = ();
  foreach my $lists (values %$rules) {
    foreach my $list (@$lists) {
      foreach my $member (@$list) {
        $dep{$member} = 1;
      }
    }
  }
  foreach my $step (keys %$steps) {
    print STDERR "Warning: step $step is not used by any rule\n"
      unless exists $dep{$step};
  }
  foreach my $rule (keys %$rules) {
    print STDERR "Warning: rule $rule is not used by any other rule\n"
      unless exists $dep{$rule} || $rule eq "all";
  }

  my %stepHmm = (); # step => hmm => 1
  # Ignore HMMs that might be matched by EC terms
  my %stepIgnoreHmm = (); # step => hmm => 1
  my %stepCurated = (); # step => curated id => info
  my %stepCurated2 = (); # step => curated2 id => info
  my %stepUniprot = (); # step => uniprot => 1
  my %stepIgnore = (); # step => curated id => info

  print STDERR "Fetching matching items\n" if $debug;
  foreach my $step (sort keys %$steps) {
    my $desc = $steps->{$step}{desc};
    my $l = $steps->{$step}{search};
    foreach my $row (@$l) {
      my ($type, $value) = @$row;

      if ($type eq "EC" || $type eq "term") {
        my $curated = FindCuratedMatching($dbhC, $value);
        print STDERR "Warning: $step\t$desc: no curated hits for $value\n"
          unless @$curated;
        foreach my $row (@$curated) {
          $stepCurated{$step}{$row->{curatedIds}} = $row;
        }
        my $curated2 = FindCurated2Matching($dbhC, $value);
        foreach my $row (@$curated2) {
          $stepCurated2{$step}{$row->{protId}} = $row;
        }
        print STDERR "Warning: $step\t$desc: no curated hits or TIGRFam for EC:$value\n"
          if $type eq "EC" && @$curated == 0 && !exists $ecTIGR{$value};
      }

      if ($type eq "EC") {
        if (exists $ecTIGR{$value}) {
          foreach my $row (@{ $ecTIGR{$value} }) {
            $stepHmm{$step}{ $row->{tigrId} } = 1;
          }
        }
      } elsif ($type eq "term") {
        ;
      } elsif ($type eq "curated") {
        my $curatedIds = CuratedIdToIds($dbhC, $value);
        die "Unknown curated id $value" unless defined $curatedIds;
        my $info = CuratedIdsToInfo($dbhC, $curatedIds);
        $stepCurated{$step}{$curatedIds} = $info;
        # Warn if the item has already been ignored (it will be skipped later on)
        if (exists $stepIgnore{$step}{$curatedIds}) {
          print STDERR "Warning: for step $step, $curatedIds was ignored and then added -- it is still ignored\n";
        }
      } elsif ($type eq "hmm") {
        $stepHmm{$step}{$value} = 1;
      } elsif ($type eq "uniprot") {
        die "Invalid uniprot identifier $value -- wrong format\n"
          unless $value =~ m/^[0-9A-Z_]+$/;
        $stepUniprot{$step}{$value} = 1;
      } elsif ($type eq "ignore") {
        my $curatedIds = CuratedIdToIds($dbhC, $value);
        if (!defined $curatedIds) {
          print STDERR "Warning: No curated id matches the ignore $value for step $step\n";
        } else {
          my $info = CuratedIdsToInfo($dbhC, $curatedIds);
          $stepIgnore{$step}{$curatedIds} = $info;
        }
      } elsif ($type eq "ignore_other") {
        # ignore curated items that match the term unless they have already matched
        my $list = FindCuratedMatching($dbhC, $value);
        print STDERR "Warning: no curated hits for $value (used in ignore_other for $step)\n"
          unless @$list > 0;
        foreach my $row (@$list) {
          $stepIgnore{$step}{$row->{curatedIds}} = $row
            unless exists $stepCurated{$step}{$row->{curatedIds}};
        }
      } elsif ($type eq "ignore_hmm") {
        $stepIgnoreHmm{$step}{$value} = 1;
      } else {
        die "Unknown attribute type $type in step $step\n";
      }
    } # end loop over search terms

    # Remove stepIgnoreHmm items from stepHmm
    foreach my $hmm (keys %{ $stepIgnoreHmm{$step} }) {
      delete $stepHmm{$step}{$hmm};
    }

    # Ensure that there are some items that match the step
    if (keys(%{ $stepHmm{$step} }) == 0 && !exists $stepUniprot{$step}) {
      die "No curated, curated2, hmm, or uniprot items for $step\n"
        unless exists $stepCurated{$step} || exists $stepCurated2{$step};
      my @ids = grep { !exists $stepIgnore{$step}{$_} } keys %{ $stepCurated{$step} };
      print STDERR "No hmm or uniprot items for $step, and the only curated items are ignored\n"
        unless @ids > 0;
      die "No item can match step $step\n"
        unless @ids > 0 || exists $stepCurated2{$step};
    }
  }

  print STDERR "Fetching UniProt sequences (not already cached)\n" if $debug;
  # Collect descriptions and sequences of uniprot ids not already known
  while (my ($step, $uniprotHash) = each %stepUniprot) {
    foreach my $uniprotId (keys %$uniprotHash) {
      next if exists $uniprot{$uniprotId};
      my ($seq, $desc) = FetchUniProtSequence($uniprotId);
      die "Cannot fetch uniprot sequence for identifier $uniprotId, is it invalid?\n" unless $seq;
      $uniprot{$uniprotId} = { 'uniprotId' => $uniprotId,
                               'desc' => $desc,
                               'seq' => $seq };
    }
  }

  print STDERR "Rewriting the UniProt cache\n" if $debug;
  open (my $fhCache, ">", $uniprotFile)
    || die "Cannot write to $uniprotFile\n";
  print $fhCache join("\t", qw{uniprotId desc seq})."\n";
  foreach my $uniprotId (sort keys %uniprot) {
    print $fhCache join("\t", $uniprotId, $uniprot{$uniprotId}{desc}, $uniprot{$uniprotId}{seq})."\n";
  }
  close($fhCache) || die "Error writing $uniprotFile";

  print STDERR "Collect metadata about HMMs and fetch them\n" if $debug;
  my %hmmInfo = (); # hmmId => hash of file, desc
  while (my ($step, $hmmHash) = each %stepHmm) {
    foreach my $hmmId (keys %$hmmHash) {
      next if exists $hmmInfo{$hmmId};
      my ($hmmDb, $acc, $desc);

      if ($hmmId =~ m/^PF\d+$/) {
        die "Unknown PFam $hmmId\n" unless exists $pfToAcc{$hmmId};
        $acc = $pfToAcc{$hmmId};
        $desc = $pfam{$pfToAcc{$hmmId}}{'name'} || die;
        $hmmDb = "$hmmDir/Pfam-A.hmm";
      } elsif ($hmmId =~ m/^TIGR\d+$/) {
        die "Unknown HMM $hmmId\n"
          unless exists $tigrinfo{$hmmId};
        my $ti = $tigrinfo{$hmmId};
        print STDERR "Warning: $hmmId ($ti->{definition}) is of type $ti->{type}\n"
            unless $ti->{type} eq "equivalog" || $ti->{type} eq "equivalog_domain";
        $acc = $hmmId;
        $hmmDb = "$hmmDir/TIGRFAMs.hmm";
        $desc = $ti->{definition};
        $desc = $ti->{geneSymbol} . ": " . $desc if $ti->{geneSymbol};
        $desc .= " (EC $ti->{ec})" if $ti->{ec};
      } else {
        die "Unknown HMM $hmmId\n";
      }
      my $hmmFile = "$outDir/$acc.hmm";
      if (-e $hmmFile) {
        print STDERR "Using existing file $hmmFile\n" if $debug;
      } else {
        my $cmd = "$hmmfetch $hmmDb $acc > $hmmFile.tmp";
        system($cmd)==0 || die "Failed to fetch hmm $acc from $hmmDb: $!";
        rename("$hmmFile.tmp", $hmmFile) || die "Failed to rename to $hmmFile";
        print STDERR "Created $hmmFile\n" if $debug;
      }
      $hmmInfo{$hmmId} = { 'file' => "$acc.hmm", 'desc' => $desc };
    }
  }

  print STDERR "Fetching curated sequences\n" if $debug;
  my %curatedFetch = (); # all curated items in stepCurated or stepIgnore
  foreach my $stepHash (\%stepCurated, \%stepIgnore) {
    foreach my $idHash (values %$stepHash) {
      foreach my $id (keys %$idHash) {
        $curatedFetch{$id} = 1;
      }
    }
  }

  my %curatedSeq = (); # curated id to sequence for relevant items
  foreach my $curatedIds (keys %curatedFetch) {
    my ($seq) = $dbhC->selectrow_array("SELECT seq FROM CuratedSeq WHERE curatedIds = ?",
                                      {}, $curatedIds);
    die "No sequence for $curatedIds" unless defined $seq && $seq ne "";
    $curatedSeq{$curatedIds} = $seq;
  }

  my $outFile = $stepsFile;
  $outFile =~ s!^.*[/]!!;
  $outFile =~ s/[.].*$//;
  $outFile = "$outDir/$outFile.query";
  open(my $fhO, ">", $outFile) || die "Cannot write to $outFile\n";

  print $fhO join("\t", qw{step type query desc file sequence})."\n";
  foreach my $step (sort keys %$steps) {
    foreach my $hmmId (sort keys %{ $stepHmm{$step} }) {
      die unless exists $hmmInfo{$hmmId};
      print $fhO join("\t", $step, "hmm", $hmmId, $hmmInfo{$hmmId}{desc}, $hmmInfo{$hmmId}{file}, "")."\n";
    }
    foreach my $uniprotId (sort keys %{ $stepUniprot{$step} }) {
      die unless exists $uniprot{$uniprotId};
      print $fhO join("\t", $step, "uniprot",
                      $uniprotId, $uniprot{$uniprotId}{desc},
                      "", $uniprot{$uniprotId}{seq})."\n";
    }
    foreach my $id (sort keys %{ $stepCurated{$step} }) {
      next if exists $stepIgnore{$step}{$id};
      my $info = $stepCurated{$step}{$id};
      die unless exists $curatedSeq{$id};
      print $fhO join("\t", $step, "curated",
                      $id, $info->{descs}, "", $curatedSeq{$id})."\n";
    }
    foreach my $id (sort keys %{ $stepCurated2{$step} }) {
      my $info = $stepCurated2{$step}{$id};
      print $fhO join("\t", $step, "curated2",
                      $id, $info->{desc}, "", $info->{seq})."\n";
    }
    foreach my $id (sort keys %{ $stepIgnore{$step} }) {
      my $info = $stepIgnore{$step}{$id};
      die unless exists $curatedSeq{$id};
      print $fhO join("\t", $step, "ignore",
                      $id, $info->{descs}, "", $curatedSeq{$id})."\n";
    }
  }
  close($fhO) || die "Error writing to $outFile\n";
  print STDERR "Wrote $outFile\n";
}

sub FindCuratedMatching($$) {
  my ($dbhC, $query) = @_;
  die "Searching for empty term"
    unless defined $query && $query ne "";
  if ($query =~ m/^[0-9][.][0-9-]+[.][0-9-]+[.][A-Za-z]?[0-9-]*$/) {
    return $dbhC->selectall_arrayref(qq{ SELECT * from ECToCurated
                                         JOIN CuratedInfo USING (curatedIds)
                                         WHERE ec = ? },
                                     { Slice => {} }, $query);
  }
  # else
  my $rows = $dbhC->selectall_arrayref("SELECT * from CuratedInfo WHERE descs LIKE ?",
                                       { Slice => {} }, "%" . $query . "%");
  return CuratedWordMatch($rows, $query, "descs");
}

sub FindCurated2Matching($$) {
  my ($dbhC, $query) = @_;
  die "Searching for empty term"
    unless defined $query && $query ne "";
  if ($query =~ m/^[0-9][.][0-9-]+[.][0-9-]+[.][A-Za-z]?[0-9-]*$/) {
    return $dbhC->selectall_arrayref(qq{ SELECT * from ECToCurated2
                                         JOIN Curated2 USING (protId)
                                         WHERE ec = ? },
                                     { Slice => {} }, $query);
  }
  #else
  my $rows = $dbhC->selectall_arrayref("SELECT * from Curated2 WHERE desc LIKE ?",
                                       { Slice => {} }, "%" . $query . "%");
  return CuratedWordMatch($rows, $query, "desc");
}

sub CuratedIdToIds($$) {
  my ($dbhC, $curatedId)  = @_;
  my ($curatedIds) = $dbhC->selectrow_array("SELECT curatedIds FROM curatedIdToIds WHERE curatedId = ?",
                                            {}, $curatedId);
  return ($curatedIds);
}

sub CuratedIdsToInfo($$) {
  my ($dbhC, $curatedIds) = @_;
  return $dbhC->selectrow_hashref("SELECT * from CuratedInfo WHERE curatedIds = ?",
                                  { Slice => {} }, $curatedIds);
}
