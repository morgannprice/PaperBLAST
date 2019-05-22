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

# This uses the pre-loaded rows, not a database handle like pbutils::CuratedMatch
# The arguments are a reference to the curated rows (containing desc and ids) and a reference to a list of ids
# Returns a reference to a list of rows
sub WordMatchRows($$);

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
    hmm:PF00001	hmm:TIGR0001
    uniprot:Q8TPT3_METAC	   curated:reanno::Miya:8499492
    ignore:CharProtDB::CH_123581   ignore_other:EC 2.2.1.6

  The step name and step description are arbitrary.
  EC numbers are expanded to matching TIGRFams as well.

  To obtain the list of queries for finding candidates by BLAST, all
  curated items matching the attributes are combined, and then those
  in the ignore list are removed. (The ignore list also affects
  scoring of reverse-hits later on.)

rules lines starting with name: followed by a list of step names or
rule names.

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
-curated curated.faa -- the file of characterized
  sequences. curated.faa.desc must also exist. By default, it looks
  for this in dir/curated.faa
-curated2 curated2.faa -- the file of curated sequences. By default,
  it looks for this in dir/curated.faa
-debug
END
;

  my ($stepsFile, $outDir, $hmmDir, $curatedFaa, $curated2Faa);
  die $usage
    unless GetOptions('steps=s' => \$stepsFile,
                      'outdir=s' => \$outDir,
                      'hmmdir=s' => \$hmmDir,
                      'curated=s' => \$curatedFaa,
                      'curated2=s' => \$curated2Faa,
                      'debug' => \$debug)
      && defined $stepsFile && defined $outDir && defined $hmmDir;
  foreach my $dir ($outDir,$hmmDir) {
    die "No such directory: $dir\n" unless -d $dir;
  }
  $curatedFaa = "$outDir/curated.faa" unless defined $curatedFaa;
  $curated2Faa = "$outDir/curated2.faa" unless defined $curated2Faa;
  @hmmIn = map "$hmmDir/$_", @hmmIn;
  foreach my $file ($stepsFile, @hmmIn, $curatedFaa, "$curatedFaa.info", $curated2Faa) {
    die "No such file: $file\n" unless -e $file;
  }
  my $hmmfetch = "$Bin/hmmfetch";
  foreach my $x ($hmmfetch) {
    die "No such executable: $x\n" unless -x $x;
  }

  my $st = ReadSteps($stepsFile);
  my $steps = $st->{steps};
  my $rules = $st->{rules};

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

  my @curatedInfo = ReadTable("$curatedFaa.info", ["ids","length","descs"]);
  my %curatedInfo = map { $_->{ids} => $_ } @curatedInfo;
  # And add individual ids
  foreach my $row (@curatedInfo) {
    my $ids = $row->{ids};
    my @ids = split /,/, $ids;
    if (@ids > 1) {
      foreach my $id (@ids) {
        $curatedInfo{$id} = $row;
      }
    }
  }

  # Load curated2. For each item, store ids, length, descs, seq
  open(my $fh2, "<", $curated2Faa) || die "Cannot read $curated2Faa\n";
  my $state2 = {};
  my @curated2 = ();
  while (my ($header,$seq) = ReadFastaEntry($fh2, $state2)) {
    my @words = split / /, $header;
    my $id = shift @words;
    die unless defined $id && $id ne "";
    push @curated2, { 'ids' => $id,
                      'length' => length($seq),
                      'descs' => join(" ", @words),
                      'seq' => $seq };
  }
  close($fh2) || die "Error reading $curated2Faa\n";
  my %curated2 = map { $_->{ids} => $_ } @curated2;
  die "Duplicate identifiers in $curated2Faa\n"
    unless scalar(keys %curated2) == scalar(@curated2);

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
  # (@curatedInfo calls them "ids" and they are often comma delimited but here they are treated as a single id)
  my %stepCurated = (); # step => curated id => 1
  my %stepCurated2 = (); # step => curated2 id => 1
  my %stepUniprot = (); # step => uniprot => 1
  my %stepIgnore = (); # step => curated id => 1

  foreach my $step (sort keys %$steps) {
    my $desc = $steps->{$step}{desc};
    my $l = $steps->{$step}{search};
    foreach my $row (@$l) {
      my ($type, $value) = @$row;

      if ($type eq "EC" || $type eq "term") {
        my $curated = WordMatchRows(\@curatedInfo, $value);
        print STDERR "Warning: $step\t$desc: no curated hits for $value\n"
          unless @$curated;
        foreach my $row (@$curated) {
          $stepCurated{$step}{$row->{ids}} = 1;
        }
        my $curated2 = WordMatchRows(\@curated2, $value);
        foreach my $row (@$curated2) {
          $stepCurated2{$step}{$row->{ids}} = 1;
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
        my $ids = $curatedInfo{$value}{ids}
          || die "No curated id matches $value\n";
        # Use the standard combined identifier
        $stepCurated{$step}{$ids} = 1;
      } elsif ($type eq "hmm") {
        $stepHmm{$step}{$value} = 1;
      } elsif ($type eq "uniprot") {
        die "Invalid uniprot identifier $value -- wrong format\n"
          unless $value =~ m/^[0-9A-Z_]+$/;
        $stepUniprot{$step}{$value} = 1;
      } elsif ($type eq "ignore") {
        if (!exists $curatedInfo{$value}) {
          print STDERR "Warning: No curated id matches the ignore $value for step $step\n";
        } else {
          my $ids = $curatedInfo{$value}{ids};
          $stepIgnore{$step}{$ids} = 1;
        }
      } elsif ($type eq "ignore_other") {
        # ignore curated items that match the term unless they have already matched
        my $list = WordMatchRows(\@curatedInfo, $value);
        print STDERR "Warning: no curated hits for $value (used in ignore_other for $step)\n"
          unless @$list > 0;
        foreach my $row (@$list) {
          $stepIgnore{$step}{$row->{ids}} = 1 unless exists $stepCurated{$step}{$row->{ids}};
        }
      } else {
        die "Unknown attribute type $type in step $step\n";
      }
    } # end loop over search terms

    # Ensure that there are some items that match the step
    if (!exists $stepHmm{$step} && !exists $stepUniprot{$step}) {
      die "No curated, curated2, hmm, or uniprot items for $step\n"
        unless exists $stepCurated{$step} || exists $stepCurated2{$step};
      my @ids = grep { !exists $stepIgnore{$step}{$_} } keys %{ $stepCurated{$step} };
      print STDERR "No hmm or uniprot items for $step, and the only curated items are ignored\n"
        unless @ids > 0;
      die "No item can match step $step\n"
        unless @ids > 0 || exists $stepCurated2{$step};
    }
  }

  # Collect descriptions and sequences of uniprot ids
  my %uniprotSeq = ();
  my %uniprotDesc = ();
  while (my ($step, $uniprotHash) = each %stepUniprot) {
    foreach my $uniprotId (keys %$uniprotHash) {
      next if exists $uniprotSeq{$uniprotId};
      my ($seq, $desc) = FetchUniProtSequence($uniprotId);
      die "Cannot fetch uniprot sequence for identifier $uniprotId, is it invalid?\n" unless $seq;
      $uniprotSeq{$uniprotId} = $seq;
      $uniprotDesc{$uniprotId} = $desc;
    }
  }

  # Collect metadata about HMMs and fetch them
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

  # Collect sequences of curated items
  my %curatedFetch = (); # all curated items in stepCurated or stepIgnore
  foreach my $stepHash (\%stepCurated, \%stepIgnore) {
    foreach my $idHash (values %$stepHash) {
      foreach my $id (keys %$idHash) {
        $curatedFetch{$id} = 1;
      }
    }
  }
  my %curatedSeq = (); # curated id to sequence for relevant items
  open (my $fhFaa, "<", $curatedFaa) || die "Cannot read $curatedFaa";
  my $state = {};
  while (my ($header, $sequence) = ReadFastaEntry($fhFaa, $state)) {
    my $id = $header; $id =~ s/ .*//;
    next unless exists $curatedFetch{$id};
    die "Duplicate sequence for $id in $curatedFaa" if exists $curatedSeq{$id};
    die "Non-matching legnths for $id -- $curatedFaa vs. $curatedFaa.info"
      unless length($sequence) == $curatedInfo{$id}{length};
    $curatedSeq{$id} = $sequence;
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
      die unless exists $uniprotSeq{$uniprotId} && exists $uniprotDesc{$uniprotId};
      print $fhO join("\t", $step, "uniprot",
                      $uniprotId, $uniprotDesc{$uniprotId},
                      "", $uniprotSeq{$uniprotId})."\n";
    }
    foreach my $id (sort keys %{ $stepCurated{$step} }) {
      next if exists $stepIgnore{$step}{$id};
      die unless exists $curatedInfo{$id} && exists $curatedSeq{$id};
      print $fhO join("\t", $step, "curated",
                      $id, $curatedInfo{$id}{descs}, "", $curatedSeq{$id})."\n";
    }
    foreach my $id (sort keys %{ $stepCurated2{$step} }) {
      die unless exists $curated2{$id};
      print $fhO join("\t", $step, "curated2",
                      $id, $curated2{$id}{descs}, "", $curated2{$id}{seq})."\n";
    }
    foreach my $id (sort keys %{ $stepIgnore{$step} }) {
      die unless exists $curatedInfo{$id};
      print $fhO join("\t", $step, "ignore",
                      $id, $curatedInfo{$id}{descs}, "", $curatedSeq{$id})."\n";
    }
  }
  close($fhO) || die "Error writing to $outFile\n";
  print STDERR "Wrote $outFile\n";
}

sub WordMatchRows($$) {
  my ($list, $query) = @_;
  die "Searching for empty term"
    unless defined $query && $query ne "";
  my $regexp = quotemeta($query);
  $regexp =~ s/\\%/.*/g;
  my @match = grep { $_->{descs} =~ m/$regexp/i } @$list;
  die "Too many hits in the curated file for $query\n" if @match >= 10000;
  return CuratedWordMatch(\@match, $query, "descs");
}
