#!/usr/bin/perl -w
# Given the relevant queries (from gapquery.pl),
# execute the searches and make a table of hits (potential candidates)
# Operates on orgs files from buildorgs.pl

use strict;
use Getopt::Long;
use DBI;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use pbutils qw{ReadTable};
use Steps qw{ReadOrgTable ReadOrgProtein ParseOrgLocus};

{
  my $nCPU = $ENV{MC_CORES} || 4;
  my @fields = qw{locusId type queryId bits locusBegin locusEnd qBegin qEnd qLength identity};
  my $fields = join(" ", @fields);
  my $maxCand = 4;
  my $maxWeak = 2;

  my $usage = <<END
Usage: gapsearch.pl -orgs prefix -set set -out hitsfile

The output file is tab-delimited with fields
$fields
where type is either hmm or blast,
queryId is curatedIds, uniprot:uniprotId, curated2:protId, or hmmId,
and identity is missing for hits of type hmm.

Optional arguments:
 -pathway pathwayId  (By default, runs all pathways)
 -dir $RealBin/../tmp/path.set
 -nCPU $nCPU -- number of CPUs to use (defaults to the MC_CORES
                environment variable, or 4)
 -verbose -- set for more verbose output to standard error
 -maxCand $maxCand -- how many candidates to keep for each step
    (by bit score)
 -maxWeak $maxWeak -- maximum number of blast-only hits of under 40%
    identity to keep for each step
 -diamond -- use diamond instead of usearch; prefix.faa.dmnd
    should be created first with diamond makedb
END
;

  my ($orgprefix, $set, $outFile, $pathSpec, $inDir, $verbose, $useDiamond);
  die $usage
    unless GetOptions('orgs=s' => \$orgprefix,
                      'set=s' => \$set,
                      'out=s' => \$outFile,
                      'pathway=s' => \$pathSpec,
                      'dir=s' => \$inDir,
                      'nCPU=i' => \$nCPU,
                      'verbose' => \$verbose,
                      'maxCand=i' => \$maxCand,
                      'maxWeak=i' => \$maxWeak,
                      'diamond' => \$useDiamond)
      && defined $orgprefix
      && defined $set
      && defined $outFile;
  die "Must request at least one CPU\n" unless $nCPU >= 1;
  $inDir = "$RealBin/../tmp/path.$set" unless defined $inDir;
  die "No such directory: $inDir\n" unless -d $inDir;
  die "maxCand must be positive\n" unless $maxCand > 0;
  die "maxWeak must be positive\n" unless $maxWeak > 0;

  my $binDir = $RealBin;
  my $hmmsearch = "$binDir/hmmsearch";
  my $usearch = "$binDir/usearch";
  my $diamond = "$binDir/diamond";
  my @exe = ($hmmsearch);
  if (defined $useDiamond) {
    push @exe, $diamond;
  } else {
    push @exe, $usearch;
  }
  foreach my $b (@exe) {
    die "No such file or not executable: $b\n" unless -x $b;
  }
  foreach my $suffix (qw{org faa}) {
    die "No such file: $orgprefix.$suffix\n" unless -e "$orgprefix.$suffix";
  }

  my $dbhS = DBI->connect("dbi:SQLite:dbname=${inDir}/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;

  my @orgs = ReadOrgTable("$orgprefix.org");
  die "The organism table $orgprefix.org has no rows\n" unless @orgs > 0;
  my %orgs = map { $_->{orgId} => $_ } @orgs;
  my $aaIn = "$orgprefix.faa";
  if (defined $useDiamond) {
    die "No such file: $aaIn.dmnd -- use\n  $diamond makedb --in $aaIn -d $aaIn.dmnd\nto create it\n"
      unless -e "$aaIn.dmnd";
  }
  if (defined $pathSpec) {
    my $pathObj = $dbhS->selectall_arrayref("SELECT * FROM Pathway WHERE pathwayId = ?",
                                           {}, $pathSpec);
    die "Pathway $pathSpec is not in ${inDir}/steps.db\n"
      unless @$pathObj == 1;
  }
  my $pathClause = "";
  $pathClause = qq{ AND pathwayId = "$pathSpec" } if defined $pathSpec;
  # Will add the queryId field to each row
  my $queries = $dbhS->selectall_arrayref(qq{ SELECT * FROM StepQuery
                                              WHERE queryType <> "ignore"
                                              $pathClause },
                                          { Slice => {} });
  my %stepQuery = (); # stepId to list of rows
  # The reported queryId is curated2:protId, uniprot:uniprotId, curatedIds or hmmId
  my %querySeq = (); # queryId to row for sequences
  my %queryHmm = (); # queryId to row for hmms
  foreach my $row (@$queries) {
    push @{ $stepQuery{ $row->{stepId} } }, $row;
    my $queryType = $row->{queryType};
    if ($queryType eq "hmm") {
      my $hmmId = $row->{hmmId};
      $queryHmm{$hmmId} = $row;
      $row->{queryId} = $hmmId;
    } elsif ($queryType eq "curated" || $queryType eq "curated2"
             || $queryType eq "uniprot"
             || $queryType eq "predicted") {
      my $id;
      $id = $row->{curatedIds} if $queryType eq "curated";
      $id = "uniprot:" . $row->{uniprotId} if $queryType eq "uniprot";
      # Hits to these are treated as lower confidence by gapsummary.pl
      $id = "curated2:" . $row->{protId} if $queryType eq "curated2";
      $id = "predicted:" . $row->{uniprotId} if $queryType eq "predicted";
      $querySeq{$id} = $row;
      $row->{queryId} = $id;
    } else {
      die "Unknown query type " . $row->{queryType};
    }
  }

  print STDERR join(" ", "Comparing", scalar(@orgs), "proteomes to",
                    scalar(keys %queryHmm), "HMMs",
                    "and", scalar(keys %querySeq), "sequences")."\n";


  my %hmmTmp = (); # hmm to tmp hits file
  my $nRunning = 0;
  foreach my $hmmId (sort keys %queryHmm) {
    $nRunning++;
    if ($nRunning > $nCPU) {
      print STDERR "Waiting\n" if defined $verbose;
      waitpid(-1,0);
      $nRunning--;
      print STDERR "Done waiting\n" if defined $verbose;
    }
    my $row = $queryHmm{$hmmId};
    my $hmmFileName = $inDir . "/" . $row->{hmmFileName};
    my $tmpfile = "/tmp/gapsearch.$$.$hmmId.domtbl";
    $hmmTmp{$hmmId} = $tmpfile;
    if (fork() == 0) {
      # the child process
      my @cmd = ($hmmsearch, "--cut_tc", "-o", "/dev/null", "--domtblout", $tmpfile, $hmmFileName, $aaIn);
      print STDERR "Running $hmmId\n" if defined $verbose;
      system(@cmd) == 0 || die "@cmd failed: $!";
      print STDERR "$hmmId finished\n" if defined $verbose;
      exit(0);
    }
  }
  if ($nRunning > 0) {
    my $kid;
    do { $kid = waitpid(-1, 0); } while $kid > 0;
  }
  print STDERR "All " . scalar(keys %hmmTmp) . " HMM analyses complete\n";

  my %hits = (); # orgId => queryId => list of hits
  # (same fields as output)
  foreach my $hmmId (sort keys %queryHmm) {
    my $nHmm = 0;
    open(my $fhIn, "<", $hmmTmp{$hmmId})
      || die "hmmsearch failed to create $hmmTmp{$hmmId}\n";
    while (my $line = <$fhIn>) {
      chomp $line;
      next if $line =~ m/^#/;
      my @F = split /\s+/, $line;
      my ($hitId, undef, $hitlen, $hmmName, undef, $hmmLen, $seqeval, $seqscore, undef, undef, undef, undef, $domeval, $domscore, undef, $qbeg, $qend, $hitbeg, $hitend) = @F;
      die "Invalid line\n$line\nin $hmmTmp{$hmmId}\n" unless defined $hitend && $hitend =~ m/^\d+$/;
      my $orgId = $hitId; $orgId =~ s/:.*//;
      die "Unknown orgId $orgId in hmm hits from $aaIn" unless exists $orgs{$orgId};
      push @{ $hits{$orgId}{$hmmId} }, [ $hitId, "hmm", $hmmId, $domscore, $hitbeg, $hitend, $qbeg, $qend, $hmmLen, "" ];
      $nHmm++;
    }
    close($fhIn) || die "Error reading $hmmTmp{$hmmId}";
    unlink($hmmTmp{$hmmId});
    print STDERR "Read $nHmm hits for $hmmId\n" if defined $verbose;
  }

  my $cfile = "/tmp/gapsearch.$$.curated";
  open(my $fhC, ">", $cfile) || die "Cannot write to $cfile\n";
  my %queryLen = ();
  foreach my $id (sort keys %querySeq) {
    my $row = $querySeq{$id};
    $queryLen{$id} = length( $row->{seq} );
    print $fhC join("", ">", $id, "\n", $row->{seq}, "\n");
  }
  close($fhC) || die "Error writing to $cfile\n";
  my $cmd;
  if (defined $useDiamond) {
    # By default, diamond returns just 25 alignments, which is not enough if comparing to
    # many genomes.
    $cmd = join(" ",
                $diamond, "blastp", "--query", $cfile, "--db", $aaIn.".dmnd",
                 "--evalue", 0.01 * scalar(@orgs), "--id", 0.3,
                 "--out", "$cfile.hits", "--very-sensitive", "--outfmt", 6,
                "--max-target-seqs", 25 * scalar(@orgs),
                 "--threads", $nCPU);
  } else {
    $cmd = join(" ",
                $usearch, "-ublast", $cfile, "-db", $aaIn,
                "-evalue", 0.01, "-id", 0.3,
                "-blast6out", "$cfile.hits",
                "-threads", $nCPU);
  }
  $cmd .= " > /dev/null 2>&1" unless defined $verbose;
  print STDERR "Running $cmd\n" if defined $verbose;
  system($cmd) == 0 || die "$cmd failed: $!\n";
  print STDERR (defined $useDiamond ? "diamond" : "ublast") . " finished\n";
  open(my $fhH, "<", "$cfile.hits") || die "usearch or diamond did not create $cfile.hits\n";
  while(my $line = <$fhH>) {
    chomp $line;
    my ($queryId, $hitId, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits)
      = split /\t/, $line;
    die "ublast or diamond found hit for unknown item $queryId\n" unless exists $queryLen{$queryId};
    $hitId =~ s/ .*//; # the first part of the header line is the locusId
    my $orgId = $hitId; $orgId =~ s/:.*//;
    die "Unknown org $orgId in hits from $aaIn" unless exists $orgs{$orgId};
    push @{ $hits{$orgId}{$queryId} }, [ $hitId, "blast", $queryId, $bits, $sbeg, $send, $qbeg, $qend,
                            $queryLen{$queryId}, $identity ];
  }
  close($fhH) || die "Error reading $cfile.hits\n";
  unlink("$cfile");
  unlink("$cfile.hits");

  open (my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
  print $fhOut join("\t", @fields)."\n";
  foreach my $orgId (sort keys %hits) {
    my $queryHits = $hits{$orgId};
    my %shortHits = (); # queryId to list of hits to keep

    # Keep just the top/good hits for each query (in this organism)
    foreach my $queryId (keys %$queryHits) {
      my @sorted = sort { $b->[3] <=> $a->[3] } @{ $queryHits->{$queryId} };
      my $old = scalar(@sorted);
      if (@sorted > $maxCand) {
        $#sorted = $maxCand-1;
        print STDERR "$queryId -- from $old to " . scalar(@sorted) . " entries\n"
          if defined $verbose;
      }
      my @keep = ();
      my $nWeak = 0;
      foreach my $row (@sorted) {
        my $identity = $row->[9];
        if ($identity eq "" || $identity >= 40) {
          push @keep, $row; # HMM hit or good hit
        } else {
          $nWeak++;
          push @keep, $row if $nWeak <= $maxWeak;
        }
      }
      print STDERR "$queryId $orgId -- from $old to " . scalar(@sorted) . " entries\n"
        if defined $verbose;
      $shortHits{$queryId} = \@keep;
    }

    # Identify the top candidates for each pathway:step
    # (But, always keep all the hmm hits)
    my %hitsByStep = ();
    foreach my $query (@$queries) {
      my $queryId = $query->{queryId};
      next unless exists $shortHits{$queryId};
      my $stepId = $query->{stepId};
      my $pathwayId = $query->{pathwayId};
      push @{ $hitsByStep{$stepId . ":::" . $pathwayId} }, @{ $shortHits{$queryId} };
    }

    my %locusKeep = ();
    foreach my $hits (values %hitsByStep) {
      my @hits = sort { $b->[3] <=> $a->[3] } @$hits; # by bit score
      my $nThisStep = 0;
      my %locusThisStep;
      foreach my $hit (@hits) {
        my $locusId = $hit->[0];
        if (!exists $locusThisStep{$locusId}) {
          $locusKeep{$locusId} = 1;
          $locusThisStep{$locusId} = 1;
          $nThisStep++;
        }
      }
    }
    print STDERR "Kept for $orgId -- " . scalar(keys %locusKeep) . " loci\n"
      if defined $verbose;

    foreach my $query (sort keys %shortHits) {
      foreach my $row (@{ $shortHits{$query} }) {
        my $locusId = $row->[0];
        print $fhOut join("\t", @$row)."\n"
          if exists $locusKeep{$locusId} || $row->[1] eq "hmm";
      }
    }
  } # end loop over orgId

  close($fhOut) || die "Error writing to $outFile\n";
}
