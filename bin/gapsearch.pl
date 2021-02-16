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
END
;

  my ($orgprefix, $set, $outFile, $pathSpec, $inDir, $verbose);
  die $usage
    unless GetOptions('orgs=s' => \$orgprefix,
                      'set=s' => \$set,
                      'out=s' => \$outFile,
                      'pathway=s' => \$pathSpec,
                      'dir=s' => \$inDir,
                      'nCPU=i' => \$nCPU,
                      'verbose' => \$verbose)
      && defined $orgprefix
      && defined $set
      && defined $outFile;
  die "Must request at least one CPU\n" unless $nCPU >= 1;
  $inDir = "$RealBin/../tmp/path.$set" unless defined $inDir;
  die "No such directory: $inDir\n" unless -d $inDir;

  my $binDir = $RealBin;
  my $hmmsearch = "$binDir/hmmsearch";
  my $usearch = "$binDir/usearch";
  foreach my $b ($hmmsearch, $usearch) {
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

  if (defined $pathSpec) {
    my $pathObj = $dbhS->selectall_arrayref("SELECT * FROM Pathway WHERE pathwayId = ?",
                                           {}, $pathSpec);
    die "Pathway $pathSpec is not in ${inDir}/steps.db\n"
      unless @$pathObj == 1;
  }
  my $pathClause = "";
  $pathClause = qq{ AND pathwayId = "$pathSpec" } if defined $pathSpec;
  my $querySeq = $dbhS->selectall_arrayref(qq{ SELECT DISTINCT queryType, curatedIds, uniprotId, protId, seq FROM StepQuery
                                               WHERE queryType IN ("curated", "curated2", "uniprot")
                                               $pathClause
                                               ORDER BY queryType, curatedIds, uniprotId, protId },
                                           { Slice => {} });
  my $queryHmm = $dbhS->selectall_arrayref(qq{ SELECT DISTINCT hmmId, hmmFileName FROM StepQuery
                                               WHERE queryType IN ("hmm")
                                               $pathClause
                                               ORDER BY hmmId },
                                          { Slice => {} });

  print STDERR join(" ", "Comparing", scalar(@orgs), "proteomes to",
                    scalar(@$queryHmm), "HMMs",
                    "and", scalar(@$querySeq), "sequences")."\n";

  open (my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
  print $fhOut join("\t", @fields)."\n";

  my %hmmTmp = (); # hmm to tmp hits file
  my $nRunning = 0;
  foreach my $row (@$queryHmm) {
    $nRunning++;
    if ($nRunning > $nCPU) {
      print STDERR "Waiting\n" if defined $verbose;
      waitpid(-1,0);
      $nRunning--;
      print STDERR "Done waiting\n" if defined $verbose;
    }
    my $hmmId = $row->{hmmId};
    my $hmmFileName = $inDir . "/" . $row->{hmmFileName};
    my $tmpfile = "/tmp/gapsearch.$$.$hmmId.domtbl";
    $hmmTmp{$hmmId} = $tmpfile;
    if (fork() == 0) {
      # the child process
      my @cmd = ($hmmsearch, "--cut_tc", "-o", "/dev/null", "--domtblout", $tmpfile, $hmmFileName, $aaIn);
      print STDERR "Running $hmmId\n" if defined $verbose;
      system(@cmd) == 0 || die "$@cmd failed: $!";
      print STDERR "$hmmId finished\n" if defined $verbose;
      exit(0);
    }
  }
  if ($nRunning > 0) {
    my $kid;
    do { $kid = waitpid(-1, 0); } while $kid > 0;
  }
  print STDERR "All " . scalar(keys %hmmTmp) . " HMM analyses complete\n";

  foreach my $row (@$queryHmm) {
    my $hmmId = $row->{hmmId};
    open(my $fhIn, "<", $hmmTmp{$hmmId})
      || die "hmmsearch failed to create $hmmTmp{$hmmId}\n";
    while (my $line = <$fhIn>) {
      chomp $line;
      next if $line =~ m/^#/;
      my @F = split /\s+/, $line;
      my ($hitId, undef, $hitlen, $hmmName, undef, $hmmLen, $seqeval, $seqscore, undef, undef, undef, undef, $domeval, $domscore, undef, $qbeg, $qend, $hitbeg, $hitend) = @F;
      print $fhOut join("\t", $hitId, "hmm", $hmmId, $domscore, $hitbeg, $hitend, $qbeg, $qend, $hmmLen, "")."\n";
    }
    unlink($hmmTmp{$hmmId});
  }

  my $cfile = "/tmp/gapsearch.$$.curated";
  my %queryLen = ();
  open(my $fhC, ">", $cfile) || die "Cannot write to $cfile\n";
  foreach my $row (@$querySeq) {
    my $id;
    if ($row->{queryType} eq "curated") {
      $id = $row->{curatedIds};
    } elsif ($row->{queryType} eq "curated2") {
      $id = "curated2:" . $row->{protId}
    } elsif ($row->{queryType} eq "uniprot") {
      $id = "uniprot:" . $row->{uniprotId};
    } else {
      die "Cannot handle query type $row->{queryType}";
    }
    $queryLen{$id} = length( $row->{seq} );
    print $fhC join("", ">", $id, "\n", $row->{seq}, "\n");
  }
  close($fhC) || die "Error writing to $cfile\n";
  my $cmd = join(" ",
                 $usearch, "-ublast", $cfile, "-db", $aaIn,
                 "-evalue", 0.01, "-id", 0.3,
                 "-blast6out", "$cfile.hits",
                 "-threads", $nCPU);
  $cmd .= " >& /dev/null" unless defined $verbose;
  print STDERR "Running $cmd\n" if defined $verbose;
  system($cmd) == 0 || die "$cmd failed: $!\n";
  print STDERR "ublast finished\n";
  open(my $fhH, "<", "$cfile.hits") || die "usearch did not create $cfile.hits\n";
  while(my $line = <$fhH>) {
    chomp $line;
    my ($id, $hitId, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits)
      = split /\t/, $line;
    die "ublast found hit for unknown item $id\n" unless exists $queryLen{$id};
    $hitId =~ s/ .*//; # the first part of the header line is the locusId
    print $fhOut join("\t", $hitId, "blast", $id, $bits, $sbeg, $send, $qbeg, $qend,
                      $queryLen{$id}, $identity)."\n";
  }
  close($fhH) || die "Error reading $cfile.hits\n";
  unlink("$cfile");
  unlink("$cfile.hits");
  close($fhOut) || die "Error writing to $outFile\n";
}
