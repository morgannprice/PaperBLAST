#!/usr/bin/perl -w
# Given the relevant queries (from gapquery.pl),
# execute the searches and make a table of hits (potential candidates)
# Operates on orgs files from buildorgs.pl

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils qw{ReadTable};
use Steps qw{ReadOrgTable ReadOrgProtein ParseOrgLocus};

{
  my $nCPU = $ENV{MC_CORES} || 4;
  my @fields = qw{locusId type curatedId bits locusBegin locusEnd cBegin cEnd cLength identity};
  my $fields = join(" ", @fields);

  my $usage = <<END
Usage: gapsearch.pl -orgs prefix -query queryfile1 ... queryfileN -out hitsfile

The output file is tab-delimited with fields
$fields
where type is either hmm or blast, curatedId is the
query field of the query file or is curated:query,
and identity is missing for hits of type hmm.

Optional arguments:
 -nCPU $nCPU -- number of CPUs to use (defaults to the MC_CORES
                environment variable, or 4)
 -verbose -- set for more verbose output to standard error
END
;

  my @queryfiles;
  my ($hmmdir, $orgprefix, $outfile, $verbose);
  die $usage
    unless GetOptions('query=s{1,}' => \@queryfiles,
                      'orgs=s' => \$orgprefix,
                      'out=s' => \$outfile,
                      'nCPU=i' => \$nCPU,
                      'verbose' => \$verbose)
      && defined $orgprefix
      && @queryfiles > 0
      && defined $outfile;
  die "Must request at least one CPU\n" unless $nCPU >= 1;
  foreach my $file (@queryfiles) {
    die "No such file: $file\n" unless -e $file;
  }
  my $bindir = $Bin;
  my $hmmsearch = "$bindir/hmmsearch";
  my $usearch = "$bindir/usearch";
  foreach my $b ($hmmsearch, $usearch) {
    die "No such file or not executable: $b\n" unless -x $b;
  }
  foreach my $suffix (qw{org faa}) {
    die "No such file: $orgprefix.$suffix\n" unless -e "$orgprefix.$suffix";
  }
  my @orgs = ReadOrgTable("$orgprefix.org");
  die "The organism table $orgprefix.org has no rows\n" unless @orgs > 0;
  my %orgs = map { $_->{orgId} => $_ } @orgs;
  my $aaIn = "$orgprefix.faa";

  my %hmmFile = (); # hmmId => filename
  my %curated = (); # curatedId (or uniprot: or curated2:) to sequence

  my @querycol = qw{step type query desc file sequence};
  foreach my $qfile (@queryfiles) {
    my $qdir;
    if ($qfile =~ m!/!) {
      $qdir = $qfile; $qdir =~ s!/[^/]*$!!;
    } else {
      $qdir = ".";
    }
    my @queries = ReadTable($qfile, \@querycol);
    foreach my $query (@queries) {
      my $id = $query->{query};
      my $type = $query->{type};
      $id = "curated2:" . $id if $type eq "curated2";
      $id = "uniprot:" . $id if $type eq "uniprot";

      if ($type eq "curated" || $type eq "curated2" || $type eq "uniprot") {
        die "Inconsistent sequence for $id in $qfile\n"
          if exists $curated{$id} && $curated{$id} ne $query->{sequence};
        $curated{$id} = $query->{sequence};
      } elsif ($type eq "hmm") {
        my $file = $query->{file};
        die "Invalid characters in hmmfile $file" unless $file =~ m/^[a-zA-Z0-9._-]+$/;
        $file = "$qdir/$file";
        die "No such file: $file\n" unless -e $file;
        die "Inconsistent file for hmm $id in $qfile\n"
          if exists $hmmFile{$id} && $hmmFile{$id} ne $file;
        $hmmFile{$id} = $file;
      } elsif ($type eq "ignore") {
        ;
      } else {
        die "Unrecognized query type $type";
      }
    }
  }
  print STDERR join(" ", "Comparing", scalar(@orgs), "proteomes to",
                    scalar(keys %hmmFile), "HMMs",
                    "and", scalar(keys %curated), "sequences")."\n";

  open (my $fhOut, ">", $outfile) || die "Cannot write to $outfile\n";
  print $fhOut join("\t", @fields)."\n";

  my %hmmTmp = (); # hmm to tmp hits file
  my $nRunning = 0;
  foreach my $hmm (sort keys %hmmFile) {
    $nRunning++;
    if ($nRunning > $nCPU) {
      print STDERR "Waiting\n" if defined $verbose;
      waitpid(-1,0);
      $nRunning--;
      print STDERR "Done waiting\n" if defined $verbose;
    }
    my $tmpfile = "/tmp/gapsearch.$$.$hmm.domtbl";
    $hmmTmp{$hmm} = $tmpfile;
    if (fork() == 0) {
      # the child process
      my @cmd = ($hmmsearch, "--cut_tc", "-o", "/dev/null", "--domtblout", $tmpfile, $hmmFile{$hmm}, $aaIn);
      print STDERR "Running $hmm\n" if defined $verbose;
      system(@cmd) == 0 || die "$@cmd failed: $!";
      print STDERR "$hmm finished\n" if defined $verbose;
      exit(0);
    }
  }
  if ($nRunning > 0) {
    my $kid;
    do { $kid = waitpid(-1, 0); } while $kid > 0;
  }
  print STDERR "All " . scalar(keys %hmmTmp) . " HMM analyses complete\n";

  foreach my $hmm (sort keys %hmmFile) {
    open(my $fhIn, "<", $hmmTmp{$hmm})
      || die "hmmsearch failed to create $hmmTmp{$hmm}\n";
    while (my $line = <$fhIn>) {
      chomp $line;
      next if $line =~ m/^#/;
      my @F = split /\s+/, $line;
      my ($hitId, undef, $hitlen, $hmmName, undef, $hmmLen, $seqeval, $seqscore, undef, undef, undef, undef, $domeval, $domscore, undef, $qbeg, $qend, $hitbeg, $hitend) = @F;
      print $fhOut join("\t", $hitId, "hmm", $hmm, $domscore, $hitbeg, $hitend, $qbeg, $qend, $hmmLen, "")."\n";
    }
    unlink($hmmTmp{$hmm});
  }

  my $cfile = "/tmp/gapsearch.$$.curated";
  open(my $fhC, ">", $cfile) || die "Cannot write to $cfile\n";
  foreach my $curatedId (sort keys %curated) {
    print $fhC join("", ">", $curatedId, "\n", $curated{$curatedId}, "\n");
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
    my ($curatedId, $hitId, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits)
      = split /\t/, $line;
    die "ublast found hit for unknown item $curatedId\n" unless exists $curated{$curatedId};
    $hitId =~ s/ .*//; # the first part of the header line is the locusId
    print $fhOut join("\t", $hitId, "blast", $curatedId, $bits, $sbeg, $send, $qbeg, $qend,
                      length($curated{$curatedId}), $identity)."\n";
  }
  close($fhH) || die "Error reading $cfile.hits\n";
  unlink("$cfile");
  unlink("$cfile.hits");
  close($fhOut) || die "Error writing to $outfile\n";
}
