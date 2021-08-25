#!/usr/bin/perl -w
use strict;
use DBI;

my ($stepsDb, $outDir) = @ARGV;
die "Usage: extractHmms.pl steps.db outdirectory\n"
  unless @ARGV == 2
  && -e $stepsDb && -d $outDir;

my $dbh = DBI->connect("dbi:SQLite:dbname=${stepsDb}","","",{ RaiseError => 1 })
  || die $DBI::errstr;

my $hmmFiles = $dbh->selectall_arrayref(qq{SELECT * from StepQuery WHERE queryType = "hmm";},
                                        { Slice => {} });
my %hmmIdToFile = ();
foreach my $row (@$hmmFiles) {
  $hmmIdToFile{ $row->{hmmId} } = $row->{hmmFileName};
}

my $hmmIds = $dbh->selectcol_arrayref("SELECT hmmId FROM HMM");
foreach my $hmmId (@$hmmIds) {
  my ($hmm) = $dbh->selectrow_array("SELECT hmm FROM HMM WHERE hmmId = ?", {}, $hmmId);
  die "No hmm data for $hmmId in $stepsDb" unless $hmm;
  my $hmmFile = $hmmIdToFile{$hmmId};
  if (!defined $hmmFile) {
    print STDERR "Warning: No query row for $hmmId but it is in the HMM table\n";
    next;
  }
  $hmmFile = "$outDir/$hmmFile";
  open(my $fh, ">", $hmmFile) || die "Cannot write to $hmmFile";
  print $fh $hmm;
  close($fh) || die "Error writing to $hmmFile";
}
print STDERR "Wrote " . scalar(@$hmmIds). " HMMs to $outDir/*.hmm\n";

