#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use DBI;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ReadFastaEntry

# List of sequences linked to at least N papers

my $minN = 3;
my $usage = <<END
interesting_seqs.pl [ -n $minN ] -dir data_directory > fasta.out

Builds a fasta file of the protein sequences that are linked to least
n papers. These links can be via curation or via text mining.

The data directory must include the sqlite database (litsearch.db)
and the fasta database (uniq.faa)
END
;

my $dir;
die $usage
  unless GetOptions('-n' => \$minN, 'dir=s' => \$dir)
  && @ARGV == 0
  && defined $dir;
die "Not a directory: $dir\n" unless -d $dir;
my $sqldb = "$dir/litsearch.db";
die "No such file: $sqldb\n" unless -e $sqldb;
my $blastdb = "$dir/litsearch.faa";
die "No such file: $blastdb\n" unless -e $blastdb;
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

# sequence_id to hash of sequence_id, duplicate_id
my $seqToUniq = $dbh->selectall_hashref("SELECT sequence_id, duplicate_id FROM SeqToDuplicate", "duplicate_id");

sub seqToUniq($) {
  my ($seq) = @_;
  return $seq unless exists $seqToUniq->{$seq};
  return $seqToUniq->{$seq}{sequence_id};
}

# maps a unique sequence to papers
my %uniqPapers = (); # uniq => pmId else doi else pmcId => 1

my $rows = $dbh->selectall_arrayref("SELECT geneId, pmId, doi, pmcId FROM GenePaper");
foreach my $row (@$rows) {
  my ($geneId, $pmId, $doi, $pmcId) = @$row;
  my $paperId = $pmId || $doi || $pmcId;
  $uniqPapers{ &seqToUniq($geneId) }{$paperId} = 1;
}

$rows = $dbh->selectall_arrayref("SELECT db, protId, pmId FROM CuratedPaper");
foreach my $row (@$rows) {
  my ($db, $protId, $pmId) = @$row;
  my $geneId = join("::", $db, $protId);
  $uniqPapers{ &seqToUniq($geneId) }{$pmId} = 1;
}

my %ids = ();
foreach my $geneId (sort keys %uniqPapers) {
  $ids{$geneId} = scalar(keys %{ $uniqPapers{$geneId} })
    if scalar(keys %{ $uniqPapers{$geneId} }) >= $minN;
}
print STDERR "Found " . scalar(keys %ids) . " sequences linked to at least $minN papers\n";

# And then fetch those sequences, using ReadFastaEntry
open(my $fh, "<", $blastdb) || die "Cannot read $blastdb";
my $nPrint = 0;
my $state = {};
while(my ($header,$sequence) = ReadFastaEntry($fh,$state)) {
  if (exists $ids{$header}) {
    print ">$header $ids{$header}\n$sequence\n";
    $nPrint++;
  }
}
close($fh) || die "Error reading $blastdb";
die "Not all sequences found: $nPrint printed" unless $nPrint == scalar(keys %ids);
