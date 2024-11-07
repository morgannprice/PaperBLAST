#!/usr/bin/perl -w
use strict;
# Given a PaperBLAST build in the testdb and work directories, build a fasta file of proteins
# that are mentioned in abstracts
use FindBin qw{$RealBin};
use Getopt::Long;

my $work = "work";
my $testdb = "testdb";
my $usage = <<END
Usage: proteinsInAbstracts.faa [ -work $work ] [ -testdb $testdb ] [ -out $testdb/abstracts.faa ]

Makes a fasta file with the sequence of each protein that is mentioned
in an abstract. Relies on the PaperBLAST database
(testdb/litsearch.db), the full set of sequences
(testdb/litsearch.faa), and the pubmed-based identifier list
(work/words.pubmed).
END
;

my $outFile;
die $usage
  unless GetOptions('work=s' => \$work,
                      'testdb=s' => \$testdb,
                      'out=s' => \$outFile)
  && @ARGV == 0;
die "No such directory: $work\n" unless -d $work;
die "No such directory: $testdb\n" unless -d $testdb;
$outFile = "$testdb/abstracts.faa" unless defined $outFile;

my @cmds = (
  qq{sqlite3 --separator '\t' $testdb/litsearch.db 'select geneId,queryTerm,pmId FROM GenePaper WHERE pmId <> "";' > $work/genePaper},
  qq{join.pl -match 1.1=2.3,1.2=2.2 -keep 2.1,2.2,2.3 $work/words.pubmed $work/genePaper > $work/genePaper.pm},
  qq{cut -f 1 $work/genePaper.pm | sort -u > $work/abstract.ids},
  qq{fastaSubset.pl $work/abstract.ids < $testdb/litsearch.faa > $outFile}
);

foreach my $cmd (@cmds) {
  system($cmd) == 0 || die "Command failed\n$cmd\n$!\n";
}
print STDERR "Wrote $outFile\n";
