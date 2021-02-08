#!/usr/bin/perl -w
# Build the sqlite3 database of curated proteins for GapMind and Curated Clusters
use strict;
use DBI;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry ReadTable SqliteImport};
sub DoSqliteCmd($$);

my $usage = <<END
buildCurateDb.pl -dir tmp/path.aa

Assumes that the directory already contains curated.faa,
curated.faa.info, curated2.faa, and hetero.tab. If you don't want to
inlude pfam hits (used by curatedClusters.cgi, but not by GapMind
itself), set pfam.hits.tab to be empty.

Optional arguments:
-curated dir/curated.faa
-curated2 dir/curated2.faa
-hetero dir/hetero.tab
-pfamhits dir/pfam.hits.tab
-out dir/curated.db
END
;

my ($dir, $curatedFile, $curated2File, $heteroFile, $pfamHitsFile, $dbFile);
die $usage
  unless GetOptions('dir=s' => \$dir,
                    'curated=s' => \$curatedFile,
                    'curated2=s' => \$curated2File,
                    'hetero=s' => \$heteroFile,
                    'pfamhits=s' => \$pfamHitsFile,
                    'out=s' => \$dbFile)
  && @ARGV == 0
  && defined $dir;

die "No such directory: $dir\n" unless -d $dir;

$curatedFile = "$dir/curated.faa" unless defined $curatedFile;
my $curatedInfoFile = "$curatedFile.info";
$curated2File = "$dir/curated2.faa" unless defined $curated2File;
$heteroFile = "$dir/hetero.tab" unless defined $heteroFile;
$pfamHitsFile = "$dir/pfam.hits.tab" unless defined $pfamHitsFile;
$dbFile = "$dir/curated.db" unless defined $dbFile;

foreach my $file ($curatedFile, $curatedInfoFile, $curated2File, $heteroFile, $pfamHitsFile) {
  die "No such file: $file\n" unless -e $file;
}

my $tmpDir = $ENV{TMP} || "/tmp";
my $tmpDbFile = "$tmpDir/buildCuratedDb.$$.db";
print STDERR "Building temporary database $tmpDbFile\n";

my $schema = "$RealBin/../lib/curated.sql";
system("sqlite3 $tmpDbFile < $schema") == 0
  || die "Error loading schema $schema into $tmpDbFile -- $!";

my @info = ReadTable($curatedInfoFile, ["ids", "length", "descs"]);
my @curatedInfo = map { $_->{descs} =~ s/\r */ /g;
                        [ $_->{ids}, $_->{length}, $_->{descs} ]
                      } @info;
SqliteImport($tmpDbFile, "CuratedInfo", \@curatedInfo);
print STDERR "Loaded CuratedInfo\n";

my @curatedSeq = ();
open(my $fhFaa, "<", $curatedFile) || die "Cannot read $curatedFile\n";
my $state = {};
while (my ($curatedIds, $seq) = ReadFastaEntry($fhFaa, $state)) {
  push @curatedSeq, [ $curatedIds, $seq ];
}
close($fhFaa) || die "Error reading $curatedFile";
SqliteImport($tmpDbFile, "CuratedSeq", \@curatedSeq);
print STDERR "Loaded CuratedSeq\n";

my @curated2 = ();
open($fhFaa, "<", "$curated2File") || die "Cannot read $curated2File\n";
$state = {};
while (my ($header, $seq) = ReadFastaEntry($fhFaa, $state)) {
  my @F = split / /, $header;
  my $protId = shift @F;
  my $desc = join(" ", @F);
  push @curated2, [ $protId, $desc, $seq ];
}
close($fhFaa) || die "Error reading $curated2File";
SqliteImport($tmpDbFile, "Curated2", \@curated2);

my @heteroIn = ReadTable($heteroFile, ["db","protId","comment"]);
my @hetero = map [ $_->{db} . "::" . $_->{protId}, $_->{comment} ], @heteroIn;
SqliteImport($tmpDbFile, "Hetero", \@hetero);
print STDERR "Loaded Hetero\n";

my @pfamHits = ();
open (my $fhHits, "<", $pfamHitsFile) || die "Cannot read $pfamHitsFile\n";
while (my $line = <$fhHits>) {
  chomp $line;
  my ($protId, $hmmName, $hmmAcc, $eval, $bits,
      $protBeg, $protEnd, $protLen,
      $hmmBeg, $hmmEnd, $hmmLen) = split /\t/, $line;
  die unless defined $hmmLen && $hmmLen =~ m/^\d+$/;
  push @pfamHits, [ $protId, $hmmName, $hmmAcc, $eval, $bits,
                    $protBeg, $protEnd, $protLen,
                    $hmmBeg, $hmmEnd, $hmmLen ];
}
SqliteImport($tmpDbFile, "CuratedPFam", \@pfamHits);
print STDERR "Loaded CuratedPFam\n";

system("cp $tmpDbFile $dbFile") == 0 || die "Copying $tmpDbFile to $dbFile failed: $!";
unlink($tmpDbFile);
print STDERR "Built curated database $dbFile\n";



