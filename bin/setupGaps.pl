#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use Getopt::Long;

my $dataDir = "$Bin/../data";
my @indFiles = qw{uniprot_sprot.dat.gz};
my $curatedParsed;
my $usage = <<END
setupGaps.pl  -ind download_directory -set setName
  Sets up tmp/path.setName

  The download directory must include these files:
@indFiles

Optional arguments:
-data data_directory -- The directory with the
  PaperBLAST database. It defaults to
  $dataDir
-skip 3 -- skip the first 3 commands
-sprot sprot.curated_parsed -- defaults to
  data_directory/../sprot.curated_parsed
-test -- print the commands instead of running them
END
;

my ($indDir, $set, $test);
my $skip = 0;
die $usage
  unless GetOptions('ind=s' => \$indDir,
                    'set=s' => \$set,
                    'data=s' => \$dataDir,
                    'sprot=s' => \$curatedParsed,
                    'skip=i' => \$skip,
                    'test' => \$test)
  && @ARGV == 0;
die $usage unless defined $indDir && defined $set;
my $outDir = "tmp/path.$set";
if (! -d $outDir) {
  print STDERR "Making directory $outDir\n";
  mkdir($outDir) || die "Cannot create $outDir\n";
}
foreach my $dir ($indDir, $dataDir) {
  die "Not a directory: $dir\n" unless -d $dir;
}
foreach my $file (@indFiles) {
  die "No such file: $indDir/$file\n"
    unless -e "$indDir/$file";
}
$curatedParsed = "$dataDir/../sprot.curated_parsed" unless defined $curatedParsed;
die "No such file: $curatedParsed -- use -sprot to change\n" unless -e $curatedParsed;

foreach my $x ("$Bin/usearch", "$Bin/blast/formatdb") {
  die "No such executable: $x\n" unless -x $x;
}

my $tmpPre = "/tmp/setupGaps.$$";
my @cmds = ( ["$Bin/curatedFaa.pl",
              "-db", "$dataDir/litsearch.db",
              "-uniq", "$dataDir/uniq.faa",
              "-out", "$outDir/curated.faa",
              "-curatedids"],
             ["$Bin/blast/formatdb",
              "-p", "T", "-o", "T",
              "-i", "$outDir/curated.faa"],
             ["$Bin/usearch",
              "-makeudb_ublast", "$outDir/curated.faa",
              "-output", "$outDir/curated.faa.udb"],
             ["zcat $indDir/uniprot_sprot.dat.gz | $Bin/sprotCuratedEc.pl > $tmpPre.faa"],
             ["$Bin/usearch",
              "-cluster_fast", "$tmpPre.faa",
              "-id", 0.6,
              "-uc", "$tmpPre.uc"],
             ["$Bin/clusterEc.pl $tmpPre.faa $tmpPre.uc > $outDir/curated2.faa"],
             ["zcat $indDir/uniprot_sprot.dat.gz | $Bin/sprotSubunit.pl -curated $curatedParsed > $tmpPre.sprot.subunits"],
             ["$Bin/findHeteromers.pl -sprot $tmpPre.sprot.subunits > $outDir/hetero.tab"],
             ["$Bin/runPfamHits.pl", $set],
             ["$Bin/buildCuratedDb.pl", "-dir", $outDir] );

die "Invalid skip: $skip\n" if $skip < 0;
while($skip > 0) {
  $skip--;
  shift @cmds;
}
foreach my $cmdlist (@cmds) {
  print STDERR join(" ", @$cmdlist)."\n";
  defined $test || system(@$cmdlist) == 0 || die "Error running @$cmdlist\n$$!n";
}
unlink("$tmpPre.faa");
unlink("$tmpPre.uc");
unlink("$tmpPre.sprot.subunits");
print STDERR "Set up gap files in $outDir\n" unless defined $test;
