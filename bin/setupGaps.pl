#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use Getopt::Long;

my $dataDir = "$Bin/../data";
my @indFiles = qw{uniprot_sprot.dat.gz};
my $usage = <<END
setupGaps.pl  -ind download_directory -out output_directory

The download directory must include these files:
@indFiles

Optional arguments:
-data data_directory -- The directory with the
  PaperBLAST database. It defaults to
  $dataDir
END
;

my ($indDir, $outDir);
die $usage
  unless GetOptions('ind=s' => \$indDir,
                    'out=s' => \$outDir,
                    'data=s' => \$dataDir)
  && @ARGV == 0;
die $usage unless defined $indDir && defined $outDir;
foreach my $dir ($indDir, $outDir, $dataDir) {
  die "Not a directory: $dir\n" unless -d $dir;
}
foreach my $file (@indFiles) {
  die "No such file: $indDir/$file\n"
    unless -e "$indDir/$file";
}

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
             ["$Bin/clusterEc.pl $tmpPre.faa $tmpPre.uc > $outDir/curated2.faa"]
           );
foreach my $cmdlist (@cmds) {
  print STDERR join(" ", @$cmdlist)."\n";
  system(@$cmdlist) == 0 || die "Error running @$cmdlist\n$$!n";
}
unlink("$tmpPre.faa");
unlink("$tmpPre.uc");
print STDERR "Set up gap files in $outDir\n";
