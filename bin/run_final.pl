#!/usr/bin/perl -w
# The last phase of the pipeline: building the database
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;

my @allsteps = qw{sprot db compare};
my $allsteps = join(",",@allsteps);
my $comparedir = "$Bin/../data";

my $usage = <<END
run_final.pl -in downloads -work work [ -test ] [ -steps ... ]
	[ -compare $comparedir ]

Given the downloads and work directories, and the output of
run_snippets.pl, builds a database in work/data, and compares it to
the existing database (from $comparedir)

END
;

my ($test, $indir, $workdir);
my $dosteps = $allsteps;

sub maybe_run($) {
    my ($cmd) = @_;
    if (defined $test) {
        print STDERR "Would run\t$cmd\n";
    } else {
        print STDERR "Running $cmd\n";
        system($cmd) == 0
            || die "Error running $cmd: $!\n";
    }
}

die $usage
    unless GetOptions('in=s' => \$indir,
                      'work=s' => \$workdir,
                      'steps=s' => \$dosteps,
                      'test' => \$test,
                      'compare' => \$comparedir)
    && @ARGV == 0
    && defined $indir && defined $workdir;
die "No such directory: $indir\n" unless -d $indir;
die "No such directory: $workdir\n" unless -d $workdir;
die "No such directory: $comparedir\n" unless -d $comparedir;

my @dosteps = split /,/, $dosteps;
my %dosteps = map { $_ => 1 } @dosteps;
my %allsteps = map { $_ => 1 } @allsteps;
foreach my $step (keys %dosteps) {
    die "Unrecognized step: $step\n" unless exists $allsteps{$step};
}

print STDERR "Test mode\n" if defined $test;

my $outdir = "$workdir/data";
&mkdir_if_needed($outdir);

if (exists $dosteps{sprot}) {
  my $sprot = "$indir/uniprot_sprot.dat.gz";
  die "No such file: $sprot\n" unless -e $sprot;
  &maybe_run("(zcat $sprot | $Bin/sprotCharacterized.pl > $workdir/sprot.char.tab) >& $workdir/sprot.char.log");
}

if (exists $dosteps{db}) {
  die "No such directory: $indir/ecocyc/data\n" unless -d "$indir/ecocyc/data";
  die "No such file: $workdir/snippets_comb\n" unless -e "$workdir/snippets_comb";
  die "No such file: $workdir/hits.papers\n" unless -e "$workdir/hits.papers";
  &maybe_run("$Bin/buildLitDb.pl -dir $outdir -snippets $workdir/snippets_comb -sprot $workdir/sprot.char.tab -ecocyc $indir/ecocyc/data $workdir/hits");
}

if (exists $dosteps{compare}) {
  my $newdb = "$outdir/litsearch.db";
  my $olddb = "$comparedir/litsearch.db";
  die "No such file: $olddb\n" unless -e $olddb;
  &maybe_run("$Bin/compare_dbs.pl -old $olddb -new $newdb -out $workdir/dbcmp");
}
