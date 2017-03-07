#!/usr/bin/perl -w
# The third phase of the pipeline: running queries against EuropePMC
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my $usage = <<END
run_search.pl -in downloads [ -work work ] [ -test ]

Given queries in work/comb.query (computed by run_terms.pl), run the
queries against EuropePMC.
END
;

my ($indir, $test);
my $workdir = "work";

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
  unless GetOptions('in=s'=> \$indir,
                    'work=s' => \$workdir,
                    'test' => \$test)
  && @ARGV == 0
  && defined $indir && defined $workdir;

print STDERR "Test mode\n" if defined $test;

my $infile = "$workdir/comb.query";
unless (-e $infile) {
  if (defined $test) {
    print STDERR "Warning: no input file $infile\n";
  } else {
    die "No input file $infile\n";
  }
}

&maybe_run("$Bin/queryEuropePMCBatch.pl -in $infile -out $workdir/epmc.part1 -dir $workdir");
if (defined $test || ! -z "$infile.fail") {
  &maybe_run("$Bin/join.pl -header 0 -match 1.1=2.3 -ignore 1.1 $workdir/epmc.part1.fail $infile > $infile.fail");
  &maybe_run("$Bin/queryEuropePMCBatch.pl -in $infile.fail -out $workdir/epmc.part2 -dir $workdir");
} else {
  &maybe_run("echo -n > $infile.fail");
  &maybe_run("echo -n > $workdir/epmc.part2");
}
&maybe_run("cat $workdir/epmc.part1 $workdir/epmc.part2 > $workdir/epmc");
print STDERR "Success\n" if ! defined $test;
