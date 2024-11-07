#!/usr/bin/perl -w
use strict;

die "Usage: fastaSubset.pl namesfile < fastafile\n" unless @ARGV==1;

my ($namesfile) = @ARGV;

my %names = ();
my %seen = (); # skip duplicates
my $nSkipDup = 0;

open(NAMES,"<",$namesfile) || die "Error reading namesfile $namesfile";
while(my $line = <NAMES>) {
    chomp $line;
    $line =~ m/^ *(\S+)[ \t]*/ || die "Cannot parse $line in namesfile $namesfile";
    my $name = $1;
    die "Duplicate name $name in namesfile $namesfile" if exists $names{$name};
    $names{$name} = 1;
}
close(NAMES) || die "Error reading namesfile";

my $writing = 0;
my $nWritten = 0;
while(my $line = <STDIN>) {
    if (substr($line,0,1) eq ">") {
	$line =~ m/^>(\S+)/ || die "Invalid header line $line in stdin";
	my $name = $1;
	$writing = 0;
	if (exists $names{$name}) {
	    if (exists $seen{$name}) {
		$nSkipDup++;
	    } else {
		$writing = 1;
	    }
	}
	if ($writing) {
	    $nWritten++;
	    $seen{$name} = 1;
	}
    }
    print $line if $writing;
}
print STDERR "Wrote $nWritten entries of " . scalar(keys %names) . " selected names in $namesfile\n";
print STDERR "Skipped $nSkipDup duplicate entries\n" if $nSkipDup > 0;
