#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $indir = ".";
my $usage = <<END
Usage: buildSiteTables.pl -db sites.db

Optional arguments:
  -indir $indir -- the directory with the input
    tab-delimited files
    (Sites.tab, HasSites.tab, PDBLigands.tab)

The schema should already have been set up (using litsearch.sql)
END
;

my $dbfile;
die $usage
  unless GetOptions('db=s' => \$dbfile,
                    'indir=s' => \$indir);
die "Not a directory: $indir\n" unless -d $indir;
die "Must specify -db:\n$usage" unless defined $dbfile && $dbfile ne "";
die "database file $dbfile must already exist (with schema loaded)\n"
  unless -s $dbfile;

open(my $fhSQL, "|-", "sqlite3", $dbfile)
  || die "Cannot run sqlite3 on $dbfile\n";
print $fhSQL <<END
.mode tabs
.import $indir/Sites.tab Site
.import $indir/HasSites.tab HasSites
.import $indir/PDBLigands.tab PDBLigand
END
;
close($fhSQL) || die "Error running loading commands for sqlite3 on $dbfile\n";
print STDERR "Loaded site tables into $dbfile from ${indir}/\n"
